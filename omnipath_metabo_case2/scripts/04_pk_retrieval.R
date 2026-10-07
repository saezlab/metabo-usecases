# User Story 1 (spec 002-case-study-2-network, FR-001-005): retrieve PK for
# every Morita et al. measured feature from the local omnipath-metabo PKN
# snapshot, and report node/edge counts against Morita et al. Figure 2.
#
# Run from omnipath_metabo_case2/, with Spatial-COSMOS-MISTy as a sibling
# checkout (quickstart.md).

suppressMessages(pkgload::load_all("../../Spatial-COSMOS-MISTy"))
source("scripts/lib/pk_helpers.R")

PKN_DIR <- "../../Spatial-COSMOS-MISTy/data/PKN"

## ---------------------------------------------------------------------
## 4.1 Load the local PKN snapshot, one category at a time (T009)
## ---------------------------------------------------------------------

# case-study-2's own category vocabulary (data-model.md) vs. the package's
# (cosmos_pkn_categories uses "enzyme_metabolite", the file is
# "cosmos_pkn_enzyme_met.csv") -- read_local_cosmos_pkn()'s per-call
# `categories`/`filenames` args bridge the two.
pkn_categories <- c("allosteric", "enzyme_met", "grn", "ppi", "receptors", "transporters")
pkn_category_arg <- c(
    allosteric = "allosteric", enzyme_met = "enzyme_metabolite", grn = "grn",
    ppi = "ppi", receptors = "receptors", transporters = "transporters"
)

pkn_edges <- do.call(rbind, lapply(pkn_categories, function(cat) {
    raw <- read_local_cosmos_pkn(
        PKN_DIR, categories = pkn_category_arg[[cat]],
        filename_pattern = "cosmos_pkn_%s.csv",
        filenames = list(enzyme_metabolite = "cosmos_pkn_enzyme_met.csv")
    )
    edges <- finalize_cosmos_pkn(raw, source = "local")
    edges$category <- cat
    as.data.frame(edges)
}))
cat("PKN loaded:", nrow(pkn_edges), "edges across", length(pkn_categories), "categories\n")
print(table(pkn_edges$category))

## ---------------------------------------------------------------------
## 4.2 Exclude edges touching a malformed node ID -- recorded, not silent
## ---------------------------------------------------------------------
#
# A small number of snapshot rows don't match the 4-class COSMOS node
# grammar parse_cosmos_node_id() expects: 52 rows (enzyme_met/transporters)
# carry a KEGG Glycan ID mislabeled as ChEBI (e.g. CHEBI:G00001 -- not a
# real ChEBI accession, which is always numeric). Root cause is upstream in
# omnipath-metabo's GEM/KEGG resource processor; not investigated here
# (2026-10-07 decision) -- this filter keeps T009 unblocked and the
# exclusion explicit, matching FR-003's exclusion-recording pattern.
valid_node_id <- function(id) {
    grepl("^Metab__CHEBI:[0-9]+(_[a-z]+)?$", id) |
        grepl("^Gene[0-9]+__orphanReac.+$", id) |
        grepl("^Gene[0-9]+__.+$", id) |
        (!grepl("__", id, fixed = TRUE) & !startsWith(id, "Metab"))
}
malformed <- !valid_node_id(pkn_edges$source) | !valid_node_id(pkn_edges$target)
if (any(malformed)) {
    cat("\nExcluding", sum(malformed), "edges touching a malformed node ID (recorded, not silent):\n")
    bad_ids <- unique(c(pkn_edges$source[malformed], pkn_edges$target[malformed]))
    bad_ids <- bad_ids[!valid_node_id(bad_ids)]
    print(utils::head(bad_ids, 10))
    dir.create("result/pk_retrieval", recursive = TRUE, showWarnings = FALSE)
    write.csv(pkn_edges[malformed, ], "result/pk_retrieval/excluded_malformed_id_edges.csv", row.names = FALSE)
    pkn_edges <- pkn_edges[!malformed, ]
}

# T006 (reused, not new): resource/mechanism/edge_category provenance for
# SC-003, joined from the *_unformat.csv reference layers.
pkn_nodes <- data.frame(source = unique(c(pkn_edges$source, pkn_edges$target)), stringsAsFactors = FALSE)
annotated <- annotate_cosmos_network(pkn_nodes, pkn_edges, pkn_unformat_dir = PKN_DIR)
pkn_edges <- annotated$edges
names(pkn_edges)[names(pkn_edges) == "resource"] <- "evidence_refs"
cat("edges with evidence_refs populated:", sum(nzchar(pkn_edges$evidence_refs)), "of", nrow(pkn_edges), "\n")

# T007 (new): snapshot provenance, joined onto every edge by category.
snapshot_provenance <- record_pkn_snapshot_provenance(PKN_DIR, categories = pkn_categories)
pkn_edges <- merge(pkn_edges, snapshot_provenance, by = "category", all.x = TRUE)

dir.create("result/pk_retrieval", recursive = TRUE, showWarnings = FALSE)
saveRDS(pkn_edges, "result/pk_retrieval/pkn_edges.rds")
cat("\nSaved result/pk_retrieval/pkn_edges.rds:", nrow(pkn_edges), "rows,", ncol(pkn_edges), "cols\n")

## ---------------------------------------------------------------------
## 4.3 Load Morita et al. measured features + per-timepoint t-stat (T010)
## ---------------------------------------------------------------------
#
# FR-012's "t-value of the ob/ob-vs-WT comparison at each timepoint" is not
# a column the paper provides directly -- S2's DEA sheets are WITHIN-
# genotype time-course stats (does WT change over time; does ob/ob change
# over time), not a BETWEEN-genotype-at-a-timepoint comparison. Computed
# here from S1's raw per-sample values (5 WT + 5 ob/ob per timepoint,
# uniform across all 5 omics layers -- confirmed 2026-10-07) via Welch's
# t-test (no assumption of equal WT/ob-ob variance).
#
# transcript/phosphorylation aren't in 00_excel_to_h5.py's .h5 conversion
# (2026-10-07 decision): read directly from xlsx here for all 5 layers --
# 04_pk_retrieval.R is R, and 00/01's .h5 output is pandas-written, awkward
# to read from R, so reading xlsx directly (via readxl, already a
# case-study-2 dependency via 02_lipidID.r) is simpler and self-contained
# for every layer, not just the two missing ones.

sample_info <- as.data.frame(readxl::read_excel(
    "data/ads2547_data_file_s1.xlsx", sheet = "Sample_information"
))

omics_sheets <- c(
    metabolome = "metabolite", proteome = "protein", transcriptome = "transcript",
    phosphoproteome = "phosphorylation", plasma_metabolome = "plasma metabolite"
)
omics_tissue <- c(
    metabolome = "liver", proteome = "liver", transcriptome = "liver",
    phosphoproteome = "liver", plasma_metabolome = "plasma"
)

#' Welch's t-statistic, vectorized over rows (features) of two sample matrices.
welch_t_rows <- function(mat_wt, mat_ob) {
    n_wt <- rowSums(!is.na(mat_wt)); n_ob <- rowSums(!is.na(mat_ob))
    mean_wt <- rowMeans(mat_wt, na.rm = TRUE); mean_ob <- rowMeans(mat_ob, na.rm = TRUE)
    var_wt <- apply(mat_wt, 1, var, na.rm = TRUE); var_ob <- apply(mat_ob, 1, var, na.rm = TRUE)
    se <- sqrt(var_wt / n_wt + var_ob / n_ob)
    t_stat <- (mean_ob - mean_wt) / se
    t_stat[n_wt < 2 | n_ob < 2] <- NA_real_
    t_stat
}

measured_features <- do.call(rbind, lapply(names(omics_sheets), function(layer) {

    raw <- as.data.frame(readxl::read_excel("data/ads2547_data_file_s1.xlsx", sheet = omics_sheets[[layer]]))
    rownames(raw) <- raw$Row
    sample_matrix <- as.matrix(raw[, setdiff(names(raw), "Row"), drop = FALSE])

    timepoints <- sort(unique(sample_info$time))
    per_timepoint <- lapply(timepoints, function(tp) {
        wt_samples <- sample_info$MouseID[sample_info$time == tp & sample_info$genotype == "WT"]
        ob_samples <- sample_info$MouseID[sample_info$time == tp & sample_info$genotype == "ob/ob"]
        data.frame(
            feature_id = rownames(sample_matrix),
            timepoint_h = tp,
            t_stat = welch_t_rows(sample_matrix[, wt_samples, drop = FALSE], sample_matrix[, ob_samples, drop = FALSE]),
            stringsAsFactors = FALSE
        )
    })
    layer_result <- do.call(rbind, per_timepoint)
    layer_result$omics_layer <- layer
    layer_result$tissue <- omics_tissue[[layer]]
    layer_result

}))

cat("\nMeasured features:", nrow(measured_features), "rows (",
    length(unique(measured_features$feature_id)), "unique feature IDs x",
    length(unique(measured_features$timepoint_h)), "timepoints ) across",
    length(omics_sheets), "omics layers\n")
print(table(measured_features$omics_layer) / length(unique(measured_features$timepoint_h)))
cat("rows with a computable t_stat:", sum(!is.na(measured_features$t_stat)), "of", nrow(measured_features), "\n")

saveRDS(measured_features, "result/pk_retrieval/measured_features.rds")
cat("Saved result/pk_retrieval/measured_features.rds\n")
