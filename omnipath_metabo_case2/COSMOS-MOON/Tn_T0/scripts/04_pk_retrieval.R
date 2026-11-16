# User Story 1 (spec 002-case-study-2-network, FR-001-005): retrieve PK for
# every Morita et al. measured feature from the local omnipath-metabo PKN
# snapshot, and report node/edge counts against Morita et al. Figure 2.
#
# Run from omnipath_metabo_case2/, with Spatial-COSMOS-MISTy as a sibling
# checkout (quickstart.md).

suppressMessages(pkgload::load_all("../../Spatial-COSMOS-MISTy"))
source("Tn_T0/scripts/lib/pk_helpers.R")

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
    dir.create("Tn_T0/result/pk_retrieval", recursive = TRUE, showWarnings = FALSE)
    write.csv(pkn_edges[malformed, ], "Tn_T0/result/pk_retrieval/excluded_malformed_id_edges.csv", row.names = FALSE)
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

dir.create("Tn_T0/result/pk_retrieval", recursive = TRUE, showWarnings = FALSE)
saveRDS(pkn_edges, "Tn_T0/result/pk_retrieval/pkn_edges.rds")
cat("\nSaved result/pk_retrieval/pkn_edges.rds:", nrow(pkn_edges), "rows,", ncol(pkn_edges), "cols\n")

## ---------------------------------------------------------------------
## 4.3 Load Morita et al. measured features + within-genotype Tn-vs-T0
## t-stat, per genotype (spec 003-case-study-2-temporal-trajectory FR-001-004)
## ---------------------------------------------------------------------
#
# Replaces 002's single per-timepoint ob/ob-vs-WT t-value (one comparison
# per timepoint, genotype encoded inside the comparison direction) with 14
# within-genotype t-values: for each genotype (WT, ob/ob) separately, each
# non-zero timepoint's 5 samples vs. that SAME genotype's 5 samples at 0h
# (research.md R1). timepoint_h == 0 has no row -- it's the baseline every
# other timepoint is compared against, not a comparison target itself
# (FR-002). welch_t_rows() is 002's own function, reused UNMODIFIED
# (research.md R2) -- it was already generic over "compare two matched
# sample groups"; only the indexing of which two groups changes here.
# Sampling is destructive/independent per timepoint (confirmed against
# sample_info during 002's work: different mice at each timepoint), so
# this is an unpaired two-sample t-test, same as 002's, not a paired one
# (research.md R3).
#
# transcript/phosphorylation aren't in 00_excel_to_h5.py's .h5 conversion
# (2026-10-07 decision, carried over unchanged): read directly from xlsx
# here for all 5 layers.

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
#' Generic over "compare two matched sample groups" -- 002's call site passed
#' (ob samples at t, WT samples at t); this feature's call site passes
#' (genotype samples at t, that same genotype's samples at t=0). Same
#' function, different indexing only (research.md R2).
welch_t_rows <- function(mat_a, mat_b) {
    n_a <- rowSums(!is.na(mat_a)); n_b <- rowSums(!is.na(mat_b))
    mean_a <- rowMeans(mat_a, na.rm = TRUE); mean_b <- rowMeans(mat_b, na.rm = TRUE)
    var_a <- apply(mat_a, 1, var, na.rm = TRUE); var_b <- apply(mat_b, 1, var, na.rm = TRUE)
    se <- sqrt(var_a / n_a + var_b / n_b)
    t_stat <- (mean_a - mean_b) / se
    t_stat[n_a < 2 | n_b < 2] <- NA_real_
    t_stat
}

# Genotype token convention (data-model.md): "WT"/"ob" in every output key
# and filename (avoids "/" from the source data's own "ob/ob" as a path
# separator); gt_raw is the literal value as it appears in sample_info.
genotype_tokens <- c(WT = "WT", ob = "ob/ob")

measured_features <- do.call(rbind, lapply(names(omics_sheets), function(layer) {

    raw <- as.data.frame(readxl::read_excel("data/ads2547_data_file_s1.xlsx", sheet = omics_sheets[[layer]]))
    rownames(raw) <- raw$Row
    sample_matrix <- as.matrix(raw[, setdiff(names(raw), "Row"), drop = FALSE])

    nonzero_timepoints <- sort(setdiff(unique(sample_info$time), 0))

    per_genotype <- lapply(names(genotype_tokens), function(gt_token) {
        gt_raw <- genotype_tokens[[gt_token]]
        baseline_samples <- sample_info$MouseID[sample_info$time == 0 & sample_info$genotype == gt_raw]
        per_timepoint <- lapply(nonzero_timepoints, function(tp) {
            tp_samples <- sample_info$MouseID[sample_info$time == tp & sample_info$genotype == gt_raw]
            data.frame(
                feature_id = rownames(sample_matrix),
                genotype = gt_token,
                timepoint_h = tp,
                t_stat = welch_t_rows(
                    sample_matrix[, tp_samples, drop = FALSE],
                    sample_matrix[, baseline_samples, drop = FALSE]
                ),
                stringsAsFactors = FALSE
            )
        })
        do.call(rbind, per_timepoint)
    })
    layer_result <- do.call(rbind, per_genotype)
    layer_result$omics_layer <- layer
    layer_result$tissue <- omics_tissue[[layer]]
    layer_result

}))
measured_features$excluded <- FALSE

cat("\nMeasured features:", nrow(measured_features), "rows (",
    length(unique(measured_features$feature_id)), "unique feature IDs x",
    length(unique(measured_features$genotype)), "genotypes x",
    length(unique(measured_features$timepoint_h)), "non-zero timepoints ) across",
    length(omics_sheets), "omics layers\n")
print(table(measured_features$omics_layer) / (length(unique(measured_features$genotype)) * length(unique(measured_features$timepoint_h))))
cat("rows with a computable t_stat:", sum(!is.na(measured_features$t_stat)), "of", nrow(measured_features), "\n")

## ---------------------------------------------------------------------
## 4.4 Record (not silently drop) excluded lipid/FFA/acyl species (T011, FR-003)
## ---------------------------------------------------------------------
#
# "lipid"/"FFAandAcyls" use the same clean Row-ID scheme as the included
# layers (e.g. "lipid;CE;C02530") and are the direct per-timepoint
# counterparts excluded here. "Lipid_all"/"Acylcarnitine_AcylCoA_all" are a
# fuller, raw species-level breakdown of this same excluded scope (the
# complex lipid-name-parsing pipeline 02_lipidID.r already targets) --
# not loaded again here, to avoid double-counting the same exclusion.
excluded_sheets <- c(lipid = "lipid", ffa_acyl = "FFAandAcyls")
excluded_features <- do.call(rbind, lapply(names(excluded_sheets), function(layer) {
    raw <- as.data.frame(readxl::read_excel("data/ads2547_data_file_s1.xlsx", sheet = excluded_sheets[[layer]]))
    data.frame(
        feature_id = raw$Row, genotype = NA_character_, timepoint_h = NA_integer_, t_stat = NA_real_,
        omics_layer = layer, tissue = "liver", excluded = TRUE,
        stringsAsFactors = FALSE
    )
}))
cat("\nExcluded (FR-003, recorded not dropped):", nrow(excluded_features),
    "lipid/FFA features (acyl-CoA/acyl-carnitine species breakdown deferred to 02_lipidID.r)\n")

measured_features <- rbind(measured_features, excluded_features)

## ---------------------------------------------------------------------
## 4.5 Resolve each non-excluded feature against the PK entity set (T012)
## ---------------------------------------------------------------------
#
# Two distinct translations, matching each layer's native ID space in the
# Row column (research.md's T012 scope): metabolome/plasma_metabolome
# carry a KEGG Compound ID (e.g. "metabolite;...;C01035") -- the PKN uses
# ChEBI, so this needs KEGGREST::keggConv() (confirmed live 2026-10-07, one
# bulk call for all unique compounds). proteome/transcriptome/
# phosphoproteome carry an Entrez ID via "mmu:<id>" (phosphosite rows also
# carry a "-<site>" suffix, stripped here) -- the PKN's kinase-substrate
# layers are protein-level only (confirmed: no site-suffixed node IDs
# anywhere in the snapshot), so phosphosites resolve to their parent
# protein's UniProt accession via org.Mm.eg.db, same as proteome/
# transcriptome; site-level specificity is lost at the PK-matching step,
# an inherent limitation of the PKN's own notation, not of this mapping.

measured_features$pk_node_id <- NA_character_

# PKN-side bare-ID lookup, reusing parse_cosmos_node_id() (T005) rather
# than re-deriving the COSMOS node grammar by hand.
pkn_node_ids <- unique(c(pkn_edges$source, pkn_edges$target))
pkn_node_info <- lapply(pkn_node_ids, parse_cosmos_node_id)
pkn_chebi_set <- unique(stats::na.omit(vapply(pkn_node_info, function(x) x$chebi_id, character(1))))
pkn_uniprot_set <- unique(stats::na.omit(vapply(pkn_node_info, function(x) x$uniprot_id, character(1))))

# Reverse lookup: bare ChEBI/UniProt -> every actual PKN node ID sharing it
# (multiple compartments for a metabolite; multiple Gene{N}__ indices for a
# protein) -- needed by T013 to expand a resolved bare ID into the PKN nodes
# it actually touches.
bare_ids <- vapply(pkn_node_info, function(x) if (!is.na(x$chebi_id)) x$chebi_id else x$uniprot_id, character(1))
bare_to_pkn_nodes <- split(pkn_node_ids, bare_ids)

# -- Metabolite layers: KEGG Compound -> ChEBI
metab_rows <- measured_features$omics_layer %in% c("metabolome", "plasma_metabolome") & !measured_features$excluded
kegg_ids <- sub("^.*;([A-Za-z0-9]+)$", "\\1", measured_features$feature_id[metab_rows])
unique_kegg <- unique(kegg_ids)
kegg_to_chebi <- tryCatch({
    conv <- KEGGREST::keggConv("chebi", paste0("cpd:", unique_kegg))
    setNames(sub("^chebi:", "", conv), sub("^cpd:", "", names(conv)))
}, error = function(e) {
    warning("KEGGREST lookup failed, metabolites will be unmapped: ", conditionMessage(e))
    character(0)
})
resolved_chebi <- unname(kegg_to_chebi[kegg_ids])
measured_features$pk_node_id[metab_rows] <- ifelse(is.na(resolved_chebi), NA_character_, paste0("CHEBI:", resolved_chebi))
cat("\nKEGG->ChEBI:", length(unique_kegg), "unique compounds,",
    sum(!is.na(kegg_to_chebi[unique_kegg])), "resolved\n")

# -- Gene/protein layers: Entrez -> UniProt (preferring a candidate already in the PKN)
gene_rows <- measured_features$omics_layer %in% c("proteome", "transcriptome", "phosphoproteome") & !measured_features$excluded
entrez_ids <- sub("^.*mmu:([0-9]+)(-[^;]*)?$", "\\1", measured_features$feature_id[gene_rows])
suppressMessages(library(org.Mm.eg.db))
uniprot_map_df <- AnnotationDbi::select(org.Mm.eg.db, keys = unique(entrez_ids), keytype = "ENTREZID", columns = "UNIPROT")
uniprot_map_df <- uniprot_map_df[!is.na(uniprot_map_df$UNIPROT), ]
entrez_to_uniprot_candidates <- split(uniprot_map_df$UNIPROT, uniprot_map_df$ENTREZID)
resolve_uniprot <- function(entrez_id) {
    candidates <- entrez_to_uniprot_candidates[[entrez_id]]
    if (is.null(candidates)) return(NA_character_)
    in_pkn <- candidates[candidates %in% pkn_uniprot_set]
    if (length(in_pkn) > 0) in_pkn[1] else candidates[1]
}
measured_features$pk_node_id[gene_rows] <- vapply(entrez_ids, resolve_uniprot, character(1))
cat("Entrez->UniProt:", length(unique(entrez_ids)), "unique genes,",
    sum(!is.na(vapply(unique(entrez_ids), resolve_uniprot, character(1)))), "resolved\n")

# -- mapping_status: distinguish "resolved to an ID, but that ID isn't a PKN
# node" from "excluded_scope" (FR-003) from "excluded by this cycle's exclusion list"
measured_features$mapping_status <- ifelse(
    measured_features$excluded, "excluded_scope",
    ifelse(
        is.na(measured_features$pk_node_id), "unmapped_no_pk_entry",
        ifelse(
            measured_features$omics_layer %in% c("metabolome", "plasma_metabolome"),
            ifelse(measured_features$pk_node_id %in% pkn_chebi_set, "mapped", "unmapped_no_pk_entry"),
            ifelse(measured_features$pk_node_id %in% pkn_uniprot_set, "mapped", "unmapped_no_pk_entry")
        )
    )
)
cat("\nmapping_status:\n")
print(table(measured_features$mapping_status))

saveRDS(measured_features, "Tn_T0/result/pk_retrieval/measured_features.rds")
cat("\nSaved result/pk_retrieval/measured_features.rds:", nrow(measured_features), "rows\n")

## ---------------------------------------------------------------------
## 4.6 Build the two per-genotype GenotypeNetworks (T013)
## ---------------------------------------------------------------------
#
# The measured feature *panel* is identical across genotypes (confirmed
# 2026-10-07: same 80 samples, same panel), so what actually makes a WT
# network differ from an ob/ob network is which features are significant
# IN THAT GENOTYPE's own time-course -- from S2's within-genotype DEA
# calls (WT_change/ob/ob_change in {NS, Up, Down}), not the between-
# genotype t_stat T010 computed (that's MOON's input, used later,
# per-timepoint not per-genotype -- see spec.md FR-013 correction).

dea_sheets <- c(
    metabolome = "metabolite", proteome = "protein", transcriptome = "transcript",
    phosphoproteome = "phosphorylation", plasma_metabolome = "plasma metabolite"
)
dea_calls <- do.call(rbind, lapply(names(dea_sheets), function(layer) {
    raw <- as.data.frame(readxl::read_excel("data/ads2547_data_file_s2.xlsx", sheet = dea_sheets[[layer]]))
    data.frame(
        feature_id = raw$Row, WT_change = raw$WT_change, ob_ob_change = raw[["ob/ob_change"]],
        stringsAsFactors = FALSE
    )
}))

feature_pk <- unique(measured_features[!measured_features$excluded, c("feature_id", "pk_node_id", "mapping_status")])
feature_pk <- merge(feature_pk, dea_calls, by = "feature_id", all.x = TRUE)

build_genotype_network <- function(change_col) {
    is_significant <- feature_pk[[change_col]] %in% c("Up", "Down") & feature_pk$mapping_status == "mapped"
    significant_bare_ids <- unique(feature_pk$pk_node_id[is_significant])
    genotype_pkn_nodes <- unique(unlist(bare_to_pkn_nodes[significant_bare_ids]))
    edges <- pkn_edges[pkn_edges$source %in% genotype_pkn_nodes | pkn_edges$target %in% genotype_pkn_nodes, ]
    nodes <- unique(c(edges$source, edges$target))
    list(nodes = nodes, edges = edges, n_significant_features = length(significant_bare_ids))
}

genotype_networks <- list(WT = build_genotype_network("WT_change"), ob_ob = build_genotype_network("ob_ob_change"))

dir.create("Tn_T0/result/networks", recursive = TRUE, showWarnings = FALSE)
for (g in names(genotype_networks)) {
    net <- genotype_networks[[g]]
    cat(sprintf(
        "\n%s network: %d significant mapped features -> %d nodes, %d edges\n",
        g, net$n_significant_features, length(net$nodes), nrow(net$edges)
    ))
    write.csv(net$edges, sprintf("Tn_T0/result/networks/%s_network_edges.csv", g), row.names = FALSE)
    saveRDS(net, sprintf("Tn_T0/result/networks/%s_network.rds", g))
}
