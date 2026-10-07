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
