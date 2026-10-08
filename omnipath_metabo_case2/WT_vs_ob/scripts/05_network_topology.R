# User Story 2 (spec 002-case-study-2-network, FR-006-010): metabolite-
# central topology with inter-organ circulation.
#
# Run from omnipath_metabo_case2/, after scripts/04_pk_retrieval.R.

suppressMessages(pkgload::load_all("../../Spatial-COSMOS-MISTy"))
source("WT_vs_ob/scripts/lib/pk_helpers.R")
suppressMessages(library(igraph))

PKN_DIR <- "../../Spatial-COSMOS-MISTy/data/PKN"
pkn_edges <- readRDS("WT_vs_ob/result/pk_retrieval/pkn_edges.rds")
measured_features <- readRDS("WT_vs_ob/result/pk_retrieval/measured_features.rds")

## ---------------------------------------------------------------------
## 5.1 Assert FR-006: metabolites appear as both source and target (T015)
## ---------------------------------------------------------------------
#
# Not newly constructed -- the PKN snapshot's receptors/allosteric/
# enzyme_met categories already carry both directions (confirmed via
# format_transporter_row()'s own fwd/rev expansion, research.md). This
# just verifies the property holds in the loaded network.

is_metab <- function(id) grepl("^Metab__", id)
metab_as_source <- unique(pkn_edges$source[is_metab(pkn_edges$source)])
metab_as_target <- unique(pkn_edges$target[is_metab(pkn_edges$target)])
metab_both <- intersect(metab_as_source, metab_as_target)
cat("T015 (FR-006): metabolite nodes as source:", length(metab_as_source),
    " as target:", length(metab_as_target),
    " as both (upstream AND downstream):", length(metab_both), "\n")
stopifnot(length(metab_both) > 0)

## ---------------------------------------------------------------------
## 5.2 Full chain metabolite->receptor->kinase->TF->GRN->enzyme<->metabolite (T016)
## ---------------------------------------------------------------------
#
# Searched constructively, category by category, rather than an
# unrestricted graph search over 202,858 edges: GRN targets and
# enzyme_met's "pre-expanded" (GEM-origin) nodes share bare-UniProt IDs
# directly (confirmed 2026-10-07: 873 overlapping nodes, exact literal-
# string matches, no Gene{N}__ prefix bridging needed) -- anchor there and
# extend backward through ppi (kinase->TF), ppi again (receptor->kinase),
# and receptors (metabolite->receptor).

receptors_edges <- pkn_edges[pkn_edges$category == "receptors", ]
ppi_edges <- pkn_edges[pkn_edges$category == "ppi", ]
grn_edges <- pkn_edges[pkn_edges$category == "grn", ]
enzyme_edges <- pkn_edges[pkn_edges$category == "enzyme_met", ]

grn_targets <- unique(grn_edges$target)

# Bare-UniProt GRN targets don't connect to enzyme_met directly -- they
# connect via a CONNECTOR edge (e.g. "Q8R3F5 -> Gene1189__Q8R3F5", mor=1,
# COSMOS's own internal identity-bridge between the bare-UniProt namespace
# and the Gene{N}__-prefixed namespace enzyme_met actually uses) to a
# Gene{N}__-prefixed relay node, and THAT relay node is what touches a
# metabolite. Discovered 2026-10-07 after the naive "bridge == enzyme_met
# source/target directly" assumption found zero reachable bridges.
connector_edges <- enzyme_edges[
    grepl("^Gene[0-9]+__", enzyme_edges$target) &
        sub("^Gene[0-9]+__|_rev$", "", enzyme_edges$target) == enzyme_edges$source,
]
bridge_proteins <- intersect(grn_targets, connector_edges$source)
cat("\nT016 (FR-007): candidate GRN-target / enzyme_met-connector bridge proteins:", length(bridge_proteins), "\n")

# Backward set-reachability per level, not a greedy single-path-per-bridge
# search (which failed: taking only the first TF/kinase/receptor candidate
# at each step and bailing on a dead end missed chains that exist via a
# DIFFERENT TF/kinase/receptor for the same bridge).
relay_nodes_with_metab <- unique(c(
    connector_edges$target[connector_edges$target %in% enzyme_edges$source[is_metab(enzyme_edges$target)]],
    connector_edges$target[connector_edges$target %in% enzyme_edges$target[is_metab(enzyme_edges$source)]]
))
level0_bridges <- unique(connector_edges$source[connector_edges$target %in% relay_nodes_with_metab])
level1_tfs <- unique(grn_edges$source[grn_edges$target %in% level0_bridges])
level2_kinases <- unique(ppi_edges$source[ppi_edges$target %in% level1_tfs])
level3_receptors <- unique(ppi_edges$source[ppi_edges$target %in% level2_kinases])
level4_metabolites <- unique(receptors_edges$source[receptors_edges$target %in% level3_receptors & is_metab(receptors_edges$source)])

cat("T016 backward reachability: ", length(level0_bridges), "connector-reachable bridges ->",
    length(level1_tfs), "TFs ->", length(level2_kinases), "kinases ->",
    length(level3_receptors), "receptors ->", length(level4_metabolites), "metabolites\n")

full_chain <- NULL
if (length(level4_metabolites) > 0) {
    start_metab <- level4_metabolites[1]
    receptor <- receptors_edges$target[receptors_edges$source == start_metab & receptors_edges$target %in% level3_receptors][1]
    kinase <- ppi_edges$target[ppi_edges$source == receptor & ppi_edges$target %in% level2_kinases][1]
    tf <- ppi_edges$target[ppi_edges$source == kinase & ppi_edges$target %in% level1_tfs][1]
    bridge <- grn_edges$target[grn_edges$source == tf & grn_edges$target %in% level0_bridges][1]
    relay <- connector_edges$target[connector_edges$source == bridge & connector_edges$target %in% relay_nodes_with_metab][1]
    end_metab <- if (relay %in% enzyme_edges$source) {
        enzyme_edges$target[enzyme_edges$source == relay & is_metab(enzyme_edges$target)][1]
    } else {
        enzyme_edges$source[enzyme_edges$target == relay & is_metab(enzyme_edges$source)][1]
    }

    full_chain <- data.frame(
        step = c("metabolite->receptor", "receptor->kinase", "kinase->TF", "TF->GRN target",
                  "GRN target->enzyme relay (connector)", "enzyme relay<->metabolite"),
        source = c(start_metab, receptor, kinase, tf, bridge, relay),
        target = c(receptor, kinase, tf, bridge, relay, end_metab),
        stringsAsFactors = FALSE
    )
}

if (is.null(full_chain)) {
    cat("T016: no full 5-hop chain found among", length(bridge_proteins), "bridge candidates (recorded, not fatal)\n")
} else {
    cat("T016: full chain found:\n")
    print(full_chain)
}

## ---------------------------------------------------------------------
## 5.3 Partial chain: direct metabolite<->enzyme, no signaling path (T017)
## ---------------------------------------------------------------------

signaling_proteins <- unique(c(receptors_edges$target, ppi_edges$source, ppi_edges$target))
direct_enzyme_edges <- enzyme_edges[
    (is_metab(enzyme_edges$source) & !(enzyme_edges$target %in% signaling_proteins) &
        !(sub("^Gene[0-9]+__", "", enzyme_edges$target) %in% signaling_proteins)) |
    (is_metab(enzyme_edges$target) & !(enzyme_edges$source %in% signaling_proteins) &
        !(sub("^Gene[0-9]+__", "", enzyme_edges$source) %in% signaling_proteins)),
]
cat("\nT017 (FR-008): direct metabolite<->enzyme edges with no signaling-pathway participant:",
    nrow(direct_enzyme_edges), "of", nrow(enzyme_edges), "enzyme_met edges\n")
stopifnot(nrow(direct_enzyme_edges) > 0)
partial_chain_example <- direct_enzyme_edges[1, ]
print(partial_chain_example)

## ---------------------------------------------------------------------
## 5.4 Liver<->blood metabolite compartment tagging + transporter edges (T018)
## ---------------------------------------------------------------------
#
# research.md R6: cosmosR's compartment codes (c/e/m/...) are subcellular,
# not organ-level -- liver/blood tagging is a case-study-2-specific
# augmentation, distinct from and in addition to that GEM compartment
# suffix. "Liver" metabolites = anything measured in the liver metabolome
# layer; "blood" = anything measured in the plasma_metabolome layer.
# Transporter edges (FR-009) use the already-loaded transporters category
# PKN edges as the mechanism connecting the two compartments for a shared
# ChEBI identity.

liver_chebi <- unique(na.omit(measured_features$pk_node_id[
    measured_features$omics_layer == "metabolome" & measured_features$mapping_status == "mapped"
]))
blood_chebi <- unique(na.omit(measured_features$pk_node_id[
    measured_features$omics_layer == "plasma_metabolome" & measured_features$mapping_status == "mapped"
]))
shared_chebi <- intersect(liver_chebi, blood_chebi)
cat("\nT018 (FR-009/010): liver-measured ChEBI:", length(liver_chebi),
    " blood-measured ChEBI:", length(blood_chebi),
    " shared (both tissues):", length(shared_chebi), "\n")

compartment_tagged_nodes <- rbind(
    data.frame(chebi_id = liver_chebi, tissue = "liver", node_id = paste0(liver_chebi, "_liver"), stringsAsFactors = FALSE),
    data.frame(chebi_id = blood_chebi, tissue = "blood", node_id = paste0(blood_chebi, "_blood"), stringsAsFactors = FALSE)
)

transporter_edges <- pkn_edges[pkn_edges$category == "transporters", ]
transporter_metab_chebi <- unique(sub("^Metab__(CHEBI:[0-9]+)(_[a-z]+)?$", "\\1",
    grep("^Metab__CHEBI:", c(transporter_edges$source, transporter_edges$target), value = TRUE)))
shared_with_transporter <- intersect(shared_chebi, transporter_metab_chebi)
cat("T018: shared liver/blood ChEBI with a transporters-category PKN edge:",
    length(shared_with_transporter), "of", length(shared_chebi), "\n")
stopifnot(length(shared_with_transporter) > 0)

liver_blood_transporter_edges <- data.frame(
    source = paste0(shared_with_transporter, "_liver"),
    target = paste0(shared_with_transporter, "_blood"),
    mor = 1L,
    category = "transporters_inter_organ",
    stringsAsFactors = FALSE
)
cat("T018: liver<->blood transporter edges added:", nrow(liver_blood_transporter_edges), "\n")

dir.create("WT_vs_ob/result/networks", recursive = TRUE, showWarnings = FALSE)
saveRDS(compartment_tagged_nodes, "WT_vs_ob/result/networks/compartment_tagged_metabolite_nodes.rds")
saveRDS(liver_blood_transporter_edges, "WT_vs_ob/result/networks/liver_blood_transporter_edges.rds")
if (!is.null(full_chain)) saveRDS(full_chain, "WT_vs_ob/result/networks/example_full_chain.rds")
saveRDS(partial_chain_example, "WT_vs_ob/result/networks/example_partial_chain.rds")

cat("\nSaved result/networks/{compartment_tagged_metabolite_nodes,liver_blood_transporter_edges}.rds\n")
