# User Story 3 (T024, FR-016): PACON pathway control analysis -- which
# pathways each high-|score| driver node controls downstream, per timepoint.
#
# Two analyses, run back to back:
#
# 1. Gene-centric PACON (pathway_control_analysis(), existing package
#    function): MSigDB C2 NABA_/KEGG_ gene sets, gene-only ORA. First pass
#    (2026-10-08) found only 1 of 80 tested pathways across all 8
#    timepoints was a real metabolic-process pathway (vs. signaling/
#    immune/ECM) -- root-caused to a structural gap, not a tunable
#    parameter: GEM edges are deliberately excluded from moon()'s own
#    scoring graph (FR-014/015 sign exemption), so most GEM-connected
#    enzyme genes never get a MOON score and are invisible to this
#    function's gene universe (confirmed: of 570 GEM-edge enzyme proteins
#    at 12h, only 11 were in the MOON-scored set) -- widening n_steps from
#    2/3 to 3/4 made no difference, as expected once the cause is known.
#
# 2. Metabolite-centric GEM-PACON (gem_pathway_control_analysis(), new
#    package function, 2026-10-08): mirrors (1) but for metabolites --
#    drivers are high-|score| metabolites, neighborhoods are walked over
#    the GEM edge graph directly (not the signaling PKN), background is
#    every metabolite in the FULL GEM edge set (not MOON-scored nodes --
#    this is the fix for (1)'s structural gap, applied metabolite-side),
#    and pathways come from a new KEGG compound->pathway table
#    (build_cosmos_metabolite_pathway_sets(), live KEGG query, cached).
#
# Operates directly on reduce_moon_network()'s ALREADY-thresholded output
# (per-timepoint primary=1.35/secondary=0.90, 06_footprint_moon.R), not the
# raw unthresholded moon_res/pruned_pkn the spatial pilot's own
# 09_pathway_control.R uses with its own separate top_node_threshold=3/
# background_threshold=1.5 -- those pilot defaults were tuned for a
# different (single-sample-per-spot) score distribution and don't transfer
# here (confirmed during the MOON threshold sweep: 3/2 collapsed most of
# this dataset's timepoints to ~0 nodes). Concretely: moon_res is
# restricted to source_original %in% pruned$nodes$source before being
# handed to pathway_control_analysis(), and pruned_pkn is the final
# (edges + reattached gem_edges) display network, not the broader n-step-
# reachable candidate graph -- so top_node_threshold=1.35 reuses the same
# primary threshold pruning already applied, and background_threshold=0
# means "background = every node that survived pruning" (they already
# cleared >= secondary=0.90), rather than re-filtering a second time on a
# different scale.
#
# Run from omnipath_metabo_case2/, after scripts/06_footprint_moon.R.

suppressMessages(pkgload::load_all("../../Spatial-COSMOS-MISTy"))
suppressMessages(library(org.Mm.eg.db))
source("WT_vs_ob/scripts/lib/pk_helpers.R")

moon_results <- readRDS("WT_vs_ob/result/moon/all_timepoints.rds")
pkn_edges <- readRDS("WT_vs_ob/result/pk_retrieval/pkn_edges.rds")
measured_features <- readRDS("WT_vs_ob/result/pk_retrieval/measured_features.rds")

## ---------------------------------------------------------------------
## ID -> label maps, built once over the whole PKN (not per-timepoint) so
## coverage is consistent across all 8 pathway-control runs.
## ---------------------------------------------------------------------

all_nodes <- unique(c(pkn_edges$source, pkn_edges$target))
uniprot_nodes <- all_nodes[grepl("^[A-Z0-9]+$", all_nodes)]  # excludes Metab__..., Gene<N>__orphanReac...
symbol_df <- AnnotationDbi::select(org.Mm.eg.db, keys = uniprot_nodes, keytype = "UNIPROT", columns = "SYMBOL")
symbol_df <- symbol_df[!is.na(symbol_df$SYMBOL) & !duplicated(symbol_df$UNIPROT), ]
uniprot_symbol_map <- stats::setNames(symbol_df$SYMBOL, symbol_df$UNIPROT)
cat("UniProt -> symbol:", length(uniprot_symbol_map), "of", length(uniprot_nodes), "PKN protein nodes\n")

met_rows <- measured_features[measured_features$omics_layer %in% c("metabolome", "plasma_metabolome") & !is.na(measured_features$pk_node_id), ]
met_parts <- strsplit(met_rows$feature_id, ";", fixed = TRUE)
met_name <- vapply(met_parts, function(p) if (length(p) >= 2) p[2] else NA_character_, character(1))
chebi_name_map <- stats::setNames(met_name, met_rows$pk_node_id)
chebi_name_map <- chebi_name_map[!duplicated(names(chebi_name_map))]
cat("ChEBI -> name:", length(chebi_name_map), "measured metabolites\n")

## ---------------------------------------------------------------------
## Pathway gene sets -- current package default (MSigDB C2, NABA_/KEGG_
## patterns). Revisit after seeing this first pass's hit rate.
## ---------------------------------------------------------------------

pathways <- build_cosmos_pathway_gene_sets()
cat("Pathway gene sets:", length(unique(pathways$source)), "pathways (NABA_/KEGG_),",
    nrow(pathways), "pathway-gene rows\n")

## ---------------------------------------------------------------------
## Per-timepoint PACON, on the already-pruned network
## ---------------------------------------------------------------------

dir.create("WT_vs_ob/result/pacon", recursive = TRUE, showWarnings = FALSE)

run_one_timepoint <- function(tp, pruned) {

    moon_res_pruned <- pruned$moon_res[pruned$moon_res$source_original %in% pruned$nodes$source, ]
    pkn_pruned <- rbind(
        pruned$edges[, c("source", "target", "interaction")],
        data.frame(source = pruned$gem_edges$source, target = pruned$gem_edges$target,
                   interaction = pruned$gem_edges$mor, stringsAsFactors = FALSE)
    )

    result <- tryCatch(
        pathway_control_analysis(
            moon_res = moon_res_pruned,
            pruned_pkn = pkn_pruned,
            pathways = pathways,
            uniprot_symbol_map = uniprot_symbol_map,
            chebi_name_map = chebi_name_map,
            top_node_threshold = 1.35,   # == 06_footprint_moon.R's primary_thresh
            background_threshold = 0,   # background = every pruned-network node (already >= secondary_thresh=0.90)
            # Widened from pilot defaults (2/3, 2026-10-08): first pass found
            # only 1 of 80 tested pathways across all 8 timepoints was a real
            # metabolic-process pathway (vs. signaling/immune/ECM) -- drivers
            # are mostly TF/kinase upstream_input nodes, whose 2-hop
            # neighborhood rarely reaches past the enzyme layer into GEM.
            # Widening the search depth lets more driver neighborhoods reach
            # actual metabolic-enzyme genes.
            n_steps = 3, n_steps_metabolite = 4
        ),
        error = function(e) {
            cat(sprintf("  timepoint %sh: FAILED (%s)\n", tp, conditionMessage(e)))
            NULL
        }
    )
    if (is.null(result)) return(NULL)

    cat(sprintf(
        "  timepoint %sh: %d drivers (of %d pruned nodes), %d pathway x driver ORA rows\n",
        tp, nrow(result$drivers), nrow(pruned$nodes), nrow(result$pathway_control)
    ))
    result
}

pacon_results <- list()
for (tp in names(moon_results)) {
    cat(sprintf("\n--- timepoint %sh ---\n", tp))
    r <- run_one_timepoint(tp, moon_results[[tp]])
    if (!is.null(r)) {
        pacon_results[[tp]] <- r
        write.csv(r$pathway_control, sprintf("WT_vs_ob/result/pacon/%sh_pathway_control.csv", tp), row.names = FALSE)
        write.csv(r$drivers, sprintf("WT_vs_ob/result/pacon/%sh_drivers.csv", tp), row.names = FALSE)
        heat <- tryCatch(
            plot_pathway_control_heatmap(r, pval_threshold = 0.01, min_hits = 5,
                                          path = "WT_vs_ob/result/pacon", plot_name = sprintf("%sh_PACON_heatmap", tp)),
            warning = function(w) { cat("  heatmap:", conditionMessage(w), "\n"); NULL }
        )
    }
}

saveRDS(pacon_results, "WT_vs_ob/result/pacon/all_timepoints.rds")
cat("\nPACON (gene-centric) complete:", length(pacon_results), "of", length(moon_results), "timepoints\n")
cat("Saved WT_vs_ob/result/pacon/{<timepoint>h_pathway_control,<timepoint>h_drivers}.csv, WT_vs_ob/result/pacon/all_timepoints.rds\n")

## ---------------------------------------------------------------------
## GEM-PACON: metabolite-centered pathway control (2026-10-08)
## ---------------------------------------------------------------------

cat("\n\n=== GEM-PACON (metabolite-centered) ===\n")

met_pathways <- build_cosmos_metabolite_pathway_sets(cache_path = "WT_vs_ob/result/pacon/metabolite_pathways_cache.rds")
cat("Metabolite pathway sets:", length(unique(met_pathways$source)), "KEGG metabolism pathways,",
    nrow(met_pathways), "pathway-compound rows\n")

# Full (pre-split, whole-PKN) GEM edge set -- NOT the per-timepoint
# reattached subset (that's only ~18 unique metabolites, too sparse a
# background/graph for ORA; confirmed via a direct comparison: full GEM
# has 2,775 unique metabolites vs. 18 in a per-timepoint reattached
# subset). Same full GEM set, shared across all 8 timepoints, so results
# are comparable across timepoints -- only the driver selection (from
# each timepoint's own moon_res) varies.
split <- split_gem_edges(pkn_edges)
full_gem_edges <- split$gem[, c("source", "target", "mor")]
cat("Full GEM edge graph:", nrow(full_gem_edges), "edges,",
    length(unique(c(full_gem_edges$source, full_gem_edges$target))), "nodes\n")

# Tried narrowing the ORA background to only measured metabolites
# (restrict_background_to -- still available on gem_pathway_control_analysis()
# for other uses), on the theory that the full GEM network's 2,775 mostly-
# unmeasured reference compounds let a handful of large generic pathways
# (every amino acid, every nucleotide) always win regardless of which
# drivers are active. Measured empirically against this data: it made
# cross-timepoint differentiation WORSE, not better (top-10 overlap rose
# to 9-10/10, from 8-10/10 on the full background) -- a smaller universe
# concentrates stats on a pathway's raw popularity even more, it doesn't
# dilute it. Reverted to the full-GEM-network background; differentiation
# is handled downstream instead (see 07b, breadth-across-timepoints view,
# not each timepoint's own absolute top-N).

run_one_timepoint_gem <- function(tp, pruned) {

    result <- tryCatch(
        gem_pathway_control_analysis(
            moon_res = pruned$moon_res,
            gem_edges = full_gem_edges,
            pathways = met_pathways,
            chebi_name_map = chebi_name_map,
            top_node_threshold = 1.35,  # == 06_footprint_moon.R's primary_thresh
            n_steps_metabolite = 2     # metabolites sharing 1 enzyme/reaction with the driver
        ),
        error = function(e) {
            cat(sprintf("  timepoint %sh: FAILED (%s)\n", tp, conditionMessage(e)))
            NULL
        }
    )
    if (is.null(result)) return(NULL)

    cat(sprintf(
        "  timepoint %sh: %d metabolite drivers, %d pathway x driver ORA rows\n",
        tp, nrow(result$drivers), nrow(result$pathway_control)
    ))
    result
}

gem_pacon_results <- list()
for (tp in names(moon_results)) {
    cat(sprintf("\n--- timepoint %sh ---\n", tp))
    r <- run_one_timepoint_gem(tp, moon_results[[tp]])
    if (!is.null(r)) {
        gem_pacon_results[[tp]] <- r
        write.csv(r$pathway_control, sprintf("WT_vs_ob/result/pacon/%sh_metabolite_pathway_control.csv", tp), row.names = FALSE)
        write.csv(r$drivers, sprintf("WT_vs_ob/result/pacon/%sh_metabolite_drivers.csv", tp), row.names = FALSE)
        heat <- tryCatch(
            plot_pathway_control_heatmap(r, pval_threshold = 0.001, min_hits = 10,
                                          path = "WT_vs_ob/result/pacon", plot_name = sprintf("%sh_GEM_PACON_heatmap", tp)),
            warning = function(w) { cat("  heatmap:", conditionMessage(w), "\n"); NULL }
        )
    }
}

saveRDS(gem_pacon_results, "WT_vs_ob/result/pacon/metabolite_all_timepoints.rds")
cat("\nGEM-PACON complete:", length(gem_pacon_results), "of", length(moon_results), "timepoints\n")
cat("Saved WT_vs_ob/result/pacon/{<timepoint>h_metabolite_pathway_control,<timepoint>h_metabolite_drivers}.csv, WT_vs_ob/result/pacon/metabolite_all_timepoints.rds\n")
