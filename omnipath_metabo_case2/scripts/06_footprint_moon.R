# User Story 3 (spec 002-case-study-2-network, FR-011-015): per-timepoint
# footprint-based TF/kinase activity inference + COSMOS-MOON pruning.
#
# One network per timepoint, not per genotype (corrected 2026-10-07):
# FR-012's single ob/ob-vs-WT t-value per timepoint is the only MOON input
# signal -- genotype is already encoded inside it.
#
# Run from omnipath_metabo_case2/, after scripts/04_pk_retrieval.R and
# scripts/05_network_topology.R.

suppressMessages(pkgload::load_all("../../Spatial-COSMOS-MISTy"))
source("scripts/lib/pk_helpers.R")

pkn_edges <- readRDS("result/pk_retrieval/pkn_edges.rds")
measured_features <- readRDS("result/pk_retrieval/measured_features.rds")

## ---------------------------------------------------------------------
## 6.1 GEM exemption (T021, FR-014/015): split before MOON ever sees it
## ---------------------------------------------------------------------

split <- split_gem_edges(pkn_edges)
pkn_for_moon <- split$non_gem[, c("source", "target", "mor")]
grn_for_moon <- pkn_edges[pkn_edges$category == "grn", c("source", "target", "mor")]
gem_edges <- split$gem[, c("source", "target", "mor")]
cat("PKN for MOON (GEM excluded):", nrow(pkn_for_moon), "edges | GEM (reattached post-pruning):", nrow(gem_edges), "edges\n")

# All subcellular compartment codes observed in the snapshot (2026-10-07),
# not just the pilot's c('e','c') -- maximizes which metabolite t_stat
# values actually reach a real PKN node.
all_compartments <- c("c", "e", "eg", "g", "i", "l", "m", "n", "r", "v", "x")

## ---------------------------------------------------------------------
## 6.2 Build per-timepoint upstream/downstream inputs (T020, research.md R5)
## ---------------------------------------------------------------------

build_moon_inputs <- function(tp) {
    tp_features <- measured_features[
        measured_features$timepoint_h == tp & !measured_features$excluded &
            measured_features$mapping_status == "mapped" & !is.na(measured_features$t_stat),
    ]
    upstream_rows <- tp_features[tp_features$omics_layer %in% c("metabolome", "plasma_metabolome"), ]
    downstream_rows <- tp_features[tp_features$omics_layer %in% c("proteome", "transcriptome", "phosphoproteome"), ]

    # Bare ChEBI (no Metab__ prefix) -- prepare_metab_inputs() does the
    # compartment-suffix expansion inside run_moon_scoring(), the same
    # mechanism the original pilot uses (not T018's _liver/_blood tags,
    # which are the FR-009/010 topology-demonstration artifact, already
    # verified separately in 05_network_topology.R).
    upstream_input <- stats::setNames(upstream_rows$t_stat, upstream_rows$pk_node_id)
    upstream_input <- upstream_input[!duplicated(names(upstream_input))]

    downstream_input <- stats::setNames(downstream_rows$t_stat, downstream_rows$pk_node_id)
    downstream_input <- downstream_input[!duplicated(names(downstream_input))]

    list(upstream_input = upstream_input, downstream_input = downstream_input)
}

## ---------------------------------------------------------------------
## 6.3 Run MOON scoring + pruning, reattach GEM edges (T021/T022)
## ---------------------------------------------------------------------

run_one_timepoint <- function(tp) {
    inputs <- build_moon_inputs(tp)
    cat(sprintf(
        "\n--- timepoint %sh: %d upstream (metabolite) + %d downstream (gene/protein) inputs ---\n",
        tp, length(inputs$upstream_input), length(inputs$downstream_input)
    ))

    moon_scoring_result <- run_moon_scoring(
        node_activities = inputs, pkn = pkn_for_moon, grn = grn_for_moon,
        n_steps = 6, statistic = "ulm", compartments = all_compartments
    )
    pruned <- reduce_moon_network(moon_scoring_result, primary_thresh = 1.5, secondary_thresh = 1)
    pruned <- reattach_gem_edges(pruned, gem_edges)
    pruned$timepoint_h <- tp

    cat(sprintf(
        "timepoint %sh: %d scored nodes -> pruned to %d nodes, %d edges (+ %d GEM edges reattached)\n",
        tp, nrow(moon_scoring_result$moon_res), nrow(pruned$nodes), nrow(pruned$edges), nrow(pruned$gem_edges)
    ))
    pruned
}

## ---------------------------------------------------------------------
## 6.4 Loop over all 8 timepoints (T023) -- not x genotypes
## ---------------------------------------------------------------------

timepoints <- sort(unique(measured_features$timepoint_h[!is.na(measured_features$timepoint_h)]))
dir.create("result/moon", recursive = TRUE, showWarnings = FALSE)

moon_results <- list()
for (tp in timepoints) {
    result <- tryCatch(run_one_timepoint(tp), error = function(e) {
        cat(sprintf("timepoint %sh: FAILED (%s) -- logged, not fatal\n", tp, conditionMessage(e)))
        NULL
    })
    if (!is.null(result)) {
        moon_results[[as.character(tp)]] <- result
        saveRDS(result, sprintf("result/moon/%sh_moon_result.rds", tp))
    }
}

cat("\nMOON runs completed:", length(moon_results), "of", length(timepoints), "timepoints\n")
saveRDS(moon_results, "result/moon/all_timepoints.rds")
cat("Saved result/moon/{<timepoint>h_moon_result,all_timepoints}.rds\n")
