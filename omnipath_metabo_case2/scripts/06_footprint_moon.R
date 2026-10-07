# User Story 3 (spec 002-case-study-2-network, FR-011-015): per-timepoint
# footprint-based TF/kinase activity inference + COSMOS-MOON pruning.
#
# One network per timepoint, not per genotype (corrected 2026-10-07):
# FR-012's single ob/ob-vs-WT t-value per timepoint is the only MOON input
# signal -- genotype is already encoded inside it.
#
# Corrected 2026-10-07 (second pass): the original version fed raw
# transcript/protein/phosphosite t-stats directly into moon() as
# downstream_input, skipping FR-011's footprint-inference step entirely.
# Per cosmosR's moon_data_pkn_mapping_principles vignette, a total
# measurement is not automatically an activity score, and a phosphosite
# must not be force-mapped onto its parent protein's PKN node unless the
# PKN represents that exact site (it doesn't here -- confirmed, kinase-
# substrate layers are protein-level only). Added a real footprint step:
# decoupleR::run_ulm() against the grn regulon (TF activity from
# transcript t-stats) and against the ppi regulon (kinase activity from
# phosphosite t-stats, aggregated to parent-protein level first, since
# run_ulm's network argument needs node IDs matching the PKN's own
# granularity). Raw transcript/protein t-stats remain legitimate as-is
# for a gene's own node (vignette's "RNA target expression" / "functional-
# readout gate" uses); raw phosphosite values are no longer injected
# directly -- they're consumed only as decoupleR's input matrix.
#
# Upstream/downstream direction (metabolite=upstream, genes=downstream)
# follows FR-007's own chain (metabolite->receptor->kinase->TF->GRN->
# enzyme<->metabolite, i.e. metabolite-as-ligand driving transcription),
# not the cosmosR vignette's "RNA+metabolomics" example, whose opposite
# direction (TF-activity=upstream, metabolite=downstream) describes a
# different causal scenario (transcriptional control of metabolism).
#
# Run from omnipath_metabo_case2/, after scripts/04_pk_retrieval.R and
# scripts/05_network_topology.R.

suppressMessages(pkgload::load_all("../../Spatial-COSMOS-MISTy"))
suppressMessages(library(decoupleR))
source("scripts/lib/pk_helpers.R")

pkn_edges <- readRDS("result/pk_retrieval/pkn_edges.rds")
measured_features <- readRDS("result/pk_retrieval/measured_features.rds")

## ---------------------------------------------------------------------
## 6.1 GEM exemption (T021, FR-014/015): split before MOON ever sees it
## ---------------------------------------------------------------------

#' Collapses a (source, target, mor) edge table to one row per pair --
#' decoupleR::run_ulm()'s `network` argument errors on repeated edges.
#' Drops (source,target) pairs where contributing resources disagree on
#' sign (ambiguous evidence); collapses true duplicates to one row. The
#' main PKN graph fed to moon() tolerates parallel/conflicting edges fine
#' (igraph traversal, not a per-source regression design matrix) -- this
#' dedup is only for the regulon inputs to run_ulm().
dedup_regulon <- function(edges) {
    key <- paste(edges$source, edges$target)
    n_signs <- ave(edges$mor, key, FUN = function(x) length(unique(x)))
    edges <- edges[n_signs == 1, ]
    edges[!duplicated(paste(edges$source, edges$target)), ]
}

split <- split_gem_edges(pkn_edges)
pkn_for_moon <- split$non_gem[, c("source", "target", "mor")]
grn_for_moon <- pkn_edges[pkn_edges$category == "grn", c("source", "target", "mor")]
ppi_for_moon <- pkn_edges[pkn_edges$category == "ppi", c("source", "target", "mor")]
grn_regulon <- dedup_regulon(grn_for_moon)
ppi_regulon <- dedup_regulon(ppi_for_moon)
gem_edges <- split$gem[, c("source", "target", "mor")]
cat("PKN for MOON (GEM excluded):", nrow(pkn_for_moon), "edges | GEM (reattached post-pruning):", nrow(gem_edges), "edges\n")
cat("grn regulon:", nrow(grn_regulon), "of", nrow(grn_for_moon), "rows (conflicting-sign pairs dropped) |",
    "ppi regulon:", nrow(ppi_regulon), "of", nrow(ppi_for_moon), "rows\n")

all_compartments <- c("c", "e", "eg", "g", "i", "l", "m", "n", "r", "v", "x")
timepoints <- sort(unique(measured_features$timepoint_h[!is.na(measured_features$timepoint_h)]))

## ---------------------------------------------------------------------
## 6.2 Footprint step (FR-011): TF activity from transcript t-stats via the
## grn regulon, kinase activity from phosphosite t-stats via the ppi
## regulon -- one run_ulm() call per layer, across all timepoints at once
## (decoupleR's mat argument takes a feature x condition matrix directly).
## ---------------------------------------------------------------------

build_tstat_matrix <- function(omics_layer, aggregate_to_parent_protein = FALSE) {
    rows <- measured_features[
        measured_features$omics_layer == omics_layer & !measured_features$excluded &
            measured_features$mapping_status == "mapped" & !is.na(measured_features$t_stat),
    ]
    if (aggregate_to_parent_protein) {
        # phosphosite pk_node_id is already the parent protein's bare UniProt
        # (T012: phosphosites resolve to their parent protein, no site-level
        # PKN node exists) -- multiple phosphosites per protein average here.
        agg <- stats::aggregate(t_stat ~ pk_node_id + timepoint_h, data = rows, FUN = mean)
        rows <- agg
    }
    feature_ids <- sort(unique(rows$pk_node_id))
    mat <- matrix(NA_real_, nrow = length(feature_ids), ncol = length(timepoints),
                  dimnames = list(feature_ids, as.character(timepoints)))
    for (i in seq_len(nrow(rows))) {
        mat[rows$pk_node_id[i], as.character(rows$timepoint_h[i])] <- rows$t_stat[i]
    }
    mat
}

transcript_mat <- build_tstat_matrix("transcriptome")
phospho_protein_mat <- build_tstat_matrix("phosphoproteome", aggregate_to_parent_protein = TRUE)
cat("\nFootprint input matrices: transcript", nrow(transcript_mat), "genes |",
    "phosphosite-aggregated-to-protein", nrow(phospho_protein_mat), "proteins\n")

tf_activity <- run_ulm(transcript_mat, grn_regulon, minsize = 5)
kinase_activity <- run_ulm(phospho_protein_mat, ppi_regulon, minsize = 5)
cat("TF activity footprint:", length(unique(tf_activity$source)), "TFs scored\n")
cat("Kinase activity footprint:", length(unique(kinase_activity$source)), "kinases scored\n")

## ---------------------------------------------------------------------
## 6.3 Build per-timepoint upstream/downstream inputs (T020, research.md R5)
## ---------------------------------------------------------------------

build_moon_inputs <- function(tp) {
    tp_features <- measured_features[
        measured_features$timepoint_h == tp & !measured_features$excluded &
            measured_features$mapping_status == "mapped" & !is.na(measured_features$t_stat),
    ]
    upstream_rows <- tp_features[tp_features$omics_layer %in% c("metabolome", "plasma_metabolome"), ]
    upstream_input <- stats::setNames(upstream_rows$t_stat, upstream_rows$pk_node_id)
    upstream_input <- upstream_input[!duplicated(names(upstream_input))]

    # Raw transcript/protein t-stats only (legitimate as-is for a gene's
    # own node; vignette's "RNA target expression" / "functional-readout
    # gate" uses). TF/kinase activity scores (FR-011, computed above) are
    # NOT injected here: run_moon_scoring()'s internal
    # filter_incohrent_TF_target() is hard-coded to expect downstream_input
    # as measured GRN-target values for its coherence check -- mixing in
    # TF-activity-keyed entries (which are GRN *sources*, not targets)
    # collapsed the network to empty after a few pruning iterations
    # (confirmed 2026-10-07). TF/kinase activity is reported as its own
    # artifact instead (6.6 below), not forced into this wrapper's slot.
    raw_rows <- tp_features[tp_features$omics_layer %in% c("proteome", "transcriptome"), ]
    downstream_input <- stats::setNames(raw_rows$t_stat, raw_rows$pk_node_id)
    downstream_input <- downstream_input[!duplicated(names(downstream_input))]

    list(upstream_input = upstream_input, downstream_input = downstream_input)
}

## ---------------------------------------------------------------------
## 6.4 Run MOON scoring + pruning, reattach GEM edges (T021/T022)
## ---------------------------------------------------------------------

run_one_timepoint <- function(tp) {
    inputs <- build_moon_inputs(tp)
    cat(sprintf(
        "\n--- timepoint %sh: %d upstream (metabolite) | %d downstream (raw transcript+protein target) ---\n",
        tp, length(inputs$upstream_input), length(inputs$downstream_input)
    ))

    moon_scoring_result <- run_moon_scoring(
        node_activities = inputs[c("upstream_input", "downstream_input")],
        pkn = pkn_for_moon, grn = grn_for_moon,
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
## 6.5 Loop over all 8 timepoints (T023) -- not x genotypes
## ---------------------------------------------------------------------

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

## ---------------------------------------------------------------------
## 6.6 Save TF/kinase activity footprint as its own artifact (FR-011)
## ---------------------------------------------------------------------
#
# Reported independently rather than injected into moon() -- see 6.3's
# comment for why mixing it into downstream_input breaks
# filter_incohrent_TF_target(). Still satisfies FR-011 as a standalone
# per-timepoint footprint-inference deliverable.

footprint_activity <- rbind(
    data.frame(node_id = tf_activity$source, timepoint_h = as.numeric(tf_activity$condition),
               score = tf_activity$score, activity_type = "TF", stringsAsFactors = FALSE),
    data.frame(node_id = kinase_activity$source, timepoint_h = as.numeric(kinase_activity$condition),
               score = kinase_activity$score, activity_type = "kinase", stringsAsFactors = FALSE)
)
saveRDS(footprint_activity, "result/moon/footprint_activity_scores.rds")
cat("Saved result/moon/footprint_activity_scores.rds:", nrow(footprint_activity), "rows (",
    length(unique(tf_activity$source)), "TFs x", length(timepoints), "timepoints +",
    length(unique(kinase_activity$source)), "kinases x", length(timepoints), "timepoints )\n")
