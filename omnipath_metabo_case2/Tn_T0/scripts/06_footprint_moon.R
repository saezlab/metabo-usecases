# User Story 2 (spec 003-case-study-2-temporal-trajectory, FR-006-010):
# per-(genotype, timepoint) footprint-based TF/kinase activity inference +
# COSMOS-MOON pruning. 14 pairs (7 non-zero timepoints x 2 genotypes), each
# with its own within-genotype Tn-vs-T0 input signal from
# 04_pk_retrieval.R -- NOT 002's single ob/ob-vs-WT-per-timepoint signal
# (8 networks total). Everything below this point mirrors 002's own
# 06_footprint_moon.R (same footprint methodology, same GEM-exemption
# split, same cofactor exclusion, same MOON direction) except where a
# comment says otherwise -- see spec 003's research.md for why each
# change was made.
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
# Upstream/downstream direction CORRECTED 2026-10-07 (third pass): per
# Morita et al.'s own framing of the starvation response (signaling/PPI
# and TF activity driving the transcriptional/metabolic response, with
# metabolite levels as the downstream consequence/readout, not the
# initiating signal), upstream = TF activity + kinase/PPI activity
# (footprint-derived), downstream = metabolite t-stats. This is also
# cosmosR's own standard "RNA/activity + metabolomics" scenario (its
# internal compress_same_children() names the downstream-input argument
# metab_input, i.e. downstream=metabolite is the library's default
# assumption) -- the earlier metabolite-upstream framing was carried over
# from the spatial pilot's ligand-upstream design and does not apply here.
# Spatial-COSMOS-MISTy's run_moon_scoring() gained a metab_side= argument
# to support this (compartment-fans whichever side is actually metabolite;
# previously hard-coded to fan upstream_input only).
#
# Run from omnipath_metabo_case2/, after scripts/04_pk_retrieval.R and
# scripts/05_network_topology.R.

suppressMessages(pkgload::load_all("../../Spatial-COSMOS-MISTy"))
suppressMessages(library(decoupleR))
source("Tn_T0/scripts/lib/pk_helpers.R")

pkn_edges <- readRDS("Tn_T0/result/pk_retrieval/pkn_edges.rds")
measured_features <- readRDS("Tn_T0/result/pk_retrieval/measured_features.rds")

## ---------------------------------------------------------------------
## 6.0 Drop generic cofactor metabolites (added 2026-10-07, via
## result/moon/network_size_per_timepoint.csv investigation): CHEBI:24636
## (proton, H+) is not a measured feature -- it's a nuclear particle that
## appears as a byproduct/cofactor in 4,201 GEM (enzyme_met) and 495
## transporter edges in the PKN, orders of magnitude more promiscuous
## than any real metabolite (the next-largest GEM hub has ~190 edges).
## When MOON's own propagation happens to push its score over the
## primary/secondary threshold (it did, at 8h/16h, via its transporter
## edges), reattach_gem_edges()'s either-endpoint rule pulls in its
## entire GEM reaction list, ballooning reattached-GEM-edge counts
## (3,553-4,217 at 8h/16h vs 220-2,220 elsewhere) without representing
## any real biological signal. Dropped from the PKN before any MOON step
## -- not just from GEM reattachment -- so it also can't inflate
## transporter-mediated reachability or appear as a pruned node at all.
## Added 2026-10-08: CHEBI:17544 (bicarbonate, HCO3-) and CHEBI:29311
## (chlorine radical, Cl.) -- both generic inorganic species (146-219 PKN
## edges each) rather than metabolites informative for this study, same
## rationale as the proton. CHEBI:29311 in particular looks like it may be
## a GEM-resource annotation artifact (a free chlorine radical has no
## obvious role in general liver metabolism, unlike HCO3- which is a
## legitimate, common carboxylation co-substrate -- but removed anyway,
## per the same "not analytically informative for this study" call as the
## proton and HCO3-).
cofactor_metabolites <- c("CHEBI:24636", "CHEBI:17544", "CHEBI:29311")  # proton (H+), bicarbonate (HCO3-), chlorine radical (Cl.)
bare_chebi <- function(node_id) sub("_[a-z]+$", "", sub("^Metab__", "", node_id))
is_cofactor_edge <- bare_chebi(pkn_edges$source) %in% cofactor_metabolites |
    bare_chebi(pkn_edges$target) %in% cofactor_metabolites
cat("Dropping", sum(is_cofactor_edge), "edges touching excluded cofactor metabolite(s):",
    paste(cofactor_metabolites, collapse = ", "), "\n")
pkn_edges <- pkn_edges[!is_cofactor_edge, ]

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

# 14 (genotype, timepoint) pairs (spec 003 FR-001/007), replacing 002's 8
# bare timepoints. gt_pairs$key ("WT_2", ..., "ob_24") is this feature's
# version of what used to be a bare timepoint string everywhere downstream
# -- same role, just genotype-qualified (data-model.md GenotypeTimepointPair).
genotypes <- sort(unique(measured_features$genotype[!is.na(measured_features$genotype)]))
timepoints <- sort(unique(measured_features$timepoint_h[!is.na(measured_features$timepoint_h)]))
gt_pairs <- expand.grid(genotype = genotypes, timepoint_h = timepoints, stringsAsFactors = FALSE)
gt_pairs$key <- paste(gt_pairs$genotype, gt_pairs$timepoint_h, sep = "_")
cat("\n", nrow(gt_pairs), "(genotype, timepoint) pairs:", paste(gt_pairs$key, collapse = ", "), "\n")

## ---------------------------------------------------------------------
## 6.2 Footprint step (FR-011/FR-006): TF activity from transcript t-stats
## via the grn regulon, kinase activity from phosphosite t-stats via the
## ppi regulon -- one run_ulm() call per layer, across all 14
## (genotype, timepoint) pairs at once (decoupleR's mat argument takes a
## feature x condition matrix directly; each pair is one "condition"
## column here, same mechanism 002 used for its 8 bare-timepoint columns).
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
        agg <- stats::aggregate(t_stat ~ pk_node_id + genotype + timepoint_h, data = rows, FUN = mean)
        rows <- agg
    }
    rows$key <- paste(rows$genotype, rows$timepoint_h, sep = "_")
    feature_ids <- sort(unique(rows$pk_node_id))
    mat <- matrix(NA_real_, nrow = length(feature_ids), ncol = nrow(gt_pairs),
                  dimnames = list(feature_ids, gt_pairs$key))
    for (i in seq_len(nrow(rows))) {
        mat[rows$pk_node_id[i], rows$key[i]] <- rows$t_stat[i]
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
## 6.3 Build per-(genotype,timepoint) upstream/downstream inputs
## (spec 003 FR-006, mirrors 002's T020/research.md R5 unchanged otherwise)
## ---------------------------------------------------------------------

build_moon_inputs <- function(gt_key) {
    genotype <- sub("_.*", "", gt_key)
    tp <- as.numeric(sub(".*_", "", gt_key))
    pair_features <- measured_features[
        measured_features$genotype == genotype & measured_features$timepoint_h == tp &
            !measured_features$excluded & measured_features$mapping_status == "mapped" & !is.na(measured_features$t_stat),
    ]

    # Upstream: TF activity + kinase/PPI activity (footprint-derived, 6.2).
    # 335 proteins score as both a GRN source (TF) and a PPI source (kinase)
    # -- the PKN has one node per UniProt ID regardless of role, so an
    # overlap needs one scalar value, not two. Averaged rather than picking
    # one arbitrarily; both are genuine activity-footprint estimates of the
    # same underlying protein's signaling output.
    tf_tp <- tf_activity[tf_activity$condition == gt_key, c("source", "score")]
    kin_tp <- kinase_activity[kinase_activity$condition == gt_key, c("source", "score")]
    combined <- rbind(tf_tp, kin_tp)
    upstream_input <- tapply(combined$score, combined$source, mean)
    upstream_input <- stats::setNames(as.numeric(upstream_input), names(upstream_input))

    # Downstream: metabolite t-stats -- the measured functional readout of
    # the signaling/transcriptional response (Morita et al.'s own framing;
    # also cosmosR's standard scenario, see comment above). Compartment-
    # fanned by run_moon_scoring(metab_side="downstream") below, same as
    # the pilot fans metabolite upstream_input.
    met_rows <- pair_features[pair_features$omics_layer %in% c("metabolome", "plasma_metabolome"), ]
    downstream_input <- stats::setNames(met_rows$t_stat, met_rows$pk_node_id)
    downstream_input <- downstream_input[!duplicated(names(downstream_input))]

    list(upstream_input = upstream_input, downstream_input = downstream_input)
}

## ---------------------------------------------------------------------
## 6.4 Run MOON scoring + pruning, reattach GEM edges
## (spec 003 FR-007-010, mirrors 002's T021/T022 unchanged otherwise)
## ---------------------------------------------------------------------
#
# Threshold sweep (research.md R6, FR-009): raised primary/secondary in
# steps across all 14 pairs' moon_scoring_result until pairs started
# collapsing to zero nodes (same method 002 used), backed off to the
# highest fully-stable rung. Swept 2026-10-08 against this feature's own
# within-genotype score distribution -- NOT copied from 002's 1.5/1.0
# (that value was swept against a different, between-genotype score
# distribution, and actually fails here: confirmed 2/14 pairs, WT_12h and
# WT_16h, collapse to 0 nodes at 1.5/1.0 on this signal). Sweep result:
# primary=1.10/secondary=0.60 is the highest threshold where all 14 pairs
# stay non-empty (397-668 nodes); 1.15/0.65 drops WT_16h to 0.
PRIMARY_THRESH <- 1.10
SECONDARY_THRESH <- 0.60

run_one_pair <- function(gt_key) {
    genotype <- sub("_.*", "", gt_key)
    tp <- as.numeric(sub(".*_", "", gt_key))
    inputs <- build_moon_inputs(gt_key)
    cat(sprintf(
        "\n--- %s: %d upstream (TF+kinase activity) | %d downstream (metabolite) ---\n",
        gt_key, length(inputs$upstream_input), length(inputs$downstream_input)
    ))

    # grn=NULL: filter_incohrent_TF_target()'s coherence check merges TF
    # scores against downstream_input keyed by the TF's GRN-regulon targets
    # (genes) -- with downstream_input now metabolite-keyed, that merge is
    # structurally empty (no-op), so it's skipped rather than left as dead
    # weight. The TF-activity footprint already incorporates the grn
    # regulon (6.2's run_ulm() call); this filter would have re-checked
    # coherence against raw RNA, which we no longer feed in.
    moon_scoring_result <- run_moon_scoring(
        node_activities = inputs[c("upstream_input", "downstream_input")],
        pkn = pkn_for_moon, grn = NULL, metab_side = "downstream",
        n_steps = 6, statistic = "ulm", compartments = all_compartments
    )
    # level0_exempt=FALSE for the same reason 002 chose it: our "level 0"
    # values are Welch's t-statistics over n=5 replicates per group, a real
    # statistical signal, not a single-sample spatial-pilot-style reading.
    #
    # Thresholds (research.md R6, spec FR-009): re-swept fresh for this
    # feature's own within-genotype Tn-vs-T0 score distribution, NOT
    # inherited from 002's primary=1.5/secondary=1.0 (that value was itself
    # specific to the between-genotype score distribution, and 002's own
    # history shows these thresholds interact with network connectivity in
    # non-obvious, not-purely-magnitude-driven ways -- see 002's documented
    # ATP/24h episode). Sweep result: primary=PRIMARY_THRESH/
    # secondary=SECONDARY_THRESH (placeholder constants defined above this
    # function -- see the sweep note right below their definition).
    pruned <- reduce_moon_network(moon_scoring_result, primary_thresh = PRIMARY_THRESH, secondary_thresh = SECONDARY_THRESH, level0_exempt = FALSE)
    pruned <- reattach_gem_edges(pruned, gem_edges)
    pruned$genotype <- genotype
    pruned$timepoint_h <- tp

    cat(sprintf(
        "%s: %d scored nodes -> pruned to %d nodes, %d edges (+ %d GEM edges reattached)\n",
        gt_key, nrow(moon_scoring_result$moon_res), nrow(pruned$nodes), nrow(pruned$edges), nrow(pruned$gem_edges)
    ))
    pruned
}

## ---------------------------------------------------------------------
## 6.5 Loop over all 14 (genotype, timepoint) pairs (spec 003 FR-007)
## ---------------------------------------------------------------------

dir.create("Tn_T0/result/moon", recursive = TRUE, showWarnings = FALSE)

moon_results <- list()
for (gt_key in gt_pairs$key) {
    result <- tryCatch(run_one_pair(gt_key), error = function(e) {
        cat(sprintf("%s: FAILED (%s) -- logged, not fatal\n", gt_key, conditionMessage(e)))
        NULL
    })
    if (!is.null(result)) {
        moon_results[[gt_key]] <- result
        saveRDS(result, sprintf("Tn_T0/result/moon/%s_moon_result.rds", gt_key))
    }
}

cat("\nMOON runs completed:", length(moon_results), "of", nrow(gt_pairs), "(genotype, timepoint) pairs\n")
saveRDS(moon_results, "Tn_T0/result/moon/all_pairs.rds")
cat("Saved Tn_T0/result/moon/{<genotype>_<timepoint>h_moon_result,all_pairs}.rds\n")

## ---------------------------------------------------------------------
## 6.6 Save TF/kinase activity footprint as its own artifact (FR-011)
## ---------------------------------------------------------------------
#
# Long-form duplicate of what 6.3 already feeds into moon() as
# upstream_input (per-source, not per-protein-averaged) -- kept as a
# standalone artifact for direct TF/kinase activity inspection without
# needing to re-derive it from a moon_results pruned-network object.

footprint_activity <- rbind(
    data.frame(node_id = tf_activity$source, genotype = sub("_.*", "", tf_activity$condition),
               timepoint_h = as.numeric(sub(".*_", "", tf_activity$condition)),
               score = tf_activity$score, activity_type = "TF", stringsAsFactors = FALSE),
    data.frame(node_id = kinase_activity$source, genotype = sub("_.*", "", kinase_activity$condition),
               timepoint_h = as.numeric(sub(".*_", "", kinase_activity$condition)),
               score = kinase_activity$score, activity_type = "kinase", stringsAsFactors = FALSE)
)
saveRDS(footprint_activity, "Tn_T0/result/moon/footprint_activity_scores.rds")
cat("Saved Tn_T0/result/moon/footprint_activity_scores.rds:", nrow(footprint_activity), "rows (",
    length(unique(tf_activity$source)), "TFs x", nrow(gt_pairs), "pairs +",
    length(unique(kinase_activity$source)), "kinases x", nrow(gt_pairs), "pairs )\n")
