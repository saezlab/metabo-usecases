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
source("WT_vs_ob/scripts/lib/pk_helpers.R")

pkn_edges <- readRDS("WT_vs_ob/result/pk_retrieval/pkn_edges.rds")
measured_features <- readRDS("WT_vs_ob/result/pk_retrieval/measured_features.rds")

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
##
## Expanded 2026-10-08, per colleague recommendation: adopted the seed ID
## list from saezlab/MetaProViz's get_exclusion_metabolites()
## (R/GetPriorKnowledge.R), which curates exactly this class of
## generic/currency species (its own "ions" and "atoms" classes -- the
## function's "small_molecules" class is KEGG-seeded; its ids below are
## the ones that also carry the "ions"/"small_molecules" tag after KEGG
## seeds were translated to CHEBI). Fetched the source directly
## (`gh api repos/saezlab/MetaProViz/contents/R/GetPriorKnowledge.R`), not
## from a summary. That function seeds 14+16 KEGG ids and ~467 native
## CHEBI ids (ions = children of CHEBI:24867, atoms = children of
## CHEBI:33250); the 27 unique KEGG seeds were translated to CHEBI via
## KEGGREST::keggConv() (same mechanism already used in 04/12 of this
## pipeline). Restricted to the subset that actually touches an edge in
## our own PKN (64 of ~486 candidates do; the rest are real CHEBI terms
## but never appear as PKN nodes here, so listing them would be a no-op)
## -- every one of those 64 was independently re-confirmed by name via a
## direct ChEBI OLS lookup (not trusted from MetaProViz's comments alone),
## e.g. water, dioxygen, CO2, H2O2, NH3/NH4+, nitric/nitrous
## oxide/nitrite, the generic "fatty acid" class term, "metal cation", and
## a long tail of bare metal atoms/ions (Na, K, Ca, Mg, Fe, Zn, Mn, Cu,
## etc.) that are GEM cofactor annotations, not measured study
## metabolites. Net effect: +6,960 edges newly dropped beyond the
## original 3-id list (12,021 total vs 5,061 before, out of 202,858).
cofactor_metabolites <- c(
    "CHEBI:24636",  # proton (H+)
    "CHEBI:17544",  # bicarbonate (HCO3-)
    "CHEBI:29311",  # chlorine radical (Cl.)
    "CHEBI:15377",  # water
    "CHEBI:15379",  # dioxygen
    "CHEBI:15858",  # bromide
    "CHEBI:16134",  # ammonia
    "CHEBI:16240",  # hydrogen peroxide
    "CHEBI:16301",  # nitrite
    "CHEBI:16382",  # iodide
    "CHEBI:16480",  # nitric oxide
    "CHEBI:16526",  # carbon dioxide
    "CHEBI:16793",  # mercury(2+)
    "CHEBI:17051",  # fluoride
    "CHEBI:17996",  # chloride
    "CHEBI:18248",  # iron atom
    "CHEBI:18276",  # dihydrogen
    "CHEBI:18291",  # manganese atom
    "CHEBI:18361",  # diphosphate(4-)
    "CHEBI:18420",  # magnesium(2+)
    "CHEBI:18421",  # superoxide
    "CHEBI:22984",  # calcium atom
    "CHEBI:25213",  # metal cation
    "CHEBI:25567",  # nitrous acid
    "CHEBI:25805",  # oxygen atom
    "CHEBI:26216",  # potassium atom
    "CHEBI:26708",  # sodium atom
    "CHEBI:27363",  # zinc atom
    "CHEBI:27563",  # arsenic atom
    "CHEBI:27568",  # selenium atom
    "CHEBI:27594",  # carbon atom
    "CHEBI:27698",  # vanadium atom
    "CHEBI:27998",  # tungsten atom
    "CHEBI:28938",  # ammonium
    "CHEBI:29033",  # iron(2+)
    "CHEBI:29034",  # iron(3+)
    "CHEBI:29035",  # manganese(2+)
    "CHEBI:29036",  # copper(2+)
    "CHEBI:29101",  # sodium(1+)
    "CHEBI:29103",  # potassium(1+)
    "CHEBI:29105",  # zinc(2+)
    "CHEBI:29108",  # calcium(2+)
    "CHEBI:29191",  # hydroxyl
    "CHEBI:29287",  # gold atom
    "CHEBI:30452",  # tellurium atom
    "CHEBI:30502",  # beryllium(2+)
    "CHEBI:30514",  # caesium atom
    "CHEBI:33336",  # lanthanum atom
    "CHEBI:33341",  # titanium atom
    "CHEBI:33364",  # platinum atom
    "CHEBI:33375",  # gadolinium atom
    "CHEBI:35104",  # strontium(2+)
    "CHEBI:35366",  # fatty acid (generic class term)
    "CHEBI:37136",  # barium(2+)
    "CHEBI:39124",  # calcium ion
    "CHEBI:48775",  # cadmium(2+)
    "CHEBI:48828",  # cobalt(2+)
    "CHEBI:49470",  # aluminium(3+)
    "CHEBI:49666",  # iridium atom
    "CHEBI:49701",  # lanthanum(3+)
    "CHEBI:49713",  # lithium(1+)
    "CHEBI:49786",  # nickel(2+)
    "CHEBI:49807",  # lead(2+)
    "CHEBI:49847",  # rubidium(1+)
    "CHEBI:84043"   # gallium(3+)
)
bare_chebi <- function(node_id) sub("_[a-z]+$", "", sub("^Metab__", "", node_id))
is_cofactor_edge <- bare_chebi(pkn_edges$source) %in% cofactor_metabolites |
    bare_chebi(pkn_edges$target) %in% cofactor_metabolites
cat("Dropping", sum(is_cofactor_edge), "edges touching", length(cofactor_metabolites),
    "excluded cofactor metabolites (3 original + 62 from MetaProViz's",
    "get_exclusion_metabolites() that touch this PKN)\n")
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

    # Upstream: TF activity + kinase/PPI activity (footprint-derived, 6.2).
    # 335 proteins score as both a GRN source (TF) and a PPI source (kinase)
    # -- the PKN has one node per UniProt ID regardless of role, so an
    # overlap needs one scalar value, not two. Averaged rather than picking
    # one arbitrarily; both are genuine activity-footprint estimates of the
    # same underlying protein's signaling output.
    tf_tp <- tf_activity[tf_activity$condition == as.character(tp), c("source", "score")]
    kin_tp <- kinase_activity[kinase_activity$condition == as.character(tp), c("source", "score")]
    combined <- rbind(tf_tp, kin_tp)
    upstream_input <- tapply(combined$score, combined$source, mean)
    upstream_input <- stats::setNames(as.numeric(upstream_input), names(upstream_input))

    # Downstream: metabolite t-stats -- the measured functional readout of
    # the signaling/transcriptional response (Morita et al.'s own framing;
    # also cosmosR's standard scenario, see comment above). Compartment-
    # fanned by run_moon_scoring(metab_side="downstream") below, same as
    # the pilot fans metabolite upstream_input.
    met_rows <- tp_features[tp_features$omics_layer %in% c("metabolome", "plasma_metabolome"), ]
    downstream_input <- stats::setNames(met_rows$t_stat, met_rows$pk_node_id)
    downstream_input <- downstream_input[!duplicated(names(downstream_input))]

    list(upstream_input = upstream_input, downstream_input = downstream_input)
}

## ---------------------------------------------------------------------
## 6.4 Run MOON scoring + pruning, reattach GEM edges (T021/T022)
## ---------------------------------------------------------------------

run_one_timepoint <- function(tp) {
    inputs <- build_moon_inputs(tp)
    cat(sprintf(
        "\n--- timepoint %sh: %d upstream (TF+kinase activity) | %d downstream (metabolite) ---\n",
        tp, length(inputs$upstream_input), length(inputs$downstream_input)
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
    # level0_exempt=TRUE (reduce_moon_network()'s default) is designed for
    # the spatial pilot's single-sample-per-spot data, where a measured
    # value has no replication/variance behind it and so is taken at face
    # value regardless of magnitude. Our "level 0" values are Welch's
    # t-statistics over n=5 replicates per group -- already a real
    # statistical signal, not a raw single reading -- so level0_exempt=FALSE
    # applies normal thresholding uniformly (decided 2026-10-07).
    # primary=1.5/secondary=1.0 (re-swept 2026-10-07 after the upstream/
    # downstream direction correction). Briefly raised to 1.75/1.25
    # (2026-10-08), then reverted back down (2026-10-08, same day) after
    # confirming it was the direct cause of ATP (CHEBI:15422) dropping out
    # of the 24h network -- ATP's raw footprint-independent score stayed
    # 8.81 throughout (way past either threshold), but at 1.75/1.25 it had
    # no sign-coherent path left to any upstream seed within n_steps=6, so
    # reduce_moon_network()'s connectivity requirement cut it entirely
    # despite the score. 1.5/1.0 restores that path. (Separately, cofactor
    # exclusion below also affects 24h's connectivity -- see that block's
    # own comment; the two issues are independent, diagnosed separately.)
    # Re-swept AGAIN 2026-10-08 after expanding the cofactor exclusion list
    # (65 ids, up from 3): 1.5/1.0 now collapses 12h to 0 nodes -- removing
    # 6,960 more edges changes connectivity, not just score magnitude, same
    # lesson as every prior resweep in this project. Checkpointed
    # run_moon_scoring() per timepoint and swept reduce_moon_network() alone
    # (cheap) at the same 1.5:1.0 ratio: 1.36/0.91 is the exact boundary
    # (12h first goes to 0 at 1.37); settled on 1.35/0.90 for a small
    # safety margin, same round-number style as the original sweep.
    pruned <- reduce_moon_network(moon_scoring_result, primary_thresh = 1.35, secondary_thresh = 0.90, level0_exempt = FALSE)
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

dir.create("WT_vs_ob/result/moon", recursive = TRUE, showWarnings = FALSE)

moon_results <- list()
for (tp in timepoints) {
    result <- tryCatch(run_one_timepoint(tp), error = function(e) {
        cat(sprintf("timepoint %sh: FAILED (%s) -- logged, not fatal\n", tp, conditionMessage(e)))
        NULL
    })
    if (!is.null(result)) {
        moon_results[[as.character(tp)]] <- result
        saveRDS(result, sprintf("WT_vs_ob/result/moon/%sh_moon_result.rds", tp))
    }
}

cat("\nMOON runs completed:", length(moon_results), "of", length(timepoints), "timepoints\n")
saveRDS(moon_results, "WT_vs_ob/result/moon/all_timepoints.rds")
cat("Saved result/moon/{<timepoint>h_moon_result,all_timepoints}.rds\n")

## ---------------------------------------------------------------------
## 6.6 Save TF/kinase activity footprint as its own artifact (FR-011)
## ---------------------------------------------------------------------
#
# Long-form duplicate of what 6.3 already feeds into moon() as
# upstream_input (per-source, not per-protein-averaged) -- kept as a
# standalone artifact for direct TF/kinase activity inspection without
# needing to re-derive it from a moon_results pruned-network object.

footprint_activity <- rbind(
    data.frame(node_id = tf_activity$source, timepoint_h = as.numeric(tf_activity$condition),
               score = tf_activity$score, activity_type = "TF", stringsAsFactors = FALSE),
    data.frame(node_id = kinase_activity$source, timepoint_h = as.numeric(kinase_activity$condition),
               score = kinase_activity$score, activity_type = "kinase", stringsAsFactors = FALSE)
)
saveRDS(footprint_activity, "WT_vs_ob/result/moon/footprint_activity_scores.rds")
cat("Saved result/moon/footprint_activity_scores.rds:", nrow(footprint_activity), "rows (",
    length(unique(tf_activity$source)), "TFs x", length(timepoints), "timepoints +",
    length(unique(kinase_activity$source)), "kinases x", length(timepoints), "timepoints )\n")
