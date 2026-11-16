# Companion figure (requested 2026-10-08, after the ATP/WT_2h episode):
# the per-pair pruned-network figures (09/13) only show molecules that
# survived MOON's connectivity pruning -- that answers "is this molecule
# part of a reconstructed signaling story," not "did this molecule
# change." Confirmed directly with ATP: at 2h, WT's raw drop (5.17 ->
# 2.60 nmol/mg, t=-3.07) is LARGER than ob/ob's (3.41 -> 2.56, t=-1.85),
# yet ATP survives pruning in ob_2h and not WT_2h -- tested n_steps up to
# 12 (double the pipeline's 6) and the pruned network barely changes
# (354->355 nodes), so it isn't a "signal just needs a longer leash"
# issue; WT_2h's ATP genuinely has no sign-coherent path to an upstream
# seed under this reconstruction. A reader who only sees the pruned-
# network figure would wrongly conclude ATP "didn't move" in WT.
#
# This script renders, per timepoint, EVERY directly-measured metabolite's
# raw Tn-vs-T0 Welch's t (both genotypes), independent of MOON pruning --
# the ground truth the network figures are built from -- with filled vs.
# open points marking whether that measurement also survived into that
# genotype's pruned network at this timepoint. Not a replacement for the
# network figures; a companion so "excluded from the figure" is never
# silently read as "no change."
#
# Run from omnipath_metabo_case2/, after scripts/06_footprint_moon.R.

suppressMessages(library(ggplot2))

TIMEPOINTS <- c(2, 4, 6, 8, 12, 16, 24)
GT_COLORS <- c(WT = "#1F6F5C", ob = "#B5631F")  # same WT/ob convention as the boss-summary artifact

measured_features <- readRDS("Tn_T0/result/pk_retrieval/measured_features.rds")
moon_results <- readRDS("Tn_T0/result/moon/all_pairs.rds")

bare_chebi <- function(node_id) sub("_[a-z]+$", "", sub("^Metab__", "", node_id))

dir.create("Tn_T0/result/networks/focus", recursive = TRUE, showWarnings = FALSE)

for (tp in TIMEPOINTS) {

    mf <- measured_features[
        measured_features$omics_layer %in% c("metabolome", "plasma_metabolome") &
            measured_features$timepoint_h == tp & !measured_features$excluded &
            measured_features$mapping_status == "mapped" & !is.na(measured_features$t_stat),
    ]
    parts <- strsplit(mf$feature_id, ";", fixed = TRUE)
    mf$name <- vapply(parts, function(p) if (length(p) >= 2) p[2] else NA_character_, character(1))

    present_in <- function(gt_key, chebi_ids) {
        pruned <- moon_results[[gt_key]]
        node_chebi <- unique(bare_chebi(pruned$nodes$source[grepl("^Metab__", pruned$nodes$source)]))
        chebi_ids %in% node_chebi
    }

    wide <- reshape(
        mf[, c("name", "pk_node_id", "genotype", "t_stat")],
        idvar = c("name", "pk_node_id"), timevar = "genotype", direction = "wide"
    )
    names(wide) <- sub("^t_stat\\.", "", names(wide))
    wide <- wide[!is.na(wide$WT) | !is.na(wide$ob), ]
    wide$WT_in_network <- present_in(paste0("WT_", tp), wide$pk_node_id)
    wide$ob_in_network <- present_in(paste0("ob_", tp), wide$pk_node_id)

    # Disambiguate duplicate display names (same metabolite name, different
    # pk_node_id -- e.g. distinct compartments/ions resolving to one label)
    # the same way 12_cytoscape_export.R does elsewhere in this pipeline.
    dup <- wide$name %in% wide$name[duplicated(wide$name)]
    wide$label <- wide$name
    wide$label[dup] <- paste0(wide$name[dup], " (", wide$pk_node_id[dup], ")")

    wide$max_abs <- pmax(abs(wide$WT), abs(wide$ob), na.rm = TRUE)
    wide <- wide[order(-wide$max_abs), ]
    wide$label <- factor(wide$label, levels = rev(wide$label))

    long <- rbind(
        data.frame(label = wide$label, genotype = "WT", score = wide$WT, in_network = wide$WT_in_network),
        data.frame(label = wide$label, genotype = "ob", score = wide$ob, in_network = wide$ob_in_network)
    )
    long <- long[!is.na(long$score), ]

    segs <- wide[!is.na(wide$WT) & !is.na(wide$ob), ]

    p <- ggplot() +
        geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
        geom_segment(data = segs, aes(x = WT, xend = ob, y = label, yend = label), color = "grey80", linewidth = 0.5) +
        geom_point(data = long, aes(x = score, y = label, color = genotype, shape = in_network), size = 2.3, stroke = 1) +
        scale_color_manual(values = GT_COLORS, name = "Genotype") +
        scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "In pruned network?",
                            labels = c(`TRUE` = "yes", `FALSE` = "no (measured only)")) +
        labs(
            title = sprintf("Case study 2, %sh -- raw Tn-vs-T0 scores for every measured metabolite", tp),
            subtitle = "Ground truth independent of MOON pruning; filled = also survived into that genotype's pruned network",
            x = "Welch's t (Tn vs T0)", y = NULL
        ) +
        theme_minimal(base_size = 11) +
        theme(
            panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
            plot.background = element_rect(fill = "white", color = NA),
            axis.text.y = element_text(size = 8)
        )

    n_rows <- nlevels(wide$label)
    out_file <- sprintf("Tn_T0/result/networks/focus/raw_scores_%sh.pdf", tp)
    ggsave(out_file, p, width = 9, height = max(6, 0.17 * n_rows + 1.5), limitsize = FALSE)
    cat(sprintf("%sh: %d measured metabolites -> %s\n", tp, n_rows, out_file))
}

cat("\nAll raw-score companion figures saved to Tn_T0/result/networks/focus/\n")
