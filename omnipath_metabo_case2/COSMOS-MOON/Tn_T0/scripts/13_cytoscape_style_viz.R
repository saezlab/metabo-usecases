# Renders the GEM-instance-merged network (scripts/12_cytoscape_export.R's
# output) as a plain scattered/force-directed layout -- NOT the layered
# TF/mRNA/PPI/enzyme/metabolite structure from 09_network_viz.R -- with
# node labels, shape by node_type (metabolite = diamond, gene_protein =
# circle), and a high-contrast edge-category palette. Orphan reaction
# placeholder nodes (is_orphan in the att table) are dropped entirely from
# this figure -- they remain in the Cytoscape export files themselves.
#
# Node FILL changed 2026-10-08 from node_type color to MOON score
# (diverging, blue = down / red = up, renormalized per pair to that
# network's own max |score|) -- requested so direction is visible on the
# main per-pair figures, not just the one-off ATP-focus figure. Enzyme
# nodes with no independent signaling-layer score (most GEM-only ones --
# see 12_cytoscape_export.R's score_for() note) render grey (NA), which is
# itself informative: it marks "this protein is in the GEM backbone but
# was never scored by MOON."
#
# Run from omnipath_metabo_case2/, after scripts/12_cytoscape_export.R.

suppressMessages(library(igraph))
suppressMessages(library(ggraph))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))

TIMEPOINTS <- as.vector(outer(c("WT", "ob"), c("2", "4", "6", "8", "12", "16", "24"), paste, sep = "_"))
# Label-density cutoff decided PROGRAMMATICALLY from each pair's own
# post-leaf-trim node count (research.md R5, spec FR-014) -- NOT a
# hardcoded per-timepoint list like 002's METABOLITE_LABELS_ONLY, which
# was tuned to that design's specific 8-timepoint node-count distribution
# and doesn't transfer to this feature's 14 pairs (likely a different
# distribution entirely, since the within-genotype score signal differs).
# Threshold chosen to match 002's own empirical boundary (full labels
# stayed readable up to ~60 nodes, metabolite-only kicked in by ~150).
DENSE_LABEL_THRESHOLD <- 200

# High-contrast categorical palette (RColorBrewer Set1, yellow dropped --
# too low-contrast on white): red/blue/green/purple/orange/brown are each
# maximally separated in hue, not just "6 distinct colors" -- important
# with many thin overlapping edges.
edge_colors <- c(
    enzyme_met = "#E41A1C", grn = "#377EB8", ppi = "#4DAF4A",
    allosteric = "#984EA3", transporters = "#FF7F00", receptors = "#A65628"
)

# Node fill switched 2026-10-08 from node_type (metabolite/gene_protein) to
# MOON score direction -- matches the ATP-focus figure
# (14_atp_focus_viz.R), requested so the same up/down-regulation reading
# applies to the main per-pair figures, not just that one-off. node_type
# is still encoded via shape (diamond/circle). Diverging scale, symmetric
# around 0 and renormalized PER FIGURE to that pair's own max |score| --
# a fixed cross-pair scale would wash out smaller-swing pairs.

dir.create("Tn_T0/result/networks/cytoscape", recursive = TRUE, showWarnings = FALSE)

for (tp in TIMEPOINTS) {

    att <- read.csv(sprintf("Tn_T0/result/networks/cytoscape/%sh_att.csv", tp))
    edges <- read.csv(sprintf("Tn_T0/result/networks/cytoscape/%sh_edge_att.csv", tp))

    # Orphan reaction placeholders (Gene<N>__orphanReac<id>[_rev], no real
    # protein behind them) are dropped entirely from this figure -- not
    # just muted -- per direct request. They remain in the Cytoscape
    # export files (0Xh_att.csv/.sif) for anyone exploring interactively.
    orphan_nodes <- att$node[att$is_orphan]
    edges <- edges[!(edges$source %in% orphan_nodes | edges$target %in% orphan_nodes), ]
    att <- att[!att$is_orphan, ]

    g <- graph_from_data_frame(edges[, c("source", "target")], directed = TRUE, vertices = att[, "node", drop = FALSE])
    idx <- match(V(g)$name, att$node)
    V(g)$node_type <- att$node_type[idx]
    V(g)$score <- att$score[idx]
    V(g)$degree <- igraph::degree(g, mode = "all")

    # Drop leaf gene/protein nodes (degree 1) -- a gene hanging off a
    # single hub with no further connections (e.g. a lone Cyp450 paralog,
    # an Acp2/Acp5/Acp6 dangling off one metabolite) adds clutter without
    # adding network structure. Metabolite leaves are kept -- they're the
    # figure's actual subject, not incidental. Single pass, not iterative:
    # a node that becomes degree-1 only after this cut stays, so this
    # doesn't cascade into stripping the network further than asked.
    leaf_genes <- V(g)$name[V(g)$node_type == "gene_protein" & V(g)$degree == 1]
    g <- delete_vertices(g, leaf_genes)
    V(g)$degree <- igraph::degree(g, mode = "all")  # recompute after the cut, for node sizing

    E(g)$category <- edges$category[match(
        paste(as_edgelist(g)[, 1], as_edgelist(g)[, 2]),
        paste(edges$source, edges$target)
    )]

    set.seed(1)
    gl <- create_layout(g, layout = "stress")
    # Compress stress's own coordinate spread (not the plot's physical
    # size) -- stress tends to leave a lot of empty space between
    # clusters; scaling coordinates toward the centroid tightens that
    # without changing relative layout/topology.
    gl$x <- gl$x * 0.6
    gl$y <- gl$y * 0.6

    label_layer <- if (vcount(g) > DENSE_LABEL_THRESHOLD) {
        geom_node_text(data = function(x) subset(x, node_type == "metabolite"), aes(label = name), repel = TRUE,
                        size = 3.6, max.overlaps = Inf, segment.size = 0.15, bg.color = "white", bg.r = 0.1)
    } else {
        geom_node_text(aes(label = name), repel = TRUE, size = 3.6,
                        max.overlaps = Inf, segment.size = 0.15, bg.color = "white", bg.r = 0.1)
    }

    score_limit <- max(abs(V(g)$score), na.rm = TRUE)

    p <- ggraph(gl) +
        geom_edge_link(aes(color = category), alpha = 0.6, width = 0.5) +
        geom_node_point(aes(shape = node_type, fill = score, size = degree), color = "black", stroke = 0.3) +
        label_layer +
        scale_edge_color_manual(values = edge_colors, name = "Edge category") +
        scale_shape_manual(values = c(metabolite = 23, gene_protein = 21), name = "Node type") +
        scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
                              limits = c(-score_limit, score_limit), na.value = "grey80",
                              name = "MOON score\n(blue = down, red = up)") +
        scale_size_continuous(range = c(2, 10), name = "Degree") +
        guides(fill = guide_colorbar(order = 1)) +
        theme_void(base_size = 13) +
        theme(plot.background = element_rect(fill = "white", color = NA)) +
        labs(title = sprintf("Case study 2, %sh -- GEM-instance-merged network", tp))

    out_file <- sprintf("Tn_T0/result/networks/cytoscape/%sh_styled.pdf", tp)
    ggsave(out_file, p, width = 13, height = 11, limitsize = FALSE)
    cat(sprintf("%sh: %d nodes, %d edges -> %s\n", tp, vcount(g), ecount(g), out_file))
}

cat("\nAll styled PDFs saved to Tn_T0/result/networks/cytoscape/\n")
