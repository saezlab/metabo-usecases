# Renders the GEM-instance-merged network (scripts/12_cytoscape_export.R's
# output) as a plain scattered/force-directed layout -- NOT the layered
# TF/mRNA/PPI/enzyme/metabolite structure from 09_network_viz.R -- with
# node labels and the same node color/shape scheme set up interactively in
# Cytoscape (metabolite = light green diamond, gene_protein = orange
# circle), plus a high-contrast edge-category palette. Orphan reaction
# placeholder nodes (is_orphan in the att table) are dropped entirely from
# this figure -- they remain in the Cytoscape export files themselves.
#
# Run from omnipath_metabo_case2/, after scripts/12_cytoscape_export.R.

suppressMessages(library(igraph))
suppressMessages(library(ggraph))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))

TIMEPOINTS <- c("0", "2", "4", "6", "8", "12", "16", "24")
# Denser timepoints (post leaf-trim: 2h=251, 12h=149, 16h=176, 24h=449
# nodes) -- gene/protein labels alone create unreadable clutter at that
# density; metabolite-only labels still convey the figure's point.
METABOLITE_LABELS_ONLY <- c("2", "12", "16", "24")

node_colors <- c(metabolite = "#ADDD8E", gene_protein = "#FEC44F")
# High-contrast categorical palette (RColorBrewer Set1, yellow dropped --
# too low-contrast on white): red/blue/green/purple/orange/brown are each
# maximally separated in hue, not just "6 distinct colors" -- important
# with many thin overlapping edges.
edge_colors <- c(
    enzyme_met = "#E41A1C", grn = "#377EB8", ppi = "#4DAF4A",
    allosteric = "#984EA3", transporters = "#FF7F00", receptors = "#A65628"
)

dir.create("WT_vs_ob/result/networks/cytoscape", recursive = TRUE, showWarnings = FALSE)

for (tp in TIMEPOINTS) {

    att <- read.csv(sprintf("WT_vs_ob/result/networks/cytoscape/%sh_att.csv", tp))
    edges <- read.csv(sprintf("WT_vs_ob/result/networks/cytoscape/%sh_edge_att.csv", tp))

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

    label_layer <- if (tp %in% METABOLITE_LABELS_ONLY) {
        geom_node_text(data = function(x) subset(x, node_type == "metabolite"), aes(label = name), repel = TRUE,
                        size = 3.6, max.overlaps = Inf, segment.size = 0.15, bg.color = "white", bg.r = 0.1)
    } else {
        geom_node_text(aes(label = name), repel = TRUE, size = 3.6,
                        max.overlaps = Inf, segment.size = 0.15, bg.color = "white", bg.r = 0.1)
    }

    p <- ggraph(gl) +
        geom_edge_link(aes(color = category), alpha = 0.6, width = 0.5) +
        geom_node_point(aes(shape = node_type, fill = node_type, size = degree), color = "black", stroke = 0.3) +
        label_layer +
        scale_edge_color_manual(values = edge_colors, name = "Edge category") +
        scale_shape_manual(values = c(metabolite = 23, gene_protein = 21), name = "Node type") +
        scale_fill_manual(values = node_colors, name = "Node type") +
        scale_size_continuous(range = c(2, 10), name = "Degree") +
        theme_void(base_size = 13) +
        theme(plot.background = element_rect(fill = "white", color = NA)) +
        labs(title = sprintf("Case study 2, %sh -- GEM-instance-merged network", tp))

    out_file <- sprintf("WT_vs_ob/result/networks/cytoscape/%sh_styled.pdf", tp)
    ggsave(out_file, p, width = 13, height = 11, limitsize = FALSE)
    cat(sprintf("%sh: %d nodes, %d edges -> %s\n", tp, vcount(g), ecount(g), out_file))
}

cat("\nAll styled PDFs saved to WT_vs_ob/result/networks/cytoscape/\n")
