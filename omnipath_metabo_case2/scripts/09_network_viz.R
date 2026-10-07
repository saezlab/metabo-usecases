# User Story 3 visualization (spec 002-case-study-2-network): render each
# timepoint's pruned mechanistic network, GEM (enzyme-metabolite) centered.
#
# Redesigned 2026-10-07 (third pass): layered top-to-bottom DAG following
# Morita et al.'s own Fig. 2B causal order -- TF -> mRNA (GRN target) ->
# PPI (kinase/signaling) -> enzyme -> metabolite -- confirmed directly
# from the paper (not the earlier force-directed/stress version, which
# had no explicit direction at all). Built fresh with igraph's Sugiyama
# layered-DAG layout (layout_with_sugiyama(), exposed in ggraph as
# layout="sugiyama"), NOT by reusing Spatial-COSMOS-MISTy's
# viz_cosmos_network() -- that function hard-codes the spatial pilot's own
# (wrong, for this study) role order and was explicitly ruled out earlier.
#
# A node's layer is its most-upstream role across the whole PKN (not just
# this timepoint's subnetwork): TF (GRN source) > mRNA (GRN target) > PPI
# (kinase/signaling node) > enzyme (GEM/transporter protein side) >
# metabolite. A protein can plausibly hold more than one role (e.g. a TF
# that is also a kinase substrate); the single most-upstream role wins,
# since a layered layout needs exactly one y-position per node.
#
# GEM is still the structural backbone for which nodes are kept: GRN/
# allosteric/PPI/transporter edges are kept only where they directly touch
# a GEM node (enzyme or metabolite); any node with no direct connection to
# the GEM backbone is dropped.
#
# Run from omnipath_metabo_case2/, after scripts/06_footprint_moon.R.

suppressMessages(library(igraph))
suppressMessages(library(ggraph))
suppressMessages(library(ggplot2))

pkn_edges <- readRDS("result/pk_retrieval/pkn_edges.rds")
moon_results <- readRDS("result/moon/all_timepoints.rds")

origin_key <- paste(pkn_edges$source, pkn_edges$target)
origin_lookup <- stats::setNames(pkn_edges$category, origin_key)
category_of <- function(source, target) {
    fwd <- origin_lookup[paste(source, target)]
    rev <- origin_lookup[paste(target, source)]
    ifelse(!is.na(fwd), fwd, rev)
}

## ---------------------------------------------------------------------
## Node role classification (TF > mRNA > PPI > enzyme > metabolite), from
## the full PKN's own edge membership -- independent of timepoint.
## ---------------------------------------------------------------------

grn_sources <- unique(pkn_edges$source[pkn_edges$category == "grn"])
grn_targets <- unique(pkn_edges$target[pkn_edges$category == "grn"])
ppi_nodes <- unique(c(pkn_edges$source[pkn_edges$category == "ppi"], pkn_edges$target[pkn_edges$category == "ppi"]))

node_layer <- function(node_names) {
    is_metab <- grepl("^Metab__", node_names)
    layer <- ifelse(
        is_metab, 5L,
        ifelse(node_names %in% grn_sources, 1L,
        ifelse(node_names %in% grn_targets, 2L,
        ifelse(node_names %in% ppi_nodes, 3L,
        4L)))  # enzyme/transporter protein fallback: in our GEM-centered
               # subnetwork, every non-metabolite node not otherwise
               # classified is there via a GEM or transporter edge.
    )
    layer
}

layer_labels <- c(`1` = "TF", `2` = "mRNA", `3` = "PPI/kinase", `4` = "enzyme", `5` = "metabolite")

category_colors <- c(
    enzyme_met = "#e41a1c", grn = "#ff7f00", ppi = "#4daf4a",
    allosteric = "#984ea3", transporters = "#377eb8", receptors = "#a65628"
)

dir.create("result/networks/viz", recursive = TRUE, showWarnings = FALSE)

for (tp in names(moon_results)) {

    pruned <- moon_results[[tp]]
    gem_edges <- pruned$gem_edges  # full, not capped -- GEM is the backbone
    gem_nodes <- unique(c(gem_edges$source, gem_edges$target))

    signed <- pruned$edges
    signed$category <- category_of(signed$source, signed$target)
    connecting <- signed[signed$source %in% gem_nodes | signed$target %in% gem_nodes, ]

    edges <- unique(rbind(
        data.frame(from = gem_edges$source, to = gem_edges$target, category = "enzyme_met", stringsAsFactors = FALSE),
        data.frame(from = connecting$source, to = connecting$target, category = connecting$category, stringsAsFactors = FALSE)
    ))

    node_names <- unique(c(edges$from, edges$to))
    node_type <- ifelse(grepl("^Metab__", node_names), "metabolite", "gene_protein")
    layer <- node_layer(node_names)
    nodes <- data.frame(name = node_names, node_type = node_type, layer = layer, stringsAsFactors = FALSE)

    g <- graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
    V(g)$degree <- igraph::degree(g, mode = "all")

    cat(sprintf("%sh: %d nodes, %d edges, computing sugiyama layered layout...\n", tp, nrow(nodes), nrow(edges)))
    set.seed(1)
    gl <- create_layout(g, layout = "sugiyama", layers = V(g)$layer)

    p <- ggraph(gl) +
        geom_edge_link(aes(color = category), alpha = 0.45, width = 0.35) +
        geom_node_point(aes(shape = node_type, size = degree), fill = "grey20", color = "black", alpha = 0.85) +
        scale_edge_color_manual(values = category_colors, name = "Edge category") +
        scale_shape_manual(values = c(metabolite = 23, gene_protein = 21), name = "Node type") +
        scale_size_continuous(range = c(1.5, 10), name = "Degree") +
        # sugiyama already places layer 1 (TF) at the highest y and layer 5
        # (metabolite) at the lowest -- i.e. y = 6 - layer -- so the axis
        # needs no reversal, just labels mapped back from y to layer name.
        scale_y_continuous(breaks = 1:5, labels = layer_labels[as.character(6 - (1:5))]) +
        theme_minimal(base_size = 14) +
        theme(
            panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
            axis.text.x = element_blank(), axis.title = element_blank(),
            plot.background = element_rect(fill = "white", color = NA)
        ) +
        labs(title = sprintf("Case study 2, %sh -- TF -> mRNA -> PPI -> enzyme -> metabolite (Morita et al. Fig. 2B order)", tp))

    out_file <- sprintf("result/networks/viz/case_study_2_%sh_layered.png", tp)
    ggsave(out_file, p, width = 18, height = 14, dpi = 150, limitsize = FALSE)
    cat("  saved:", out_file, "\n")
}

cat("\nAll layered DAG plots saved to result/networks/viz/\n")
