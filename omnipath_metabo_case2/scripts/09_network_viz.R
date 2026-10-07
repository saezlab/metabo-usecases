# User Story 3 visualization (spec 002-case-study-2-network): render each
# timepoint's pruned mechanistic network, GEM (enzyme-metabolite) centered.
#
# Redesigned 2026-10-07 (second pass): dropped Spatial-COSMOS-MISTy's
# viz_cosmos_network() entirely -- it hard-codes a strict linear pipeline
# order (Metabolite -> ... -> GRN -> mRNA expression, cosmos_network_
# role_levels) built for the spatial pilot's single-direction upstream-
# ligand/downstream-target framing. Our GEM-centered view isn't that: GEM
# is a hub multiple other categories (GRN, allosteric, PPI, transporter)
# attach to from different sides, not a one-way cascade. Standard
# force-directed/stress graph layout (igraph + ggraph), letting actual
# topology -- not a manually assigned role -- determine node position, is
# the right tool here. "stress" layout (vs "fr") was chosen after visual
# comparison: it separates individual hub-and-spoke structures into
# distinguishable petals rather than one overlapping mass.
#
# GEM is the structural backbone: GRN/allosteric/PPI/transporter edges are
# kept only where they directly touch a GEM node (enzyme or metabolite);
# any node with no direct connection to the GEM backbone is dropped.
#
# Run from omnipath_metabo_case2/, after scripts/06_footprint_moon.R.

suppressMessages(pkgload::load_all("../../Spatial-COSMOS-MISTy"))
suppressMessages(library(igraph))
suppressMessages(library(ggraph))
suppressMessages(library(ggplot2))
source("scripts/lib/pk_helpers.R")

pkn_edges <- readRDS("result/pk_retrieval/pkn_edges.rds")
moon_results <- readRDS("result/moon/all_timepoints.rds")

origin_key <- paste(pkn_edges$source, pkn_edges$target)
origin_lookup <- stats::setNames(pkn_edges$category, origin_key)
category_of <- function(source, target) {
    fwd <- origin_lookup[paste(source, target)]
    rev <- origin_lookup[paste(target, source)]
    ifelse(!is.na(fwd), fwd, rev)
}

category_colors <- c(
    enzyme_met = "#e31a1c", grn = "#ff7f00", ppi = "#b15928",
    allosteric = "#33a02c", transporters = "#1f78b4", receptors = "#6a3d9a"
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
    nodes <- data.frame(name = node_names, node_type = node_type, stringsAsFactors = FALSE)

    g <- graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
    V(g)$degree <- igraph::degree(g, mode = "all")

    cat(sprintf("%sh: %d nodes, %d edges, computing stress layout...\n", tp, nrow(nodes), nrow(edges)))
    set.seed(1)
    layout <- create_layout(g, layout = "stress")

    p <- ggraph(layout) +
        geom_edge_link(aes(color = category), alpha = 0.25, width = 0.3) +
        geom_node_point(aes(shape = node_type, size = degree), fill = "grey30", color = "black", alpha = 0.7) +
        scale_edge_color_manual(values = category_colors, name = "Edge category") +
        scale_shape_manual(values = c(metabolite = 23, gene_protein = 21), name = "Node type") +
        scale_size_continuous(range = c(0.5, 6), name = "Degree") +
        theme_void() +
        labs(title = sprintf("Case study 2, %sh -- GEM-centered subnetwork (stress layout)", tp)) +
        theme(plot.background = element_rect(fill = "white", color = NA))

    out_file <- sprintf("result/networks/viz/case_study_2_%sh_stress.png", tp)
    ggsave(out_file, p, width = 16, height = 16, dpi = 150, limitsize = FALSE)
    cat("  saved:", out_file, "\n")
}

cat("\nAll GEM-centered stress-layout plots saved to result/networks/viz/\n")
