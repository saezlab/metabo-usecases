# User Story 3 visualization (spec 002-case-study-2-network): render each
# timepoint's pruned mechanistic network, clustered by functional role
# (Metabolite / Transporter-Receptor / Transport / GEM / PPI / GRN / mRNA
# expression), via Spatial-COSMOS-MISTy's viz_cosmos_network().
#
# Run from omnipath_metabo_case2/, after scripts/06_footprint_moon.R.

suppressMessages(pkgload::load_all("../../Spatial-COSMOS-MISTy"))
source("scripts/lib/pk_helpers.R")

PKN_DIR <- "../../Spatial-COSMOS-MISTy/data/PKN"
moon_results <- readRDS("result/moon/all_timepoints.rds")

dir.create("result/networks/viz", recursive = TRUE, showWarnings = FALSE)

for (tp in names(moon_results)) {

    pruned <- moon_results[[tp]]
    kept_nodes <- pruned$nodes$source

    # GEM edges: cap to one representative edge per already-kept node
    # (2026-10-07) -- "either endpoint kept" (reattach_gem_edges()'s own
    # rule, correct for the analysis artifact) brings in thousands of new
    # context nodes and balloons the plot; "both endpoints kept" gives zero
    # GEM edges (metabolite endpoints are almost never already in the
    # signaling-backbone set). Capping shows the GEM relationship class
    # without either extreme.
    gem_touching <- pruned$gem_edges[pruned$gem_edges$source %in% kept_nodes | pruned$gem_edges$target %in% kept_nodes, ]
    gem_touching$anchor <- ifelse(gem_touching$source %in% kept_nodes, gem_touching$source, gem_touching$target)
    gem_for_viz <- gem_touching[!duplicated(gem_touching$anchor), c("source", "target", "mor")]
    names(gem_for_viz)[names(gem_for_viz) == "mor"] <- "interaction"

    all_edges <- unique(rbind(pruned$edges[, c("source", "target", "interaction")], gem_for_viz[, c("source", "target", "interaction")]))
    nodes_df <- data.frame(source = unique(c(all_edges$source, all_edges$target)), stringsAsFactors = FALSE)
    annotated <- annotate_cosmos_network(nodes_df, all_edges, pkn_unformat_dir = PKN_DIR)

    # Merge MOON's own type (upstream_input/level0/other) onto the
    # annotated nodes -- classify_cosmos_node_roles() reads this column.
    type_lookup <- stats::setNames(pruned$nodes$type, pruned$nodes$source)
    annotated$nodes$type <- type_lookup[annotated$nodes$source]
    annotated$nodes$type[is.na(annotated$nodes$type)] <- "other"

    cat(sprintf("%sh: %d nodes (%d GEM edges capped-in), rendering...\n", tp, nrow(annotated$nodes), nrow(gem_for_viz)))

    result <- viz_cosmos_network(
        nodes = annotated$nodes, edges = annotated$edges,
        label_mode = "none",
        plot_name = sprintf("case_study_2_%sh_network", tp),
        save_plot = "png",
        path = "result/networks/viz",
        width = 16,
        print_plot = FALSE
    )
    cat("  saved:", result$saved_files, "\n")
}

cat("\nAll timepoint network plots saved to result/networks/viz/\n")
