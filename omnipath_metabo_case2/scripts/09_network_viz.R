# User Story 3 visualization (spec 002-case-study-2-network): render each
# timepoint's pruned mechanistic network, GEM (enzyme-metabolite)
# centered -- GRN/allosteric/PPI/transporter edges are kept only where
# they directly touch a GEM node (enzyme or metabolite); any node with no
# direct connection to the GEM backbone is dropped. Decided 2026-10-07
# after the edge-category stats showed GEM is the single largest category
# at every timepoint (2,474-4,094 edges, bigger than GRN), while receptors
# are genuinely rare (1-5 edges) -- the uncentered viz_cosmos_network()
# attempt buried this structure under the much larger mRNA-expression
# cluster.
#
# Run from omnipath_metabo_case2/, after scripts/06_footprint_moon.R.

suppressMessages(pkgload::load_all("../../Spatial-COSMOS-MISTy"))
source("scripts/lib/pk_helpers.R")

PKN_DIR <- "../../Spatial-COSMOS-MISTy/data/PKN"
pkn_edges <- readRDS("result/pk_retrieval/pkn_edges.rds")
moon_results <- readRDS("result/moon/all_timepoints.rds")

# File-origin category lookup (which cosmos_pkn_<category>.csv an edge
# came from), independent of whether annotate_cosmos_network() later finds
# a resource/evidence match for it -- same lookup used for the edge-
# category-by-origin stats.
origin_key <- paste(pkn_edges$source, pkn_edges$target)
origin_lookup <- stats::setNames(pkn_edges$category, origin_key)
category_of <- function(source, target) {
    fwd <- origin_lookup[paste(source, target)]
    rev <- origin_lookup[paste(target, source)]
    ifelse(!is.na(fwd), fwd, rev)
}

dir.create("result/networks/viz", recursive = TRUE, showWarnings = FALSE)

for (tp in names(moon_results)) {

    pruned <- moon_results[[tp]]
    gem_edges <- pruned$gem_edges  # full, not capped -- GEM is the backbone now
    gem_nodes <- unique(c(gem_edges$source, gem_edges$target))

    signed <- pruned$edges
    signed$category <- category_of(signed$source, signed$target)
    non_gem_touching_gem <- signed[signed$source %in% gem_nodes | signed$target %in% gem_nodes, ]

    gem_for_annot <- gem_edges
    names(gem_for_annot)[names(gem_for_annot) == "mor"] <- "interaction"
    all_edges <- unique(rbind(
        non_gem_touching_gem[, c("source", "target", "interaction")],
        gem_for_annot[, c("source", "target", "interaction")]
    ))

    nodes_df <- data.frame(source = unique(c(all_edges$source, all_edges$target)), stringsAsFactors = FALSE)
    annotated <- annotate_cosmos_network(nodes_df, all_edges, pkn_unformat_dir = PKN_DIR)

    # Override the resource-matched edge_category (only ~36% coverage for
    # GEM edges -- the rest fall into "pending" and get misclassified into
    # classify_cosmos_node_roles()'s generic fallback bucket) with the
    # reliable file-origin category instead, mapped to the same vocabulary
    # cosmos_edge_category_map uses internally.
    origin_to_display <- c(
        enzyme_met = "metabolic enzyme", grn = "GRN", ppi = "PPI",
        allosteric = "allosteric regulation", transporters = "transport", receptors = "ligand-receptor"
    )
    edge_origin <- ifelse(
        paste(all_edges$source, all_edges$target) %in% paste(gem_for_annot$source, gem_for_annot$target),
        "enzyme_met", category_of(all_edges$source, all_edges$target)
    )
    annotated$edges$edge_category <- unname(origin_to_display[edge_origin])

    type_lookup <- stats::setNames(pruned$nodes$type, pruned$nodes$source)
    annotated$nodes$type <- type_lookup[annotated$nodes$source]
    annotated$nodes$type[is.na(annotated$nodes$type)] <- "other"

    cat(sprintf(
        "%sh: %d GEM-backbone nodes -> %d total nodes (+%d connecting non-GEM), %d edges (%d GEM + %d connecting)\n",
        tp, length(gem_nodes), nrow(annotated$nodes), nrow(annotated$nodes) - length(gem_nodes),
        nrow(annotated$edges), nrow(gem_for_annot), nrow(non_gem_touching_gem)
    ))

    result <- viz_cosmos_network(
        nodes = annotated$nodes, edges = annotated$edges,
        label_mode = "none",
        plot_name = sprintf("case_study_2_%sh_gem_centered", tp),
        save_plot = "png",
        path = "result/networks/viz",
        width = 16,
        print_plot = FALSE
    )
    cat("  saved:", result$saved_files, "\n")
}

cat("\nAll GEM-centered timepoint network plots saved to result/networks/viz/\n")
