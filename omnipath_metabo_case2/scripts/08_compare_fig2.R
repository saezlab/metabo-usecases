# User Story 1 (spec 002-case-study-2-network, SC-001): node/edge-count
# comparison between our reconstructed per-genotype networks and the
# published Morita et al. network(s).
#
# Run from omnipath_metabo_case2/, after scripts/04_pk_retrieval.R.
#
# Attribution note (2026-10-07): only the supplementary-materials PDF is
# available in this project (paper/scisignal.ads2547_sm.pdf); its captions
# are "fig. S1"-"fig. S3" etc., not the main-text "Figure 2" spec.md's
# Intro names. Data File S5 ("the global starvation responsive transomic
# network") and Data File S6 ("the starvation responsive metabolic
# network") are the two candidate published networks with exact,
# genotype-specific edge-level ground truth (a WT/ob boolean per edge, not
# a visual estimate) -- both are reported below, clearly labeled, rather
# than guessing which one is literally "Figure 2" without the main text.

genotype_networks <- list(
    WT = readRDS("result/networks/WT_network.rds"),
    ob_ob = readRDS("result/networks/ob_ob_network.rds")
)

published_network_counts <- function(data_file, label) {
    edges <- as.data.frame(readxl::read_excel(
        sprintf("data/ads2547_data_file_%s.xlsx", data_file), sheet = "edge"
    ))
    wt_edges <- edges[edges$WT == TRUE, ]
    ob_edges <- edges[edges$ob == TRUE, ]
    data.frame(
        published_network = label,
        genotype = c("WT", "ob_ob"),
        published_node_count = c(
            length(unique(c(wt_edges$source_id, wt_edges$target_id))),
            length(unique(c(ob_edges$source_id, ob_edges$target_id)))
        ),
        published_edge_count = c(nrow(wt_edges), nrow(ob_edges)),
        stringsAsFactors = FALSE
    )
}

published <- rbind(
    published_network_counts("s5", "Data File S5 (global transomic network)"),
    published_network_counts("s6", "Data File S6 (metabolic network)")
)

ours <- data.frame(
    genotype = names(genotype_networks),
    our_node_count = vapply(genotype_networks, function(n) length(n$nodes), integer(1)),
    our_edge_count = vapply(genotype_networks, function(n) nrow(n$edges), integer(1)),
    stringsAsFactors = FALSE
)

comparison <- merge(published, ours, by = "genotype")
comparison$node_fold_increase <- round(comparison$our_node_count / comparison$published_node_count, 1)
comparison$edge_fold_increase <- round(comparison$our_edge_count / comparison$published_edge_count, 1)
comparison <- comparison[order(comparison$published_network, comparison$genotype), ]

cat("SC-001: our reconstructed network vs. the two candidate published networks\n\n")
print(comparison[, c("published_network", "genotype", "published_node_count", "our_node_count",
                      "node_fold_increase", "published_edge_count", "our_edge_count", "edge_fold_increase")])

dir.create("result", recursive = TRUE, showWarnings = FALSE)
write.csv(comparison, "result/fig2_comparison.csv", row.names = FALSE)
cat("\nSaved result/fig2_comparison.csv\n")
