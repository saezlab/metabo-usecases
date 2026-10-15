# Ad hoc follow-up figure (requested 2026-10-08): ob_2h shows a strong ATP
# (CHEBI:15422) signal in its pruned network (score -1.85, level0/measured)
# that is completely absent from WT_2h's pruned network -- confirmed by
# direct lookup against Tn_T0/result/moon/all_pairs.rds, not assumed.
#
# Node/edge construction follows 09_network_viz.R's pattern (category_of(),
# gem_nodes/connecting-edges), restricted to the GEM edges that touch an
# ATP compartment node, with GEM-only enzymes (no non-GEM/"upstream"
# connection) dropped -- fixed 2026-10-08 to match on bare UniProt, since
# GEM edges key the protein side as a reaction-instance id (Gene3__
# Q9QXG4) while ppi/grn/allosteric/transporters/receptors edges key the
# same protein by bare UniProt (Q9QXG4); a naive same-string check matches
# nothing.
#
# Rendering REDONE 2026-10-08 (second pass) to match the main per-pair
# figures' own visual language instead of 09's layered DAG: force-directed
# "stress" layout, GEM reaction-instance merging, and symbol/name labels,
# all reused verbatim from 12_cytoscape_export.R / 13_cytoscape_style_viz.R
# (collapse_gem_instances(), label_for(), the KEGG name-fallback cache).
# ATP's own compartment nodes stay distinct (Metab__CHEBI:15422_c, _m, _n,
# ... are different nodes, not merged) -- label_for() already keeps
# metabolite compartment suffixes separate; only enzyme GEM-instances get
# collapsed, same split 12 already makes.
#
# Node fill is now the MOON score (diverging scale) instead of node_type,
# so direction (up/down) is visible at a glance -- the thing this whole
# figure exists to answer for ATP. Node_type is still encoded via shape.
#
# Run from omnipath_metabo_case2/, after scripts/06_footprint_moon.R and
# scripts/12_cytoscape_export.R (reuses its KEGG name cache).

suppressMessages(library(igraph))
suppressMessages(library(ggraph))
suppressMessages(library(ggplot2))
suppressMessages(library(org.Mm.eg.db))

ATP_CHEBI <- "CHEBI:15422"
PAIR_KEY <- "ob_2"

pkn_edges <- readRDS("Tn_T0/result/pk_retrieval/pkn_edges.rds")
moon_results <- readRDS("Tn_T0/result/moon/all_pairs.rds")
measured_features <- readRDS("Tn_T0/result/pk_retrieval/measured_features.rds")
pruned <- moon_results[[PAIR_KEY]]

## ---------------------------------------------------------------------
## Same category_of() as 09/12.
## ---------------------------------------------------------------------

origin_key <- paste(pkn_edges$source, pkn_edges$target)
origin_lookup <- stats::setNames(pkn_edges$category, origin_key)
category_of <- function(source, target) {
    fwd <- origin_lookup[paste(source, target)]
    rev <- origin_lookup[paste(target, source)]
    ifelse(!is.na(fwd), fwd, rev)
}

## ---------------------------------------------------------------------
## GEM edges restricted to ATP's own compartment nodes (kept distinct).
## ---------------------------------------------------------------------

atp_nodes <- grep(ATP_CHEBI, unique(c(pkn_edges$source, pkn_edges$target)), value = TRUE, fixed = TRUE)
gem_edges <- pruned$gem_edges
gem_edges <- gem_edges[gem_edges$source %in% atp_nodes | gem_edges$target %in% atp_nodes, ]
gem_nodes <- unique(c(gem_edges$source, gem_edges$target))

## ---------------------------------------------------------------------
## Upstream-connection filter, matched on bare UniProt (see header note),
## fanned back out to every GEM-instance node for that protein.
## ---------------------------------------------------------------------

enzymes_in_gem <- setdiff(gem_nodes, atp_nodes)
bare_uniprot <- function(x) sub("_rev$", "", sub("^Gene[0-9]+__", "", x))
bare_of <- bare_uniprot(enzymes_in_gem)
bare_to_instances <- split(enzymes_in_gem, bare_of)

signed <- pruned$edges
signed$category <- category_of(signed$source, signed$target)
touch_set <- unique(c(atp_nodes, bare_of))
connecting_raw <- signed[signed$source %in% touch_set | signed$target %in% touch_set, ]

expand_side <- function(id) if (id %in% names(bare_to_instances)) bare_to_instances[[id]] else id
connecting <- do.call(rbind, Map(function(s, t, cat) {
    expand.grid(from = expand_side(s), to = expand_side(t), category = cat, stringsAsFactors = FALSE)
}, connecting_raw$source, connecting_raw$target, connecting_raw$category))
if (is.null(connecting)) connecting <- data.frame(from = character(), to = character(), category = character())

has_upstream <- bare_of %in% c(connecting_raw$source, connecting_raw$target)
drop_enzymes <- enzymes_in_gem[!has_upstream]
cat(sprintf("ATP-linked enzymes: %d total, %d with no upstream connection (dropped)\n",
            length(enzymes_in_gem), length(drop_enzymes)))

edges <- unique(rbind(
    data.frame(from = gem_edges$source, to = gem_edges$target, category = "enzyme_met", stringsAsFactors = FALSE),
    connecting
))
edges <- edges[!(edges$category == "enzyme_met" & (edges$from %in% drop_enzymes | edges$to %in% drop_enzymes)), ]

## ---------------------------------------------------------------------
## Label/score lookups -- same construction as 12_cytoscape_export.R
## (collapse_gem_instances(), label_for(), KEGG name-fallback cache).
## ---------------------------------------------------------------------

all_nodes <- unique(c(pkn_edges$source, pkn_edges$target))
uniprot_nodes <- all_nodes[grepl("^[A-Z0-9]+$", all_nodes)]
symbol_df <- AnnotationDbi::select(org.Mm.eg.db, keys = uniprot_nodes, keytype = "UNIPROT", columns = "SYMBOL")
symbol_df <- symbol_df[!is.na(symbol_df$SYMBOL) & !duplicated(symbol_df$UNIPROT), ]
uniprot_symbol_map <- stats::setNames(symbol_df$SYMBOL, symbol_df$UNIPROT)

met_rows <- measured_features[measured_features$omics_layer %in% c("metabolome", "plasma_metabolome") & !is.na(measured_features$pk_node_id), ]
met_parts <- strsplit(met_rows$feature_id, ";", fixed = TRUE)
met_name <- vapply(met_parts, function(p) if (length(p) >= 2) p[2] else NA_character_, character(1))
chebi_name_map <- stats::setNames(met_name, met_rows$pk_node_id)
chebi_name_map <- chebi_name_map[!duplicated(names(chebi_name_map))]

kegg_chebi_cache <- "Tn_T0/result/networks/cytoscape/kegg_chebi_name_cache.rds"
kegg_chebi_name_map <- readRDS(kegg_chebi_cache)  # built by 12_cytoscape_export.R

collapse_gem_instances <- function(node_id) {
    m <- regmatches(node_id, regexec("^Gene[0-9]+__([A-Z][A-Z0-9]+)(?:_rev)?$", node_id, perl = TRUE))[[1]]
    if (length(m) == 2) return(paste0("GEMenz__", m[2]))
    node_id
}

label_for <- function(node_id) {
    if (grepl("^Metab__", node_id)) {
        chebi <- sub("_[a-z]+$", "", sub("^Metab__", "", node_id))
        compartment <- sub("^.*_([a-z]+)$", "\\1", node_id)
        name <- if (chebi %in% names(chebi_name_map)) chebi_name_map[[chebi]]
                else if (chebi %in% names(kegg_chebi_name_map)) kegg_chebi_name_map[[chebi]]
                else chebi
        return(paste0(name, "_", compartment))
    }
    if (grepl("^GEMenz__", node_id)) {
        uniprot <- sub("^GEMenz__", "", node_id)
        return(if (uniprot %in% names(uniprot_symbol_map)) uniprot_symbol_map[[uniprot]] else uniprot)
    }
    if (node_id %in% names(uniprot_symbol_map)) return(uniprot_symbol_map[[node_id]])
    node_id
}

score_lookup <- stats::setNames(pruned$nodes$score, pruned$nodes$source)
score_for <- function(node_id) {
    if (node_id %in% names(score_lookup)) return(unname(score_lookup[[node_id]]))
    if (grepl("^GEMenz__", node_id)) {
        uniprot <- sub("^GEMenz__", "", node_id)
        if (uniprot %in% names(score_lookup)) return(unname(score_lookup[[uniprot]]))
    }
    NA_real_
}

## ---------------------------------------------------------------------
## Collapse GEM reaction-instance nodes to one per protein (12's own
## merge), then build the final node/edge tables.
## ---------------------------------------------------------------------

edges$from <- vapply(edges$from, collapse_gem_instances, character(1))
edges$to <- vapply(edges$to, collapse_gem_instances, character(1))
edges <- unique(edges)

node_names <- unique(c(edges$from, edges$to))
node_type <- ifelse(grepl("^Metab__", node_names), "metabolite", "gene_protein")
label <- vapply(node_names, label_for, character(1))
score <- vapply(node_names, score_for, numeric(1))
nodes <- data.frame(name = node_names, node_type = node_type, label = unname(label), score = unname(score), stringsAsFactors = FALSE)

cat(sprintf("Final figure: %d nodes, %d edges\n", nrow(nodes), nrow(edges)))

## ---------------------------------------------------------------------
## Render -- 13_cytoscape_style_viz.R's force-directed "stress" layout and
## shape scheme, node fill now the MOON score (diverging) instead of
## node_type, so direction is visible directly.
## ---------------------------------------------------------------------

category_colors <- c(
    enzyme_met = "#E41A1C", grn = "#377EB8", ppi = "#4DAF4A",
    allosteric = "#984EA3", transporters = "#FF7F00", receptors = "#A65628"
)

g <- graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
V(g)$degree <- igraph::degree(g, mode = "all")

set.seed(1)
gl <- create_layout(g, layout = "stress")
gl$x <- gl$x * 0.7
gl$y <- gl$y * 0.7

score_limit <- max(abs(nodes$score), na.rm = TRUE)

p <- ggraph(gl) +
    geom_edge_link(aes(color = category), alpha = 0.55, width = 0.4,
                    arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
                    end_cap = circle(3, "mm")) +
    geom_node_point(aes(shape = node_type, fill = score, size = degree), color = "black", stroke = 0.3) +
    geom_node_text(aes(label = label), repel = TRUE, size = 3.3, max.overlaps = Inf,
                    segment.size = 0.15, bg.color = "white", bg.r = 0.1) +
    scale_edge_color_manual(values = category_colors, name = "Edge category") +
    scale_shape_manual(values = c(metabolite = 23, gene_protein = 21), name = "Node type") +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
                          limits = c(-score_limit, score_limit), na.value = "grey80",
                          name = "MOON score\n(blue = down, red = up)") +
    scale_size_continuous(range = c(2.5, 10), name = "Degree") +
    guides(fill = guide_colorbar(order = 1)) +
    theme_void(base_size = 13) +
    theme(plot.background = element_rect(fill = "white", color = NA)) +
    labs(title = sprintf("Case study 2, %s -- ATP neighborhood", PAIR_KEY),
         subtitle = "Enzymes touching ATP with an upstream (non-GEM) connection only; fill = MOON score direction")

dir.create("Tn_T0/result/networks/focus", recursive = TRUE, showWarnings = FALSE)
out_file <- sprintf("Tn_T0/result/networks/focus/ATP_%s_enzyme_ppi.pdf", PAIR_KEY)
ggsave(out_file, p, width = 13, height = 11, limitsize = FALSE)
cat("\nSaved:", out_file, "\n")
