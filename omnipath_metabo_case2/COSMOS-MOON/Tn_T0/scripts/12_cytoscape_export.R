# Exports the GEM-centered pruned network (same network rendered by
# 09_network_viz.R) as Cytoscape-ready files: a .sif (network topology)
# and an _att.csv (node attribute table), for requested timepoints.
#
# .sif: standard 3-column Cytoscape format (source<TAB>interaction<TAB>
# target). `interaction` is the edge's sign (1 = activation, -1 =
# inhibition) -- GEM (enzyme_met) edges included, per the GEM-centered
# filtering already used for the figures (full GEM backbone + any other-
# category edge that touches a GEM node).
#
# _att.csv: one row per node -- score/level/type from MOON's own ATT
# table, plus node_type (metabolite/gene_protein), label (gene symbol or
# metabolite name), and layer (TF/mRNA/PPI/enzyme/metabolite, same
# classification 09_network_viz.R uses for the layered figure) so the
# same visual scheme can be rebuilt in Cytoscape (Style > node color by
# `layer`, shape by `node_type`).
#
# A separate _edge_att.csv carries edge `category` (enzyme_met/grn/ppi/
# allosteric/transporters/receptors) for edge-color styling -- SIF's own
# `interaction` column only has room for the sign, not both.
#
# Run from omnipath_metabo_case2/, after scripts/06_footprint_moon.R.

suppressMessages(library(org.Mm.eg.db))

TIMEPOINTS <- as.vector(outer(c("WT", "ob"), c("2", "4", "6", "8", "12", "16", "24"), paste, sep = "_"))

pkn_edges <- readRDS("Tn_T0/result/pk_retrieval/pkn_edges.rds")
moon_results <- readRDS("Tn_T0/result/moon/all_pairs.rds")
measured_features <- readRDS("Tn_T0/result/pk_retrieval/measured_features.rds")

## ---------------------------------------------------------------------
## Label/role lookups (same as 09_network_viz.R / 07_pathway_control.R)
## ---------------------------------------------------------------------

origin_key <- paste(pkn_edges$source, pkn_edges$target)
origin_lookup <- stats::setNames(pkn_edges$category, origin_key)
category_of <- function(source, target) {
    fwd <- origin_lookup[paste(source, target)]
    rev <- origin_lookup[paste(target, source)]
    ifelse(!is.na(fwd), fwd, rev)
}

grn_sources <- unique(pkn_edges$source[pkn_edges$category == "grn"])
grn_targets <- unique(pkn_edges$target[pkn_edges$category == "grn"])
ppi_nodes <- unique(c(pkn_edges$source[pkn_edges$category == "ppi"], pkn_edges$target[pkn_edges$category == "ppi"]))
node_layer <- function(node_names) {
    is_metab <- grepl("^Metab__", node_names)
    ifelse(
        is_metab, "metabolite",
        ifelse(node_names %in% grn_sources, "TF",
        ifelse(node_names %in% grn_targets, "mRNA",
        ifelse(node_names %in% ppi_nodes, "PPI/kinase",
        "enzyme")))
    )
}

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

# chebi_name_map only covers this study's ~141 directly-measured
# metabolites -- most network nodes are MOON-propagated/GEM-connected
# metabolites never measured here at all, so they'd otherwise stay as raw
# CHEBI ids (confirmed: up to 73% of a timepoint's metabolite nodes, e.g.
# 12h). Fall back to KEGG's own compound name (via the same compound<->
# chebi conversion build_cosmos_metabolite_pathway_sets() already uses)
# for anything chebi_name_map doesn't cover -- resolves ~22 of 30
# previously-unresolved ids across all 8 timepoints; a handful of generic
# ions/isotopes with no 1:1 KEGG compound counterpart stay as raw ChEBI.
kegg_chebi_cache <- "Tn_T0/result/networks/cytoscape/kegg_chebi_name_cache.rds"
if (file.exists(kegg_chebi_cache)) {
    kegg_chebi_name_map <- readRDS(kegg_chebi_cache)
} else {
    suppressMessages(library(KEGGREST))
    cpd_names <- keggList("compound")
    conv <- keggConv("chebi", "compound")
    chebi_id <- sub("^chebi:", "CHEBI:", unname(conv))
    cpd_id <- sub("^cpd:", "", names(conv))
    name <- sub(";.*$", "", cpd_names[cpd_id])  # first synonym only
    kegg_chebi_name_map <- stats::setNames(name, chebi_id)
    kegg_chebi_name_map <- kegg_chebi_name_map[!is.na(kegg_chebi_name_map) & !duplicated(names(kegg_chebi_name_map))]
    saveRDS(kegg_chebi_name_map, kegg_chebi_cache)
}
cat("Name coverage: measured panel", length(chebi_name_map), "| KEGG fallback", length(kegg_chebi_name_map), "\n")

# A handful of ids still miss both of the above -- not because they're
# unnamed, but because ChEBI records protonation-state variants (e.g.
# arachidonic acid vs. its conjugate base arachidonate) as *separate* ids,
# and KEGG's own compound<->chebi conversion only points at one of them.
# keggConv("chebi","cpd:C00219") -> chebi:15843 (arachidonic acid), but
# this PKN's GEM node uses CHEBI:32395 (arachidonate) -- confirmed via
# direct ChEBI OLS lookup, not guessed. Manual aliases for cases found
# this way, checked individually (not a general proton-state resolver).
manual_chebi_name_overrides <- c(
    "CHEBI:32395" = "Arachidonate"
)

#' Collapses every `Gene<N>__<UniProt>[_rev]` GEM reaction-instance id for
#' the same protein into one merged id. omnipath-metabo's GEM builder gives
#' a multi-step enzyme (e.g. Fasn/P19096, fatty acid synthase, which
#' catalyzes a whole sequence of condensation/reduction/dehydration steps)
#' a separate `Gene<N>__` node per reaction step it participates in -- ~48
#' distinct nodes for Fasn alone in the 0h network. That granularity is
#' real and meaningful for the underlying GEM/MOON analysis (kept as-is
#' everywhere else in the pipeline), but is just visual clutter for a
#' node-level Cytoscape view, where "this protein shows up as an enzyme
#' node" is what's useful, not which specific reaction instance. Forward
#' and `_rev` (reverse-direction) instances of the same protein are merged
#' together too. Deliberately does NOT merge into the bare-UniProt
#' (non-GEM, GRN/PPI signaling-layer) node for the same protein, if one
#' exists -- that node represents a different role (regulated gene /
#' signaling node vs. catalytic activity) and merging across them would
#' blur the TF/mRNA/PPI/enzyme/metabolite layering. `Gene<N>__orphanReac*`
#' placeholders (no real protein behind them) are left untouched.
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
                else if (chebi %in% names(manual_chebi_name_overrides)) manual_chebi_name_overrides[[chebi]]
                else chebi
        return(paste0(name, "_", compartment))  # compartment suffix kept -- different compartments are different nodes
    }
    if (grepl("^GEMenz__", node_id)) {  # post-collapse_gem_instances() merged id
        uniprot <- sub("^GEMenz__", "", node_id)
        return(if (uniprot %in% names(uniprot_symbol_map)) uniprot_symbol_map[[uniprot]] else uniprot)
    }
    # Gene<N>__<UniProt>[_rev] not yet collapsed (shouldn't normally reach
    # here -- kept as a fallback), or Gene<N>__orphanReac<id>[_rev]
    # (reaction placeholder, no real gene -- resolve only if there's a
    # UniProt to resolve).
    m <- regmatches(node_id, regexec("^Gene[0-9]+__([A-Z][A-Z0-9]+)(_rev)?$", node_id))[[1]]
    if (length(m) == 3) {
        uniprot <- m[2]
        rev_suffix <- m[3]
        symbol <- if (uniprot %in% names(uniprot_symbol_map)) uniprot_symbol_map[[uniprot]] else uniprot
        return(paste0(symbol, rev_suffix))
    }
    if (node_id %in% names(uniprot_symbol_map)) return(uniprot_symbol_map[[node_id]])
    node_id
}

#' Renames every node to a readable label, disambiguating collisions (two
#' different raw IDs resolving to the same label, e.g. two UniProt ids
#' annotated with the same gene symbol) by appending the raw id -- so
#' Cytoscape's SIF import never silently merges two distinct network nodes
#' that happen to share a display name.
unique_labels_for <- function(node_ids) {
    raw_labels <- vapply(node_ids, label_for, character(1))
    dup <- raw_labels %in% raw_labels[duplicated(raw_labels)]
    final <- raw_labels
    final[dup] <- paste0(raw_labels[dup], " (", node_ids[dup], ")")
    stats::setNames(final, node_ids)
}

## ---------------------------------------------------------------------
## Per-timepoint export
## ---------------------------------------------------------------------

dir.create("Tn_T0/result/networks/cytoscape", recursive = TRUE, showWarnings = FALSE)

for (tp in TIMEPOINTS) {

    pruned <- moon_results[[tp]]
    gem_edges <- pruned$gem_edges
    gem_nodes <- unique(c(gem_edges$source, gem_edges$target))

    signed <- pruned$edges
    signed$category <- category_of(signed$source, signed$target)
    connecting <- signed[signed$source %in% gem_nodes | signed$target %in% gem_nodes, ]

    edges <- unique(rbind(
        data.frame(source = gem_edges$source, target = gem_edges$target, interaction = gem_edges$mor,
                   category = "enzyme_met", stringsAsFactors = FALSE),
        data.frame(source = connecting$source, target = connecting$target, interaction = connecting$interaction,
                   category = connecting$category, stringsAsFactors = FALSE)
    ))

    # Merge GEM reaction-instance nodes sharing the same protein (see
    # collapse_gem_instances() doc) before computing node names/labels --
    # a multi-step enzyme's ~40+ instance nodes become one.
    edges$source <- vapply(edges$source, collapse_gem_instances, character(1))
    edges$target <- vapply(edges$target, collapse_gem_instances, character(1))
    edges <- unique(edges)

    node_names <- unique(c(edges$source, edges$target))
    node_type <- ifelse(grepl("^Metab__", node_names), "metabolite", "gene_protein")
    layer <- node_layer(node_names)
    name_map <- unique_labels_for(node_names)  # raw id -> readable, collision-safe label
    score_lookup <- stats::setNames(pruned$nodes$score, pruned$nodes$source)
    level_lookup <- stats::setNames(pruned$nodes$level, pruned$nodes$source)
    type_lookup <- stats::setNames(pruned$nodes$type, pruned$nodes$source)
    # GEM reaction-instance nodes are collapsed to "GEMenz__<UniProt>" above
    # (collapse_gem_instances()), but pruned$nodes$source never contains
    # that prefix -- GEM is excluded from MOON's own scoring graph (FR-014/
    # 015 sign exemption), so a collapsed enzyme id matches nothing in
    # score_lookup and silently comes back NA for every enzyme node (found
    # 2026-10-08 while adding score-based node coloring). Falls back to the
    # bare UniProt's own score, if that same protein is independently
    # scored elsewhere in the network (e.g. as a PPI/kinase node) --
    # otherwise stays NA (no signaling-layer score exists for it).
    score_for <- function(id) {
        if (id %in% names(score_lookup)) return(unname(score_lookup[[id]]))
        if (grepl("^GEMenz__", id)) {
            uniprot <- sub("^GEMenz__", "", id)
            if (uniprot %in% names(score_lookup)) return(unname(score_lookup[[uniprot]]))
        }
        NA_real_
    }
    # Reaction placeholder with no real protein behind it (Gene<N>__
    # orphanReac<id>[_rev]) -- not biologically meaningful on its own, just
    # a GEM-model bookkeeping node. Flagged so the viz script can mute it
    # (no label, 50% alpha) rather than cluttering the figure.
    is_orphan <- grepl("^Gene[0-9]+__orphanReac", node_names)

    att <- data.frame(
        node = name_map[node_names],   # readable, Cytoscape-facing id
        raw_id = node_names,           # COSMOS node id (post-GEM-instance-merge for enzymes; see collapse_gem_instances())
        node_type = node_type,
        layer = layer,
        is_orphan = is_orphan,
        score = vapply(node_names, score_for, numeric(1)),
        level = level_lookup[node_names],
        moon_type = type_lookup[node_names],  # upstream_input / level0 / other / NA (GEM-only, no MOON score)
        stringsAsFactors = FALSE
    )

    edges_named <- edges
    edges_named$source <- name_map[edges$source]
    edges_named$target <- name_map[edges$target]
    # Cytoscape's own internal edge key (what it names each edge after
    # importing the .sif) is the literal string "source (interaction)
    # target" -- e.g. "Acp3 (1) Adenosine_c". The edge table import step
    # needs a column matching that exact format to key on; three separate
    # source/interaction/target columns are not enough; Table Column for
    # Edge Key -> "edge_key" in the import dialog.
    edges_named$edge_key <- paste0(edges_named$source, " (", edges_named$interaction, ") ", edges_named$target)

    sif_file <- sprintf("Tn_T0/result/networks/cytoscape/%sh.sif", tp)
    write.table(edges_named[, c("source", "interaction", "target")], sif_file,
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

    att_file <- sprintf("Tn_T0/result/networks/cytoscape/%sh_att.csv", tp)
    write.csv(att, att_file, row.names = FALSE)

    edge_att_file <- sprintf("Tn_T0/result/networks/cytoscape/%sh_edge_att.csv", tp)
    write.csv(edges_named[, c("edge_key", "source", "target", "interaction", "category")], edge_att_file, row.names = FALSE)

    cat(sprintf("%sh: %d nodes, %d edges -> %s, %s, %s\n",
                tp, nrow(att), nrow(edges), sif_file, att_file, edge_att_file))
}

cat("\nIn Cytoscape: File > Import > Network from File (the .sif), then\n")
cat("File > Import > Table from File (the _att.csv, key column 'node') to\n")
cat("attach node attributes -- style node color by 'layer', shape by\n")
cat("'node_type'. Import the _edge_att.csv as Edge Table Columns, key\n")
cat("column 'edge_key' (matches Cytoscape's own 'source (interaction)\n")
cat("target' edge name format) -- for edge color by 'category'.\n")
