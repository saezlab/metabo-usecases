# Case-study-2-specific PK helpers (spec 002-case-study-2-network, T007/T008).
# Everything else (PKN loading, node-type classification, edge provenance,
# MOON scoring/pruning, PACON) reuses Spatial-COSMOS-MISTy functions directly
# -- see research.md R1-R7. These two are the only genuinely new pieces.

#' Records which local PKN snapshot files (and their state) a run used.
#'
#' FR-019: the pipeline is pinned to the local, already-exported snapshot
#' rather than a live query, so reproducibility across reruns depends on
#' knowing exactly which files (and whether their content changed) a given
#' run read. mtime alone isn't enough (a checkout can touch a file without
#' changing its content), so this also records an md5 of each file.
#'
#' @param pkn_dir Directory holding the `cosmos_pkn_<category>.csv` files
#'   (`Spatial-COSMOS-MISTy/data/PKN`).
#' @param categories Character vector of PKN categories to record.
#' @return A data.frame, one row per category: `category`, `pkn_snapshot_file`,
#'   `pkn_snapshot_mtime` (ISO 8601), `pkn_snapshot_md5`.
record_pkn_snapshot_provenance <- function(
    pkn_dir,
    categories = c("allosteric", "enzyme_met", "grn", "ppi", "receptors", "transporters")
) {

    files <- file.path(pkn_dir, sprintf("cosmos_pkn_%s.csv", categories))
    missing <- files[!file.exists(files)]
    if (length(missing) > 0L) {
        stop(
            sprintf("case-study-2: missing PKN snapshot file(s): %s", paste(missing, collapse = ", ")),
            call. = FALSE
        )
    }

    info <- file.info(files)
    data.frame(
        category = categories,
        pkn_snapshot_file = files,
        pkn_snapshot_mtime = format(info$mtime, "%Y-%m-%dT%H:%M:%S%z"),
        pkn_snapshot_md5 = unname(tools::md5sum(files)),
        stringsAsFactors = FALSE
    )

}


#' Splits a combined PKEdge table into GEM and non-GEM edges.
#'
#' FR-015 / research.md R3: GEM (enzyme-metabolite) edge sign is always +1
#' and doesn't follow the same up/downstream activation-sign rule as
#' GRN/PPI/kinase-substrate edges, so GEM edges must not be scored by
#' `moon()`'s sign-weighted propagation. Nothing in the reused
#' `pkn_to_meta_network()` carves them out -- this must happen before the
#' combined PKN reaches `run_moon_scoring()`.
#'
#' @param pkn_edges A PKEdge data.frame with a `category` column (e.g. from
#'   T009's per-category loading loop).
#' @return `list(non_gem, gem)`, both PKEdge data.frames.
split_gem_edges <- function(pkn_edges) {

    is_gem <- pkn_edges$category == "enzyme_met"
    list(
        non_gem = pkn_edges[!is_gem, , drop = FALSE],
        gem = pkn_edges[is_gem, , drop = FALSE]
    )

}


#' Reattaches GEM edges to a pruned MOON network, unconstrained by sign coherence.
#'
#' FR-015: GEM edges are retained on evidence/connectivity grounds only, not
#' evaluated against the sign-coherence rule. "Connectivity grounds" means
#' anchored to the mechanistic hypothesis network, not both-endpoints-already-
#' scored: a GEM edge's metabolite side (substrate/product) is rarely itself a
#' MOON-scored node (GEM nodes aren't propagated through moon()'s sign-weighted
#' scoring at all, per research.md R3), so requiring both endpoints already
#' kept would filter out nearly every GEM edge and defeat the point of
#' exempting them. A GEM edge is kept if *either* endpoint survived pruning,
#' then added as a separate `gem_edges` field rather than folded into `edges`
#' (data-model.md's PrunedMechanisticNetwork).
#'
#' @param pruned_network A [reduce_moon_network()] result (`list(nodes,
#'   edges, moon_res, pruned_pkn)`).
#' @param gem_edges The `gem` element of [split_gem_edges()]'s output.
#' @return `pruned_network` with a `gem_edges` field added.
reattach_gem_edges <- function(pruned_network, gem_edges) {

    # Anchor to the final, threshold-reduced node set (reduce_moon_network()'s
    # `nodes`, i.e. reduced$ATT) -- NOT `pruned_pkn` (the broad, unthresholded
    # sign-coherent candidate network moon() searched over). Using `pruned_pkn`
    # here would reattach GEM edges against every node in that much larger
    # candidate set, defeating the point of pruning down to a compact
    # mechanistic hypothesis in the first place.
    kept_nodes <- unique(pruned_network$nodes$source)
    relevant_gem <- gem_edges[
        gem_edges$source %in% kept_nodes | gem_edges$target %in% kept_nodes,
        , drop = FALSE
    ]

    pruned_network$gem_edges <- relevant_gem
    pruned_network

}
