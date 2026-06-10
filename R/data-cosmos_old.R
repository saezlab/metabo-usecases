#' Load the old COSMOS prior-knowledge network
#'
#' Reads the vendored `inst/extdata/cosmos/meta_network.RData` —
#' the pre-OmniPath-Metabo COSMOS PKN serialised from
#' \code{cosmosR/data/meta_network.RData} (v1.18.1). Returns
#' edge counts per interaction-type category (FR-011, FR-011a)
#' suitable for direct use in the grouped-bar comparison renderer
#' in \code{\link{cosmos_old_vs_new}}.
#'
#' The old PKN contains ONLY edges derived from metabolic reactions.
#' The comparison figure (FR-011a) must make the scope mismatch
#' explicit: categories present in the new PKN but absent here appear
#' with an empty old-PKN bar.
#'
#' @return A tibble with columns \code{interaction_type} and
#'     \code{n_edges}, ordered by \code{n_edges} descending.
#'     Carries a \code{"source_path"} attribute (vendored file path)
#'     and a \code{"fingerprint"} attribute (MD5 hex digest) for
#'     sidecar recording.
#'
#' @importFrom rlang abort
#' @importFrom dplyr count rename arrange desc
#' @importFrom tibble as_tibble
#' @importFrom digest digest
#' @export
cosmos_old_pkn <- function() {

    rdata_path <- system.file(
        "extdata", "cosmos", "meta_network.RData",
        package = "metabo.figures"
    )

    if (!nzchar(rdata_path) || !file.exists(rdata_path)) {
        rlang::abort(paste0(
            "Vendored meta_network.RData not found. ",
            "Re-vendor per inst/extdata/cosmos/SOURCE.md."
        ))
    }

    env <- new.env(parent = emptyenv())
    load(rdata_path, envir = env)

    # The cosmosR object is named `meta_network`; it is a data.frame
    # with at minimum columns: source, interaction (sign), target.
    # An optional `interaction_type` column carries the category used
    # in FR-011a; fall back to inferring from edge count if absent.
    obj_names <- ls(envir = env)
    if (!"meta_network" %in% obj_names) {
        rlang::abort(sprintf(
            "meta_network.RData does not contain an object named meta_network. Found: %s",
            paste(obj_names, collapse = ", ")
        ))
    }

    pkn <- tibble::as_tibble(get("meta_network", envir = env))

    if ("interaction_type" %in% names(pkn)) {
        counts <- pkn |>
            dplyr::count(interaction_type) |>
            dplyr::rename(n_edges = n) |>
            dplyr::arrange(dplyr::desc(n_edges))
    } else {
        # Older cosmosR builds use integer sign only — all edges are
        # metabolic reactions by the scope of the old PKN.
        counts <- tibble::tibble(
            interaction_type = "metabolic_reactions",
            n_edges          = nrow(pkn)
        )
    }

    fingerprint <- digest::digest(file = rdata_path, algo = "md5")
    attr(counts, "source_path") <- rdata_path
    attr(counts, "fingerprint") <- fingerprint

    counts
}
