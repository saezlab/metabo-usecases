#' Load the old COSMOS prior-knowledge network from cosmosR
#'
#' Installs \pkg{cosmosR} from \code{saezlab/cosmosR} on GitHub if not
#' already present, then loads the \code{meta_network} dataset bundled
#' with the package.  Returns edge counts per interaction-type category
#' (FR-011, FR-011a) suitable for the grouped-bar comparison renderer.
#'
#' The old PKN contains ONLY edges derived from metabolic reactions.
#' The comparison figure (FR-011a) must make the scope mismatch
#' explicit: categories present in the new PKN but absent here appear
#' with an empty old-PKN bar.
#'
#' @return A tibble with columns \code{interaction_type} and
#'     \code{n_edges}, ordered by \code{n_edges} descending.
#'     Carries a \code{"source_pkg"} attribute (package name + version)
#'     and a \code{"fingerprint"} attribute (package version string) for
#'     sidecar recording.
#'
#' @importFrom rlang abort
#' @importFrom dplyr count rename arrange desc
#' @importFrom tibble as_tibble tibble
#' @export
cosmos_old_pkn <- function() {

    if (!requireNamespace("devtools", quietly = TRUE)) {
        install.packages("devtools")
    }

    if (!requireNamespace("cosmosR", quietly = TRUE)) {
        devtools::install_github("saezlab/cosmosR")
    }

    env <- new.env(parent = emptyenv())
    utils::data("meta_network", package = "cosmosR", envir = env)

    if (!"meta_network" %in% ls(envir = env)) {
        rlang::abort(
            "cosmosR::meta_network dataset not found after loading the package."
        )
    }

    pkn <- tibble::as_tibble(get("meta_network", envir = env))

    if ("interaction_type" %in% names(pkn)) {
        counts <- pkn |>
            dplyr::count(interaction_type) |>
            dplyr::rename(n_edges = n) |>
            dplyr::arrange(dplyr::desc(n_edges))
    } else {
        # Older cosmosR builds carry integer sign only — all edges are
        # metabolic reactions by the scope of the old PKN.
        counts <- tibble::tibble(
            interaction_type = "metabolic_reactions",
            n_edges          = nrow(pkn)
        )
    }

    pkg_version <- as.character(utils::packageVersion("cosmosR"))
    attr(counts, "source_pkg")  <- paste0("cosmosR@", pkg_version)
    attr(counts, "fingerprint") <- pkg_version

    counts
}
