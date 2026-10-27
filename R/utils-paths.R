#' Path to a file in the repository
#'
#' Resolves a path relative to the repository root, found by walking up
#' from the working directory to the folder holding \code{manuscript.yaml}.
#' The pipeline runs from the repository root, where the result is the
#' plain relative path (so provenance records stay free of absolute
#' paths); from elsewhere, e.g. \code{tests/testthat/}, the root is
#' prepended.
#'
#' @param ... Path components below the repository root.
#'
#' @return Character scalar path.
#' @keywords internal
#' @noRd
repo_path <- function(...) {

    dir <- normalizePath(getwd(), mustWork = FALSE)
    depth <- 0L

    while (!file.exists(file.path(dir, "manuscript.yaml"))) {
        parent <- dirname(dir)
        if (identical(parent, dir)) {
            return(file.path(...))
        }
        dir <- parent
        depth <- depth + 1L
    }

    if (depth == 0L) {
        file.path(...)
    } else {
        file.path(dir, ...)
    }
}
