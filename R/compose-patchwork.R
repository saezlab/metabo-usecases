#' Compose ggplot panels into a single figure via patchwork
#'
#' Applies capital-letter panel labels (FR-017 / FR-024), explicit
#' layout widths/heights for deterministic rendering (research.md R-10),
#' and the shared style typography. Use this for R-only composites
#' (every panel produced by ggplot in the same R session); for mixed
#' R+TikZ+manual composites use \code{\link{compose_mixed_source}}.
#'
#' @param panels A named list of ggplot objects; names become the panel
#'     labels' alphabetical ordering source (first panel = A, etc.).
#'     Pass an unnamed list for default \code{LETTERS}-ordering.
#' @param layout Optional patchwork \code{plot_layout} args as a list
#'     (e.g. \code{list(widths = c(1, 2), heights = c(1, 1))}).
#' @param tag_levels Character: passed to \code{patchwork::plot_annotation}.
#'     Default \code{'A'} = capital letters per FR-017.
#'
#' @return A patchwork composite object suitable for \code{ggsave}.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p1 <- ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_bw_metabo(89)
#' p2 <- ggplot(mtcars, aes(hp, mpg)) + geom_point() + theme_bw_metabo(89)
#' compose_patchwork(list(p1, p2))
#' }
#'
#' @importFrom patchwork wrap_plots plot_annotation plot_layout
#' @importFrom rlang abort
#' @importFrom logger log_info
#' @export
compose_patchwork <- function(panels,
                              layout     = NULL,
                              tag_levels = "A") {

    if (length(panels) == 0L) {
        rlang::abort("compose_patchwork() requires at least one panel")
    }

    composite <- patchwork::wrap_plots(panels) +
        patchwork::plot_annotation(tag_levels = tag_levels)

    if (!is.null(layout)) {
        composite <- composite + do.call(patchwork::plot_layout, layout)
    }

    logger::log_info("Composed {length(panels)} panels via patchwork")
    composite
}
