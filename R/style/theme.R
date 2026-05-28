#' The metabo.figures shared ggplot2 theme
#'
#' \code{theme_bw}-derived theme used by every pipeline plot. Reduced
#' gridlines (FR-017), sans-serif Arial/Helvetica-compatible family,
#' minimum point sizes enforced per FR-017a (≥ 6 pt for labels /
#' panel / legend / body text; ≥ 5 pt for tick labels). The
#' \code{width_mm} argument lets the theme adjust panel-level
#' line weights proportionally to the rendered width.
#'
#' Used together with \code{\link{readability_check}} (SC-008) which
#' enforces the same minimum point sizes at render time.
#'
#' @param width_mm Numeric: target physical figure width in mm.
#'     \code{89} for single-column, \code{180} for double-column per
#'     the Nature/Bioinformatics envelope (FR-017a).
#'
#' @return A ggplot2 theme object.
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) +
#'     geom_point() +
#'     theme_bw_metabo(width_mm = 89)
#'
#' @importFrom ggplot2 theme_bw theme element_text element_line element_blank
#' @importFrom ggplot2 element_rect rel
#' @export
theme_bw_metabo <- function(width_mm = 89) {

    sizes <- font_sizes()

    base_line <- ifelse(width_mm >= 180, 0.4, 0.3)
    grid_line <- ifelse(width_mm >= 180, 0.25, 0.2)

    ggplot2::theme_bw(base_size = sizes$body) +
        ggplot2::theme(
            text             = ggplot2::element_text(
                family = "sans", colour = "black"
            ),
            axis.title       = ggplot2::element_text(
                size = sizes$axis_label
            ),
            axis.text        = ggplot2::element_text(
                size = sizes$tick
            ),
            axis.ticks       = ggplot2::element_line(
                linewidth = base_line
            ),
            axis.line        = ggplot2::element_blank(),
            panel.border     = ggplot2::element_rect(
                fill = NA, colour = "black", linewidth = base_line
            ),
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_line(
                colour = "grey88", linewidth = grid_line
            ),
            plot.title       = ggplot2::element_text(
                size = sizes$panel_title, face = "bold"
            ),
            plot.subtitle    = ggplot2::element_text(
                size = sizes$body
            ),
            legend.title     = ggplot2::element_text(
                size = sizes$legend, face = "bold"
            ),
            legend.text      = ggplot2::element_text(
                size = sizes$legend
            ),
            strip.text       = ggplot2::element_text(
                size = sizes$panel_title, face = "bold"
            ),
            strip.background = ggplot2::element_rect(
                fill = "grey95", colour = NA
            )
        )
}


#' Minimum point sizes enforced by the style module
#'
#' Implements the FR-017a envelope. The values are the floor — the
#' readability gate (\code{\link{readability_check}}) fails any plot
#' that renders below these at the target figure width.
#'
#' @return A named list \code{axis_label}, \code{panel_title},
#'     \code{legend}, \code{tick}, \code{body}.
#'
#' @export
font_sizes <- function() {
    list(
        axis_label  = 6L,
        panel_title = 6L,
        legend      = 6L,
        body        = 6L,
        tick        = 5L
    )
}


#' Standard panel widths in millimetres
#'
#' Returns the Nature/Bioinformatics-family envelope from FR-017a.
#'
#' @return A named list \code{single_col = 89}, \code{double_col = 180}.
#'
#' @export
panel_widths_mm <- function() {
    list(
        single_col = 89L,
        double_col = 180L
    )
}
