#' Render the FR-007b coverage profile panel
#'
#' Combined line graph: x = number of resources, y = number of items
#' supported by at least that many resources. One coloured line per
#' variant (Entities, Molecular entities, Structures). Log10 y because
#' the count drops by 4-5 orders of magnitude from N = 1 to N = 22.
#'
#' @param data Tibble from \code{\link{coverage_profile}} called with
#'     \code{variant = "all"}.
#' @param width_mm Numeric.
#'
#' @return A ggplot.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs theme
#' @importFrom ggplot2 element_text scale_y_log10 scale_x_continuous
#' @importFrom ggplot2 scale_colour_manual expansion
#' @importFrom rlang .data
#' @export
plot_coverage_profile <- function(data, width_mm = 180L) {

    variant_levels <- c("entities", "molecular_entities", "structures")
    variant_labels <- c("Entities", "Molecular entities", "Structures")
    names(variant_labels) <- variant_levels

    data$variant <- factor(data$variant, levels = variant_levels)

    cols <- stats::setNames(
        palette_n(length(variant_levels), unknown = FALSE),
        variant_levels
    )

    ggplot2::ggplot(
        data,
        ggplot2::aes(
            x      = .data$n_resources,
            y      = .data$n_items_ge_n,
            colour = .data$variant,
            group  = .data$variant
        )
    ) +
        ggplot2::geom_line(linewidth = 1.2) +
        ggplot2::geom_point(size = 2.5) +
        ggplot2::scale_y_log10(
            labels = scales::label_number(
                scale_cut = scales::cut_short_scale()
            )
        ) +
        ggplot2::scale_x_continuous(
            # n_resources is an integer count; force integer breaks +
            # minor gridlines so the eye reads the curve at 1/2/3/...
            breaks       = scales::breaks_width(1L),
            minor_breaks = NULL,
            expand       = ggplot2::expansion(mult = c(0.02, 0.02))
        ) +
        ggplot2::scale_colour_manual(
            values = cols,
            labels = variant_labels,
            name   = NULL
        ) +
        ggplot2::labs(
            x = "Number of resources (>= N)",
            y = "Items present in >= N resources (log scale)"
        ) +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(
            # 2.5x scale-up so labels read at composite-figure scale.
            legend.position = "bottom",
            legend.text     = ggplot2::element_text(size = 13),
            legend.key.size = grid::unit(7, "mm"),
            axis.text       = ggplot2::element_text(size = 13),
            axis.title      = ggplot2::element_text(size = 15)
        )
}
