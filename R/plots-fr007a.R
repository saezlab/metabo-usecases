#' Render the FR-007a faceted resource overview
#'
#' Three independent per-band ggplots stacked vertically via
#' \code{patchwork}: large (>= 100K) on top, medium (10K-100K) in
#' the middle, small (< 10K) at the bottom. Each band-plot has its
#' own x-axis scale per facet column, its own y axis (only the
#' resources in that band, in descending order of total). The 6
#' facet-column strip labels appear only on the top plot. The
#' per-band band label sits in a left strip on each plot. A single
#' shared legend at the bottom of the composition (collected via
#' \code{patchwork::plot_layout(guides = "collect")}).
#'
#' Within each panel, every resource has TWO adjacent horizontal
#' bars at the same x-length: top = shared/unique split, bottom =
#' major-class breakdown.
#'
#' @param data Tibble from \code{\link{fr007a_overview}}.
#' @param facet_order Character vector: facet column order.
#' @param resource_order_facet Character: which facet's totals drive
#'     the within-band resource ordering. Default \code{"Entities"}.
#' @param width_mm Numeric: target physical width.
#'
#' @return A \code{patchwork} object (one combined ggplot).
#'
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual labs theme
#' @importFrom ggplot2 facet_grid vars element_text element_blank
#' @importFrom ggplot2 element_rect scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 expansion position_stack margin
#' @importFrom dplyr filter group_by summarise arrange mutate desc
#' @importFrom dplyr left_join
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom rlang .data
#' @export
plot_fr007a_overview <- function(data,
                                 facet_order = c(
                                     "Entities", "Associations",
                                     "Interactions", "Identifiers",
                                     "Structures", "Literature"
                                 ),
                                 resource_order_facet = "Entities",
                                 width_mm = 320L) {

    # NSE workaround
    facet <- magnitude_band <- NULL

    facet_order <- intersect(facet_order, unique(data$facet))
    data$facet <- factor(data$facet, levels = facet_order)

    ranking_facet <- if (resource_order_facet %in% facet_order) {
        resource_order_facet
    } else {
        facet_order[[1L]]
    }

    fill_map <- fr007a_fill_map(data)
    data$category <- factor(data$category, levels = fill_map$levels)

    bands <- levels(data$magnitude_band)
    bands <- bands[bands %in% unique(data$magnitude_band)]
    n_bands <- length(bands)

    plots <- vector("list", n_bands)
    heights <- integer(n_bands)
    for (i in seq_along(bands)) {
        band <- bands[[i]]
        band_data <- data[as.character(data$magnitude_band) == band, ,
                          drop = FALSE]
        plots[[i]] <- fr007a_band_plot(
            band_data     = band_data,
            band_label    = band,
            fill_map      = fill_map,
            ranking_facet = ranking_facet,
            show_strip_x  = i == 1L,
            show_axis_x   = i == n_bands,
            width_mm      = width_mm
        )
        # Height proportional to resource count in the band.
        heights[[i]] <- length(unique(band_data$resource_label))
    }

    patchwork::wrap_plots(plots, ncol = 1L) +
        patchwork::plot_layout(heights = heights, guides = "collect") &
        ggplot2::theme(legend.position = "bottom")
}


#' One band's plot inside the FR-007a composition
#'
#' \code{facet_grid(magnitude_band ~ facet)} so the left strip
#' carries the band label and each facet column gets its own x
#' scale via \code{scales = "free_x"}.
#'
#' @param band_data Tibble subset to one magnitude_band.
#' @param band_label Character: band display name.
#' @param fill_map List from \code{\link{fr007a_fill_map}}.
#' @param ranking_facet Character: facet whose totals drive the
#'     within-band ordering.
#' @param show_strip_x Logical: whether to render the top facet
#'     column strips.
#' @param show_axis_x Logical: whether to render the bottom x-axis
#'     tick labels.
#' @param width_mm Numeric.
#'
#' @return A ggplot.
#'
#' @keywords internal
#' @noRd
fr007a_band_plot <- function(band_data,
                             band_label,
                             fill_map,
                             ranking_facet,
                             show_strip_x,
                             show_axis_x,
                             width_mm) {

    # NSE
    facet <- resource_label <- n <- bar_type <- category <- NULL
    y_pos <- magnitude_band <- total <- NULL

    rank_in_band <- band_data %>%
        dplyr::filter(.data$facet == ranking_facet,
                      .data$bar_type == "shared_unique") %>%
        dplyr::group_by(.data$resource_label) %>%
        dplyr::summarise(
            total = sum(as.numeric(.data$n)),
            .groups = "drop"
        ) %>%
        dplyr::arrange(dplyr::desc(.data$total))

    ordered_lbl <- rank_in_band$resource_label
    extras <- setdiff(unique(band_data$resource_label), ordered_lbl)
    ordered_lbl <- c(ordered_lbl, extras)

    band_data <- band_data %>%
        dplyr::mutate(
            rank  = match(.data$resource_label, ordered_lbl),
            y_pos = -.data$rank +
                ifelse(.data$bar_type == "shared_unique", 0.22, -0.22)
        )

    breaks <- -seq_along(ordered_lbl)
    labels <- ordered_lbl

    p <- ggplot2::ggplot(
        band_data,
        ggplot2::aes(
            y    = .data$y_pos,
            x    = .data$n,
            fill = .data$category
        )
    ) +
        ggplot2::geom_col(
            position = ggplot2::position_stack(reverse = TRUE),
            width    = 0.4
        ) +
        ggplot2::scale_fill_manual(values = fill_map$hex,
                                   breaks = fill_map$levels,
                                   drop   = FALSE) +
        ggplot2::scale_x_continuous(
            expand = ggplot2::expansion(mult = c(0, 0.05)),
            labels = scales::label_number(
                scale_cut = scales::cut_short_scale()
            )
        ) +
        ggplot2::scale_y_continuous(
            breaks = breaks,
            labels = labels,
            expand = ggplot2::expansion(add = 0.5)
        ) +
        ggplot2::facet_grid(
            rows   = ggplot2::vars(.data$magnitude_band),
            cols   = ggplot2::vars(.data$facet),
            scales = "free_x",
            switch = NULL
        ) +
        ggplot2::labs(x = NULL, y = NULL, fill = NULL) +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(
            legend.position    = "bottom",
            legend.text        = ggplot2::element_text(size = 5),
            legend.key.size    = grid::unit(2.5, "mm"),
            legend.box.spacing = grid::unit(1, "mm"),
            axis.text.y        = ggplot2::element_text(size = 5),
            axis.text.x        = ggplot2::element_text(size = 4.5),
            strip.text.x       = if (show_strip_x) {
                ggplot2::element_text(face = "bold", size = 6)
            } else {
                ggplot2::element_blank()
            },
            strip.text.y       = ggplot2::element_text(
                face = "bold", size = 5, angle = 90
            ),
            strip.background.x = if (show_strip_x) {
                ggplot2::element_rect(fill = "#F0F0F0", colour = NA)
            } else {
                ggplot2::element_blank()
            },
            strip.background.y = ggplot2::element_rect(
                fill = "#F0F0F0", colour = NA
            ),
            panel.spacing.x    = grid::unit(2, "mm"),
            panel.spacing.y    = grid::unit(1, "mm")
        )

    # Hide bottom x-axis on non-bottom plots so axes don't repeat
    # between rows; keep ticks visible everywhere for readability.
    if (!show_axis_x) {
        p <- p + ggplot2::theme(
            axis.text.x  = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank()
        )
    }

    p
}


#' Build the FR-007a category → hex colour mapping
#'
#' \code{shared} and \code{unique} get a fixed light / dark pair;
#' major-class values cycle through the lead palette
#' deterministically.
#'
#' @keywords internal
#' @noRd
fr007a_fill_map <- function(data) {

    su_levels <- c("unique", "shared")
    su_hex    <- c(unique = "#1B5E73", shared = "#A6D8E5")

    class_levels <- unique(as.character(
        data[as.character(data$bar_type) == "major_class", "category",
             drop = TRUE]
    ))

    palette <- palette_n(min(length(class_levels), length(palette_lead())))
    palette <- rep(palette, length.out = length(class_levels))
    class_hex <- stats::setNames(palette, class_levels)

    list(
        hex    = c(su_hex, class_hex),
        levels = c(su_levels, class_levels)
    )
}
