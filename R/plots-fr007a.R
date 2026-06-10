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
#' @importFrom ggnewscale new_scale_fill
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

    # Per-facet major-class colour map (cycle the lead palette so the
    # same value across facets gets the same colour).
    facet_class_colours <- lapply(facet_order, function(f) {
        classes <- unique(as.character(
            data$category[data$facet == f &
                          data$bar_type == "major_class"]
        ))
        if (length(classes) == 0L) {
            return(stats::setNames(character(0L), character(0L)))
        }
        cols <- palette_n(min(length(classes), length(palette_lead())))
        cols <- rep(cols, length.out = length(classes))
        stats::setNames(cols, classes)
    })
    names(facet_class_colours) <- facet_order

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
            band_data           = band_data,
            band_label          = band,
            facet_order         = facet_order,
            facet_class_colours = facet_class_colours,
            ranking_facet       = ranking_facet,
            show_strip_x        = i == 1L,
            # Show x-axis tick labels on EVERY band so the per-band
            # x scale (each row has its own) is readable.
            show_axis_x         = TRUE,
            width_mm            = width_mm
        )
        # Height proportional to resource count in the band.
        heights[[i]] <- length(unique(band_data$resource_label))
    }

    # Suppress legends on all band plots — we'll extract each scale's
    # legend individually and arrange them in a grid below the plots.
    plots <- lapply(plots, function(p) {
        p + ggplot2::theme(legend.position = "none")
    })

    # Build one legend-only plot per scale, extract its guide-box,
    # then arrange the 7 legend grobs so the shared/unique legend
    # sits at the top of the legend area and the 6 per-facet
    # legends form a row that aligns with the plot's 6 facet
    # columns (each legend directly below its facet's plot column).
    su_legend <- fr007a_extract_legend(
        name   = "Occurrence across resources",
        values = c(unique = "#1B5E73", shared = "#A6D8E5"),
        ncol   = 2L
    )
    facet_legends <- lapply(facet_order, function(f) {
        cols <- facet_class_colours[[f]]
        if (length(cols) == 0L) {
            return(patchwork::plot_spacer())
        }
        grob <- fr007a_extract_legend(
            name   = paste0(f, " types"),
            values = cols,
            ncol   = if (length(cols) <= 4L) 1L else 2L
        )
        patchwork::wrap_elements(full = grob)
    })

    facet_legend_row <- patchwork::wrap_plots(
        facet_legends, nrow = 1L
    )

    # Heights: bands + shared/unique row + per-facet legends row.
    # Per-facet legends are taller because they have more entries.
    su_h    <- max(1L, as.integer(min(heights) * 0.25))
    facet_h <- max(2L, as.integer(min(heights) * 0.5))

    patchwork::wrap_plots(
        c(plots, list(
            patchwork::wrap_elements(full = su_legend),
            facet_legend_row
        )),
        ncol = 1L
    ) +
        patchwork::plot_layout(heights = c(heights, su_h, facet_h))
}


#' Build a legend-only ggplot for one fill scale and return its grob
#'
#' Renders a tiny dummy plot with the requested scale, then extracts
#' the bottom guide-box via \code{cowplot::get_plot_component}. The
#' guide is configured with the title above the entries (per user's
#' wish) and entries laid out in 1-2 columns depending on how many
#' there are.
#'
#' @param name Character: legend title.
#' @param values Named character vector: category → hex colour.
#' @param ncol Integer: number of columns the legend entries lay out
#'     into.
#'
#' @return A grob.
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_manual
#' @importFrom ggplot2 guide_legend theme_void theme element_text
#' @importFrom ggplot2 element_blank margin guides
#' @importFrom cowplot get_plot_component
#' @keywords internal
#' @noRd
fr007a_extract_legend <- function(name, values, ncol = 2L) {

    df <- data.frame(
        x = seq_along(values),
        y = 1L,
        c = factor(names(values), levels = names(values))
    )

    p <- ggplot2::ggplot(df, ggplot2::aes(.data$x, .data$y, fill = .data$c)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_manual(
            name   = name,
            values = values,
            limits = names(values),
            drop   = FALSE
        ) +
        ggplot2::guides(
            fill = ggplot2::guide_legend(
                ncol           = ncol,
                title.position = "top",
                title.hjust    = 0
            )
        ) +
        ggplot2::theme_void() +
        ggplot2::theme(
            legend.position    = "bottom",
            legend.title       = ggplot2::element_text(
                size = 5, face = "bold", hjust = 0,
                margin = ggplot2::margin(0, 0, 1, 0)
            ),
            legend.text        = ggplot2::element_text(size = 4.5),
            legend.key.size    = grid::unit(2.5, "mm"),
            legend.spacing.y   = grid::unit(0.5, "mm"),
            legend.box.spacing = grid::unit(0, "mm"),
            plot.margin        = ggplot2::margin(0, 1, 0, 1, "mm")
        )

    if (utils::packageVersion("cowplot") >= "1.1.3") {
        cowplot::get_plot_component(p, "guide-box-bottom", return_all = FALSE)
    } else {
        cowplot::get_legend(p)
    }
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
                             facet_order,
                             facet_class_colours,
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

    # Shared/unique data: factor order with unique first so
    # position_stack(reverse = FALSE) puts unique at x = 0 (closer
    # to the y axis) and shared stacks outward — consistent across
    # every panel.
    data_su <- band_data[
        as.character(band_data$bar_type) == "shared_unique", , drop = FALSE
    ]
    data_su$category <- factor(data_su$category,
                               levels = c("unique", "shared"))

    p <- ggplot2::ggplot() +
        # ---- Shared/unique layer + scale -------------------------------
        ggplot2::geom_col(
            data        = data_su,
            mapping     = ggplot2::aes(
                y     = .data$y_pos,
                x     = .data$n,
                fill  = .data$category,
                group = interaction(.data$resource_label, .data$bar_type)
            ),
            position    = ggplot2::position_stack(reverse = FALSE),
            width       = 0.4,
            orientation = "y"
        ) +
        ggplot2::scale_fill_manual(
            name   = "Occurrence across resources",
            values = c(unique = "#1B5E73", shared = "#A6D8E5"),
            breaks = c("unique", "shared"),
            limits = c("unique", "shared"),
            drop   = FALSE
        )

    # ---- Per-facet major-class layers + scales ----------------------
    # Always add ALL facets' scales (even when this band has no
    # data for one) so the ggnewscale aesthetic-rename chain has the
    # SAME length across all 3 band plots → patchwork's guides =
    # "collect" sees identical scale identities and dedups across
    # bands.
    for (f in facet_order) {
        class_colours <- facet_class_colours[[f]]
        if (length(class_colours) == 0L) next

        data_mc <- band_data[
            as.character(band_data$facet) == f &
                as.character(band_data$bar_type) == "major_class", ,
            drop = FALSE
        ]
        # When this band has no data for the facet, use an empty
        # tibble with the correct columns so the geom + scale still
        # add to the plot.
        if (nrow(data_mc) == 0L) {
            data_mc <- band_data[FALSE, , drop = FALSE]
            data_mc$category <- factor(character(0L),
                                       levels = names(class_colours))
        } else {
            data_mc$category <- factor(data_mc$category,
                                       levels = names(class_colours))
        }

        p <- p +
            ggnewscale::new_scale_fill() +
            ggplot2::geom_col(
                data        = data_mc,
                mapping     = ggplot2::aes(
                    y     = .data$y_pos,
                    x     = .data$n,
                    fill  = .data$category,
                    group = interaction(.data$resource_label,
                                        .data$bar_type)
                ),
                position    = ggplot2::position_stack(reverse = FALSE),
                width       = 0.4,
                orientation = "y"
            ) +
            ggplot2::scale_fill_manual(
                name   = paste0(f, " types"),
                values = class_colours,
                breaks = names(class_colours),
                limits = names(class_colours),
                drop   = FALSE
            )
    }

    p <- p +
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
