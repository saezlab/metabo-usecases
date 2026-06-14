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

    # Per-facet class-only maps (for the horizontal-facets +
    # per-facet-legend layout used by plot_fr007a_total).
    facets <- unique(as.character(data$facet))
    facet_class_hex <- stats::setNames(
        lapply(facets, function(f) {
            cls <- unique(as.character(
                data[as.character(data$facet) == f &
                         as.character(data$bar_type) == "major_class",
                     "category", drop = TRUE]
            ))
            cls <- cls[!is.na(cls) & nzchar(cls)]
            pal <- palette_n(min(length(cls), length(palette_lead())))
            pal <- rep(pal, length.out = length(cls))
            c(su_hex, stats::setNames(pal, cls))
        }),
        facets
    )
    facet_class_levels <- stats::setNames(
        lapply(facets, function(f) {
            cls <- unique(as.character(
                data[as.character(data$facet) == f &
                         as.character(data$bar_type) == "major_class",
                     "category", drop = TRUE]
            ))
            cls <- cls[!is.na(cls) & nzchar(cls)]
            c(su_levels, cls)
        }),
        facets
    )

    list(
        hex                = c(su_hex, class_hex),
        levels             = c(su_levels, class_levels),
        facet_hex          = facet_class_hex,
        facet_levels       = facet_class_levels,
        shared_unique      = su_levels
    )
}


#' Render the FR-007a tiny Total-only variant
#'
#' Compact alternative to \code{\link{plot_fr007a_overview}}: instead
#' of the band-resource matrix, show only the \code{Total} row across
#' all six facets, arranged HORIZONTALLY (one facet per column).
#' Each facet shows two adjacent horizontal bars stacked along the
#' y-axis — the \code{shared / unique} bar at the TOP and the
#' \code{major_class} stacked breakdown at the BOTTOM. Inside the
#' shared/unique bar, \code{unique} (dark) sits at the BASE of the
#' bar (next to the y-axis) and \code{shared} (light) stacks on top.
#'
#' Each facet carries its own legend below it (Session 2026-06-14:
#' "per-facet legends with scale titles") so a reader can decode the
#' colour scale facet-by-facet without scanning a single global key.
#' The six sub-plots are assembled with \pkg{patchwork} in one row.
#'
#' This variant is what lands in Figure 2 (post-2026-06-14 six-figure
#' renumbering); the full faceted-resource overview from
#' \code{\link{plot_fr007a_overview}} stays as a supplementary
#' artifact (\code{fr007a-overview-supplementary.pdf}).
#'
#' @param data Tibble from \code{\link{fr007a_overview}}.
#' @param facet_order Character vector: facet column order.
#' @param width_mm Numeric: target physical width of the combined
#'     six-facet row (each sub-plot gets \code{width_mm / 6}).
#'
#' @return A patchwork composite.
#'
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual labs theme
#' @importFrom ggplot2 element_text element_blank
#' @importFrom ggplot2 scale_x_continuous scale_y_discrete expansion
#' @importFrom ggplot2 guides guide_legend position_stack
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom rlang .data
#' @export
plot_fr007a_total <- function(data,
                              facet_order = c(
                                  "Entities", "Associations",
                                  "Interactions", "Identifiers",
                                  "Structures", "Literature"
                              ),
                              width_mm = 180L) {

    fill_map <- fr007a_fill_map(data)

    total <- data[as.character(data$resource) == "Total", , drop = FALSE]
    facet_order <- intersect(facet_order, unique(total$facet))

    # bar_type levels: major_class at BOTTOM, shared_unique on TOP.
    # The y-axis is discrete and the FIRST level appears at the
    # bottom of the panel — so listing major_class first puts the
    # categories bar below the unique/shared bar (Session 2026-06-14
    # bar-order swap).
    total$bar_type <- factor(
        total$bar_type,
        levels = c("major_class", "shared_unique")
    )

    sub_plots <- lapply(facet_order, function(f) {
        d <- total[as.character(total$facet) == f, , drop = FALSE]
        levels_f <- fill_map$facet_levels[[f]]
        hex_f    <- fill_map$facet_hex[[f]]
        d$category <- factor(d$category, levels = levels_f)

        ggplot2::ggplot(
            d,
            ggplot2::aes(
                x    = .data$n,
                y    = .data$bar_type,
                fill = .data$category
            )
        ) +
            # reverse=TRUE puts the FIRST factor level (e.g. "unique")
            # at the BASE of the stacked bar (x=0, next to the
            # y-axis) instead of the default last-first stacking.
            ggplot2::geom_col(
                position = ggplot2::position_stack(reverse = TRUE)
            ) +
            ggplot2::scale_fill_manual(
                values = hex_f,
                name   = f,
                breaks = levels_f
            ) +
            ggplot2::scale_x_continuous(
                labels = scales::label_number(
                    scale_cut = scales::cut_short_scale()
                ),
                expand = ggplot2::expansion(mult = c(0, 0.04))
            ) +
            ggplot2::scale_y_discrete(labels = NULL) +
            ggplot2::labs(x = NULL, y = NULL, title = f) +
            theme_bw_metabo(width_mm = width_mm / length(facet_order)) +
            ggplot2::guides(
                fill = ggplot2::guide_legend(
                    title.position = "top",
                    title.hjust    = 0.5,
                    ncol           = 1L
                )
            ) +
            ggplot2::theme(
                plot.title         = ggplot2::element_text(
                    size = 12, face = "bold", hjust = 0.5
                ),
                legend.position    = "bottom",
                legend.title       = ggplot2::element_text(
                    size = 10, face = "bold"
                ),
                legend.text        = ggplot2::element_text(size = 9),
                legend.key.size    = grid::unit(3, "mm"),
                legend.box.margin  = ggplot2::margin(t = 2, r = 0,
                                                    b = 0, l = 0),
                strip.background   = ggplot2::element_blank(),
                axis.ticks.y       = ggplot2::element_blank(),
                axis.text.y        = ggplot2::element_blank(),
                axis.text.x        = ggplot2::element_text(size = 9),
                panel.spacing.y    = grid::unit(2, "mm")
            )
    })

    patchwork::wrap_plots(sub_plots, nrow = 1L) +
        patchwork::plot_layout(guides = "keep")
}
