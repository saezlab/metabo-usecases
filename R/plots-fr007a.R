#' Render the FR-007a faceted resource overview
#'
#' One ggplot with \code{facet_grid(rows = magnitude_band, cols =
#' facet, scales = "free", space = "free_y")}: 6 facet columns,
#' 3 rows (large / medium / small) so the ChEMBL-sized tail does
#' not squash everyone onto an invisible scale. Within each panel
#' every resource has TWO adjacent horizontal bars at the same
#' x-length: top = shared/unique, bottom = major-class breakdown.
#'
#' Each (magnitude_band, facet) panel auto-zooms its own x axis via
#' \code{scales = "free"} — the small band shows counts at 10x
#' lower magnitude than the large band so its bars are not crushed.
#' Each band's panel height is proportional to the resource count
#' in that band via \code{space = "free_y"}, so bars look the same
#' physical width across bands.
#'
#' Resource ordering: rank globally as (band ASC, total DESC). The
#' Total row gets rank 1 (top of the figure). All resources have a
#' unique \code{global_rank}, so the y-axis breaks vector contains
#' no duplicates — fixes v3/v4 where per-band ranks collided in the
#' shared scale.
#'
#' @param data Tibble from \code{\link{fr007a_overview}}.
#' @param facet_order Character vector: facet column order.
#' @param resource_order_facet Character: which facet's totals drive
#'     within-band ordering. Default \code{"Entities"}.
#' @param width_mm Numeric: target physical width.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual
#' @importFrom ggplot2 facet_grid vars labs theme
#' @importFrom ggplot2 element_text element_blank element_rect
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous expansion
#' @importFrom ggplot2 position_stack
#' @importFrom dplyr group_by summarise arrange filter mutate desc ungroup
#' @importFrom dplyr left_join bind_rows select distinct row_number
#' @importFrom tibble tibble
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
    facet <- resource <- bar_type <- category <- n <- NULL
    resource_label <- magnitude_band <- global_rank <- y_pos <- NULL
    total <- NULL

    facet_order <- intersect(facet_order, unique(data$facet))
    data$facet <- factor(data$facet, levels = facet_order)

    ranking_facet <- if (resource_order_facet %in% facet_order) {
        resource_order_facet
    } else {
        facet_order[[1L]]
    }

    # Globally-unique rank — large-band resources first (descending
    # total), then medium, then small. The magnitude_band factor
    # levels (set in the data layer) drive the band ordering.
    rank_tbl <- data %>%
        dplyr::filter(.data$facet == ranking_facet,
                      .data$bar_type == "shared_unique") %>%
        dplyr::group_by(.data$resource_label, .data$magnitude_band) %>%
        dplyr::summarise(
            total = sum(as.numeric(.data$n)),
            .groups = "drop"
        ) %>%
        dplyr::arrange(.data$magnitude_band, dplyr::desc(.data$total)) %>%
        dplyr::mutate(global_rank = dplyr::row_number()) %>%
        dplyr::select(.data$resource_label, .data$magnitude_band,
                      .data$global_rank)

    # Resources that don't appear in the ranking facet (e.g. ontology-
    # only contributors for Associations). Drop them in at the tail
    # of their band so nothing is silently lost.
    extras_lbl <- setdiff(data$resource_label, rank_tbl$resource_label)
    if (length(extras_lbl) > 0L) {
        extras_band <- vapply(extras_lbl, function(lbl) {
            as.character(unique(
                data$magnitude_band[data$resource_label == lbl]
            ))[[1L]]
        }, character(1L))
        rank_tbl <- dplyr::bind_rows(
            rank_tbl,
            tibble::tibble(
                resource_label = extras_lbl,
                magnitude_band = extras_band,
                global_rank    = max(rank_tbl$global_rank) +
                                 seq_along(extras_lbl)
            )
        )
    }

    data <- data %>%
        dplyr::left_join(
            dplyr::select(rank_tbl, .data$resource_label,
                          .data$global_rank),
            by = "resource_label"
        ) %>%
        dplyr::mutate(
            # Continuous y so we can offset the two bar types around
            # each resource's integer slot. Negative so that rank 1
            # (Total) sits at the TOP of the y axis (largest y value).
            y_pos = -.data$global_rank +
                ifelse(.data$bar_type == "shared_unique", 0.22, -0.22)
        )

    fill_map <- fr007a_fill_map(data)
    data <- data %>%
        dplyr::mutate(
            category = factor(.data$category, levels = fill_map$levels)
        )

    # Scale's breaks/labels are GLOBAL (one entry per resource at its
    # unique global_rank). Each panel auto-clips to its data's y
    # range via scales = "free", so only the labels for that band's
    # resources appear in each band-row.
    breaks <- -rank_tbl$global_rank
    labels <- rank_tbl$resource_label

    ggplot2::ggplot(
        data,
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
                                   breaks = fill_map$levels) +
        ggplot2::scale_x_continuous(
            expand = ggplot2::expansion(mult = c(0, 0.05)),
            labels = scales::label_number(scale_cut = scales::cut_short_scale())
        ) +
        ggplot2::scale_y_continuous(
            breaks = breaks,
            labels = labels,
            expand = ggplot2::expansion(add = 0.5)
        ) +
        ggplot2::facet_grid(
            rows   = ggplot2::vars(.data$magnitude_band),
            cols   = ggplot2::vars(.data$facet),
            scales = "free",
            space  = "free_y"
        ) +
        ggplot2::labs(x = "Count", y = NULL, fill = NULL) +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(
            legend.position    = "bottom",
            legend.text        = ggplot2::element_text(size = 5),
            legend.key.size    = grid::unit(2.5, "mm"),
            legend.box.spacing = grid::unit(1, "mm"),
            axis.text.y        = ggplot2::element_text(size = 5),
            axis.text.x        = ggplot2::element_text(size = 4.5),
            strip.text.x       = ggplot2::element_text(face = "bold", size = 6),
            strip.text.y       = ggplot2::element_text(face = "bold", size = 5,
                                                       angle = 0),
            strip.background   = ggplot2::element_rect(
                fill = "#F0F0F0", colour = NA
            ),
            panel.spacing.x    = grid::unit(2, "mm"),
            panel.spacing.y    = grid::unit(1, "mm")
        )
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
