#' Render the FR-007a faceted resource overview
#'
#' One ggplot with \code{facet_grid(rows = magnitude_band, cols =
#' facet, scales = "free", space = "free_y")}: 6 facet columns,
#' 2–3 rows split by resource magnitude so a single mega-resource
#' (ChEMBL with 2.5 M entities) does not squash everyone else onto
#' an invisible scale. Within each panel, every resource has TWO
#' adjacent horizontal bars at the same length: top = shared/unique
#' split, bottom = major-class breakdown.
#'
#' @param data Tibble from \code{\link{fr007a_overview}}: must
#'     include \code{resource_label} and \code{magnitude_band}
#'     columns (the data layer adds them).
#' @param facet_order Character vector: facet column order.
#' @param resource_order_facet Character: which facet's totals drive
#'     the ranking within each band. Default \code{"Entities"}.
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
#' @importFrom dplyr left_join
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
    resource_label <- magnitude_band <- band_rank <- y_pos <- NULL
    total <- NULL

    facet_order <- intersect(facet_order, unique(data$facet))
    data$facet <- factor(data$facet, levels = facet_order)

    # Resource ranking within each magnitude band: descending by
    # ranking-facet shared/unique total.
    ranking_facet <- if (resource_order_facet %in% facet_order) {
        resource_order_facet
    } else {
        facet_order[[1L]]
    }

    rank_tbl <- data %>%
        dplyr::filter(.data$facet == ranking_facet,
                      .data$bar_type == "shared_unique") %>%
        dplyr::group_by(.data$resource_label, .data$magnitude_band) %>%
        dplyr::summarise(total = sum(.data$n), .groups = "drop") %>%
        dplyr::group_by(.data$magnitude_band) %>%
        dplyr::arrange(dplyr::desc(.data$total), .by_group = TRUE) %>%
        dplyr::mutate(band_rank = dplyr::row_number()) %>%
        dplyr::ungroup() %>%
        dplyr::select(.data$resource_label, .data$band_rank)

    # Fallback for resources that don't appear in the ranking facet
    # (e.g. ontology-only resources for Associations). Append them at
    # the end of their band with arbitrary order.
    extras <- setdiff(data$resource_label, rank_tbl$resource_label)
    if (length(extras) > 0L) {
        rank_tbl <- dplyr::bind_rows(
            rank_tbl,
            tibble::tibble(
                resource_label = extras,
                band_rank      = max(rank_tbl$band_rank,
                                     na.rm = TRUE) + seq_along(extras)
            )
        )
    }

    data <- data %>%
        dplyr::left_join(rank_tbl, by = "resource_label") %>%
        dplyr::mutate(
            # Within-panel two-bar layout: top = shared_unique, bottom
            # = major_class. Continuous y so we can offset bars by
            # +/-0.22 around each resource's integer slot. Bands flip
            # the sign so the highest-rank resource sits at the TOP
            # of the panel (largest y).
            y_pos = -.data$band_rank
                + ifelse(.data$bar_type == "shared_unique", 0.22, -0.22)
        )

    fill_map <- fr007a_fill_map(data)
    data <- data %>%
        dplyr::mutate(
            category = factor(.data$category, levels = fill_map$levels)
        )

    # Per-band y-axis labels — derived from band_rank → resource_label.
    # scale_y_continuous can't carry different labels per panel row,
    # so we use a custom labeller via scale_y_continuous(breaks=,
    # labels=) computed from the FULL data and rely on
    # `space="free_y"` to drop unused rows in each panel.
    label_tbl <- rank_tbl %>%
        dplyr::distinct(.data$resource_label, .data$band_rank) %>%
        dplyr::arrange(.data$band_rank)

    breaks <- -label_tbl$band_rank
    labels <- label_tbl$resource_label

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
            expand = ggplot2::expansion(mult = c(0, 0.05))
        ) +
        ggplot2::scale_y_continuous(
            breaks = breaks,
            labels = labels,
            expand = ggplot2::expansion(add = 0.6)
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
#' Produces a named character vector keyed by every category value
#' present in the data. \code{shared} and \code{unique} get a fixed
#' light / dark pair; major-class values cycle through the lead
#' palette deterministically (same value across facets gets the same
#' colour).
#'
#' @param data Tibble from \code{\link{fr007a_overview}}.
#'
#' @return A list with \code{hex} (named character vector) and
#'     \code{levels} (categories in factor order).
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
