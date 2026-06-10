#' Render the FR-007a faceted resource overview
#'
#' One ggplot with \code{facet_grid(rows = bar_type, cols = facet,
#' scales = "free_x")}: top row = shared/unique split, bottom row =
#' major-class breakdown; six facet columns (Entities, Associations,
#' Interactions, Identifiers, Structures, Literature). Each panel is
#' a horizontal stacked bar chart with one bar per resource. The y
#' axis is shared across all 12 panels so resources line up vertically.
#'
#' Resources are ordered by their total in the ranking facet
#' (\code{Entities} by default) with a leading \code{Total} row.
#' Resource display labels come from \code{resources.resource_short}
#' via the \code{resource_label} column the data layer adds.
#'
#' @param data Tibble from \code{\link{fr007a_overview}} — MUST
#'     include a \code{resource_label} column (added by
#'     \code{\link{resources_label_map}} inside the combiner).
#' @param facet_order Character vector: facet column order.
#'     Default Entities → Literature per the spec.
#' @param resource_order_facet Character: which facet's totals
#'     determine the resource ordering. Default \code{"Entities"}.
#' @param width_mm Numeric: target physical width.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual
#' @importFrom ggplot2 facet_grid vars labeller label_value labs theme
#' @importFrom ggplot2 element_text element_blank element_rect
#' @importFrom ggplot2 scale_x_continuous expansion
#' @importFrom ggplot2 position_stack
#' @importFrom dplyr group_by summarise arrange filter mutate desc
#' @importFrom dplyr ungroup
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
    facet <- resource <- bar_type <- category <- n <- total <- NULL
    resource_label <- NULL

    facet_order <- intersect(facet_order, unique(data$facet))
    data$facet <- factor(data$facet, levels = facet_order)
    data$bar_type <- factor(
        data$bar_type,
        levels = c("shared_unique", "major_class"),
        labels = c("shared / unique", "by major class")
    )

    # Resource ordering — descending by Total in the selected facet,
    # "Total" first.
    ranking_facet <- if (resource_order_facet %in% facet_order) {
        resource_order_facet
    } else {
        facet_order[[1L]]
    }
    rank_tbl <- data %>%
        dplyr::filter(.data$facet == ranking_facet,
                      .data$bar_type == "shared / unique") %>%
        dplyr::group_by(.data$resource_label) %>%
        dplyr::summarise(total = sum(.data$n), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(.data$total))

    ordered_labels <- c(
        "Total",
        setdiff(rank_tbl$resource_label[rank_tbl$resource_label != "Total"],
                "Total")
    )
    extras <- sort(setdiff(unique(data$resource_label), ordered_labels))
    ordered_labels <- c(ordered_labels, extras)

    # ggplot draws low y at the bottom; reversing puts Total at the
    # top of the panel.
    data <- data %>%
        dplyr::mutate(
            resource_label = factor(
                .data$resource_label, levels = rev(ordered_labels)
            )
        )

    fill_map <- fr007a_fill_map(data)
    data <- data %>%
        dplyr::mutate(
            category = factor(.data$category, levels = fill_map$levels)
        )

    ggplot2::ggplot(
        data,
        ggplot2::aes(
            y    = .data$resource_label,
            x    = .data$n,
            fill = .data$category
        )
    ) +
        ggplot2::geom_col(position = ggplot2::position_stack(reverse = TRUE),
                          width = 0.75) +
        ggplot2::scale_fill_manual(values = fill_map$hex,
                                   breaks = fill_map$levels) +
        ggplot2::scale_x_continuous(
            expand = ggplot2::expansion(mult = c(0, 0.05))
        ) +
        ggplot2::facet_grid(
            rows = ggplot2::vars(.data$bar_type),
            cols = ggplot2::vars(.data$facet),
            scales = "free_x"
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
            strip.text.y       = ggplot2::element_text(face = "bold", size = 5),
            strip.background   = ggplot2::element_rect(
                fill = "#F0F0F0", colour = NA
            ),
            panel.spacing.x    = grid::unit(2, "mm"),
            panel.spacing.y    = grid::unit(2, "mm")
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
        data[as.character(data$bar_type) == "by major class", "category",
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
