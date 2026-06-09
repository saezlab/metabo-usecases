#' Render the FR-007a faceted resource overview
#'
#' One ggplot with \code{facet_wrap(~ facet)} per the spec
#' (Session 2026-06-02). Each (facet, resource) pair gets two
#' horizontal stacked bars:
#'   1. top bar — \code{shared} / \code{unique} (light/dark of the
#'      facet's shade);
#'   2. bottom bar — major-class breakdown coloured from the lead
#'      palette via the category registry.
#'
#' Resources are ordered by their size in the first facet
#' (\code{Entities} by default) with a leading \code{Total} row.
#' Per FR-009b/c, when the major-class palette would exceed the
#' 10-12 colour cap, the long tail is already collapsed at the data
#' layer (e.g. Identifiers Misc).
#'
#' @param data Tibble from \code{\link{fr007a_overview}}.
#' @param facet_order Character vector: explicit facet order.
#'     Defaults to spec order (Entities, Associations, Interactions,
#'     Identifiers, Structures, Literature).
#' @param resource_order_facet Character: which facet's totals
#'     determine the resource ordering. Default \code{"Entities"}.
#' @param width_mm Numeric: target physical width.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_col coord_flip scale_fill_manual
#' @importFrom ggplot2 facet_wrap labs theme element_text element_blank
#' @importFrom ggplot2 scale_x_discrete scale_y_continuous expansion
#' @importFrom ggplot2 position_stack
#' @importFrom dplyr group_by summarise arrange filter mutate
#' @importFrom dplyr ungroup desc bind_rows
#' @importFrom rlang .data
#' @importFrom tibble tibble
#' @export
plot_fr007a_overview <- function(data,
                                 facet_order = c(
                                     "Entities", "Associations",
                                     "Interactions", "Identifiers",
                                     "Structures", "Literature"
                                 ),
                                 resource_order_facet = "Entities",
                                 width_mm = 180L) {

    # NSE workaround
    facet <- resource <- bar_type <- category <- n <- y_pos <- NULL
    fill_hex <- y_label <- NULL

    facet_order <- intersect(facet_order, unique(data$facet))
    data$facet <- factor(data$facet, levels = facet_order)

    # Resource ordering — descending by Total in the selected facet,
    # "Total" first.
    ranking_facet <- if (resource_order_facet %in% facet_order) {
        resource_order_facet
    } else {
        facet_order[[1L]]
    }
    rank_tbl <- data %>%
        dplyr::filter(.data$facet == ranking_facet,
                      .data$bar_type == "shared_unique") %>%
        dplyr::group_by(.data$resource) %>%
        dplyr::summarise(total = sum(.data$n), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(.data$total))

    ordered_resources <- c(
        "Total",
        setdiff(rank_tbl$resource[rank_tbl$resource != "Total"],
                "Total")
    )
    # Resources that appear in some other facet but not in the
    # ranking facet — keep them at the bottom in alphabetical order
    # so nothing is silently dropped.
    extras <- sort(setdiff(unique(data$resource), ordered_resources))
    ordered_resources <- c(ordered_resources, extras)

    data <- data %>%
        dplyr::mutate(
            resource = factor(.data$resource, levels = rev(ordered_resources))
        )

    # Build a single hex fill column.
    # shared/unique colours come from a fixed pair (dark blue +
    # light blue). Major-class colours come from
    # category_colour("fr007a_<facet>", value); register on first
    # use.
    fill_map <- fr007a_fill_map(data)

    data <- data %>%
        dplyr::mutate(
            fill_key = paste(.data$bar_type, .data$category, sep = "/")
        )

    # Stack ordering: shared on top of unique (shared = lighter).
    data <- data %>%
        dplyr::mutate(
            category = factor(.data$category, levels = fill_map$levels)
        )

    # Bar offsets: shared/unique above major_class.
    data <- data %>%
        dplyr::mutate(
            bar_offset = ifelse(.data$bar_type == "shared_unique", 0.18, -0.18)
        )

    ggplot2::ggplot(
        data,
        ggplot2::aes(
            y    = as.numeric(.data$resource) + .data$bar_offset,
            x    = .data$n,
            fill = .data$category
        )
    ) +
        ggplot2::geom_col(
            position = ggplot2::position_stack(),
            width    = 0.32
        ) +
        ggplot2::scale_fill_manual(values = fill_map$hex,
                                   breaks = fill_map$levels) +
        ggplot2::scale_y_continuous(
            breaks = seq_along(levels(data$resource)),
            labels = levels(data$resource),
            expand = ggplot2::expansion(add = 0.6)
        ) +
        ggplot2::scale_x_continuous(
            expand = ggplot2::expansion(mult = c(0, 0.05))
        ) +
        ggplot2::facet_wrap(~ facet, scales = "free_x", nrow = 2L) +
        ggplot2::labs(x = "Count", y = NULL, fill = NULL) +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(
            legend.position    = "bottom",
            legend.text        = ggplot2::element_text(size = 5),
            legend.key.size    = grid::unit(2.5, "mm"),
            panel.spacing.y    = grid::unit(2, "mm"),
            strip.text         = ggplot2::element_text(face = "bold")
        )
}


#' Build the FR-007a category → hex colour mapping
#'
#' Produces a named character vector keyed by every category value
#' present in the data. \code{shared} / \code{unique} get a fixed
#' light / dark pair; major-class values cycle through the lead
#' palette, with already-registered values in the category-colour
#' registry honoured first (so resource colours match across
#' figures).
#'
#' @param data Tibble from \code{\link{fr007a_overview}} after the
#'     facet factor has been set.
#'
#' @return A list with \code{hex} (named character vector) and
#'     \code{levels} (categories in factor order).
#'
#' @importFrom dplyr distinct
#' @keywords internal
#' @noRd
fr007a_fill_map <- function(data) {

    # Shared/unique — fixed dark = unique, light = shared
    su <- c(
        unique = "#1B5E73",
        shared = "#A6D8E5"
    )

    # Major-class values across all facets, ordered by first
    # appearance.
    class_values <- data[data$bar_type == "major_class", "category",
                         drop = TRUE]
    class_values <- unique(as.character(class_values))

    # Cycle the lead palette deterministically; same value across
    # facets gets the same colour.
    palette <- palette_n(min(length(class_values),
                             length(palette_lead())))
    palette <- rep(palette, length.out = length(class_values))
    class_hex <- stats::setNames(palette, class_values)

    list(
        hex    = c(su, class_hex),
        levels = c(names(su), class_values)
    )
}
