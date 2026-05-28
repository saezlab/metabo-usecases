#' Panel B — entities by resource
#'
#' Horizontal bar chart, one bar per resource, coloured via the
#' \code{resources} category registry. Uniform bar widths (FR-019).
#' Caller MUST register every resource in the registry first (a
#' missing entry fails fast per the spec Edge Case).
#'
#' @param data Tibble from \code{\link{entities_by_resource}}.
#' @param width_mm Numeric: target physical width for the panel.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_col coord_flip scale_fill_manual
#' @importFrom ggplot2 labs scale_y_continuous expansion
#' @importFrom rlang .data
#' @export
plot_entities_by_resource <- function(data, width_mm = 89L) {

    assert_category_known("resources", data$resource)
    fills <- category_colour("resources", data$resource)

    ggplot2::ggplot(
        data,
        ggplot2::aes(
            x    = stats::reorder(.data$resource, .data$n_entities),
            y    = .data$n_entities,
            fill = .data$resource
        )
    ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::scale_fill_manual(values = setNames(fills, data$resource)) +
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = c(0, 0.05))
        ) +
        ggplot2::coord_flip() +
        ggplot2::labs(x = NULL, y = "Entities") +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(legend.position = "none")
}


#' Panel C — interactions by resource
#'
#' Same shape as Panel B; reuses the resource colour mapping.
#'
#' @param data Tibble from \code{\link{interactions_by_resource}}.
#' @param width_mm Numeric.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_col coord_flip scale_fill_manual
#' @importFrom ggplot2 labs scale_y_continuous expansion theme
#' @importFrom rlang .data
#' @export
plot_interactions_by_resource <- function(data, width_mm = 89L) {

    assert_category_known("resources", data$resource)
    fills <- category_colour("resources", data$resource)

    ggplot2::ggplot(
        data,
        ggplot2::aes(
            x    = stats::reorder(.data$resource, .data$n_relations),
            y    = .data$n_relations,
            fill = .data$resource
        )
    ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::scale_fill_manual(values = setNames(fills, data$resource)) +
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = c(0, 0.05))
        ) +
        ggplot2::coord_flip() +
        ggplot2::labs(x = NULL, y = "Interactions") +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(legend.position = "none")
}


#' Panel D — interactions by type
#'
#' Vertical bar chart coloured by interaction-type registry. Uniform
#' bar widths matching panels B and C.
#'
#' @param data Tibble from \code{\link{interactions_by_type}}.
#' @param width_mm Numeric.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual labs
#' @importFrom ggplot2 scale_y_continuous expansion theme element_text
#' @importFrom rlang .data
#' @export
plot_interactions_by_type <- function(data, width_mm = 89L) {

    assert_category_known("interaction_types", data$interaction_type)
    fills <- category_colour(
        "interaction_types", data$interaction_type
    )

    ggplot2::ggplot(
        data,
        ggplot2::aes(
            x    = stats::reorder(.data$interaction_type, -.data$n),
            y    = .data$n,
            fill = .data$interaction_type
        )
    ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::scale_fill_manual(
            values = setNames(fills, data$interaction_type)
        ) +
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = c(0, 0.05))
        ) +
        ggplot2::labs(x = NULL, y = "Interactions") +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(
            legend.position = "none",
            axis.text.x     = ggplot2::element_text(angle = 30, hjust = 1)
        )
}


#' Panel E — annotation classes by resource
#'
#' @param data Tibble from \code{\link{annotation_classes_by_resource}}.
#' @param width_mm Numeric.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_col coord_flip scale_fill_manual
#' @importFrom ggplot2 labs scale_y_continuous expansion theme
#' @importFrom rlang .data
#' @export
plot_annotation_classes_by_resource <- function(data, width_mm = 89L) {

    assert_category_known("resources", data$resource)
    fills <- category_colour("resources", data$resource)

    ggplot2::ggplot(
        data,
        ggplot2::aes(
            x    = stats::reorder(.data$resource, .data$n_classes),
            y    = .data$n_classes,
            fill = .data$resource
        )
    ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::scale_fill_manual(
            values = setNames(fills, data$resource)
        ) +
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = c(0, 0.05))
        ) +
        ggplot2::coord_flip() +
        ggplot2::labs(x = NULL, y = "Annotation classes") +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(legend.position = "none")
}


#' Panel F — ontology terms by ontology
#'
#' @param data Tibble from \code{\link{ontology_terms_by_ontology}}.
#' @param width_mm Numeric.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual labs
#' @importFrom ggplot2 scale_y_continuous expansion theme element_text
#' @importFrom rlang .data
#' @export
plot_ontology_terms_by_ontology <- function(data, width_mm = 89L) {

    fills <- palette_n(nrow(data), unknown = FALSE)

    ggplot2::ggplot(
        data,
        ggplot2::aes(
            x    = stats::reorder(.data$ontology, -.data$n_terms),
            y    = .data$n_terms,
            fill = .data$ontology
        )
    ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::scale_fill_manual(values = setNames(fills, data$ontology)) +
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = c(0, 0.05))
        ) +
        ggplot2::labs(x = NULL, y = "Ontology terms") +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(
            legend.position = "none",
            axis.text.x     = ggplot2::element_text(angle = 30, hjust = 1)
        )
}
