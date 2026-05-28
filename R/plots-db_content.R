#' Panel B -- entities by resource
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
#' @importFrom ggplot2 ggplot aes geom_col coord_flip scale_fill_identity
#' @importFrom ggplot2 labs scale_y_continuous expansion theme
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @export
plot_entities_by_resource <- function(data, width_mm = 89L) {

    # NSE workaround
    resource <- n_entities <- fill_hex <- NULL

    data <- dplyr::mutate(
        data,
        fill_hex = category_colour("resources", resource)
    )

    ggplot2::ggplot(
        data,
        ggplot2::aes(
            x    = stats::reorder(.data$resource, .data$n_entities),
            y    = .data$n_entities,
            fill = .data$fill_hex
        )
    ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::scale_fill_identity() +
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = c(0, 0.05))
        ) +
        ggplot2::coord_flip() +
        ggplot2::labs(x = NULL, y = "Entities") +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(legend.position = "none")
}


#' Panel C -- interactions by resource
#'
#' Same shape as Panel B; reuses the resource colour mapping.
#'
#' @param data Tibble from \code{\link{interactions_by_resource}}.
#' @param width_mm Numeric.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_col coord_flip scale_fill_identity
#' @importFrom ggplot2 labs scale_y_continuous expansion theme
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @export
plot_interactions_by_resource <- function(data, width_mm = 89L) {

    resource <- n_relations <- fill_hex <- NULL

    data <- dplyr::mutate(
        data,
        fill_hex = category_colour("resources", resource)
    )

    ggplot2::ggplot(
        data,
        ggplot2::aes(
            x    = stats::reorder(.data$resource, .data$n_relations),
            y    = .data$n_relations,
            fill = .data$fill_hex
        )
    ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::scale_fill_identity() +
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = c(0, 0.05))
        ) +
        ggplot2::coord_flip() +
        ggplot2::labs(x = NULL, y = "Interactions") +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(legend.position = "none")
}


#' Panel D -- interactions by type
#'
#' Vertical bar chart coloured by interaction-type registry. Uniform
#' bar widths matching panels B and C.
#'
#' @param data Tibble from \code{\link{interactions_by_type}}.
#' @param width_mm Numeric.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_identity labs
#' @importFrom ggplot2 scale_y_continuous expansion theme element_text
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @export
plot_interactions_by_type <- function(data, width_mm = 89L) {

    interaction_type <- n <- fill_hex <- NULL

    data <- dplyr::mutate(
        data,
        fill_hex = category_colour("interaction_types", interaction_type)
    )

    ggplot2::ggplot(
        data,
        ggplot2::aes(
            x    = stats::reorder(.data$interaction_type, -.data$n),
            y    = .data$n,
            fill = .data$fill_hex
        )
    ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::scale_fill_identity() +
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


#' Panel E -- annotation classes by resource
#'
#' Empty-data case emits a single-bar placeholder so the panel layout
#' remains stable when the snapshot has no annotation relations yet
#' (spec Edge Case: zero-row queries get an explicit placeholder).
#'
#' @param data Tibble from \code{\link{annotation_classes_by_resource}}.
#' @param width_mm Numeric.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_col coord_flip scale_fill_identity
#' @importFrom ggplot2 labs scale_y_continuous expansion theme
#' @importFrom dplyr mutate
#' @importFrom tibble tibble
#' @importFrom rlang .data
#' @export
plot_annotation_classes_by_resource <- function(data, width_mm = 89L) {

    resource <- n_classes <- fill_hex <- NULL

    if (nrow(data) == 0L) {
        data <- tibble::tibble(
            resource   = "(none)",
            n_classes  = 0L,
            fill_hex   = "#BEBEBE"
        )
    } else {
        data <- dplyr::mutate(
            data,
            fill_hex = category_colour("resources", resource)
        )
    }

    ggplot2::ggplot(
        data,
        ggplot2::aes(
            x    = stats::reorder(.data$resource, .data$n_classes),
            y    = .data$n_classes,
            fill = .data$fill_hex
        )
    ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::scale_fill_identity() +
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = c(0, 0.05))
        ) +
        ggplot2::coord_flip() +
        ggplot2::labs(x = NULL, y = "Annotation classes") +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(legend.position = "none")
}


#' Panel F -- ontology terms by ontology
#'
#' Ontologies are not currently tracked in the category-colour
#' registry; cycle through the lead palette deterministically.
#'
#' @param data Tibble from \code{\link{ontology_terms_by_ontology}}.
#' @param width_mm Numeric.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_identity labs
#' @importFrom ggplot2 scale_y_continuous expansion theme element_text
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @export
plot_ontology_terms_by_ontology <- function(data, width_mm = 89L) {

    ontology <- n_terms <- fill_hex <- NULL

    n_rows <- max(nrow(data), 1L)
    fills <- palette_n(n_rows, unknown = (n_rows >= length(palette_lead())))

    data <- dplyr::mutate(data, fill_hex = fills[seq_len(nrow(data))])

    ggplot2::ggplot(
        data,
        ggplot2::aes(
            x    = stats::reorder(.data$ontology, -.data$n_terms),
            y    = .data$n_terms,
            fill = .data$fill_hex
        )
    ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::scale_fill_identity() +
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
