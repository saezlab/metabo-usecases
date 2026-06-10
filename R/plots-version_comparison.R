#' Build a placeholder Figure 3 panel when data are unavailable
#'
#' @param title Character scalar.
#' @param subtitle Character scalar.
#' @param width_mm Numeric.
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes geom_col labs scale_fill_identity theme
#' @importFrom tibble tibble
#' @keywords internal
#' @noRd
empty_fig03_panel <- function(title, subtitle, width_mm = 89L) {
    ggplot2::ggplot(
        tibble::tibble(label = subtitle, value = 0, fill_hex = '#BEBEBE'),
        ggplot2::aes(x = label, y = value, fill = fill_hex)
    ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::scale_fill_identity() +
        ggplot2::labs(title = title, x = NULL, y = NULL) +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(legend.position = 'none')
}


#' Figure 3A -- overall resource coverage
#'
#' @param data Harmonized MPI rows.
#' @param width_mm Numeric target width.
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes geom_col labs scale_fill_manual theme
#' @importFrom rlang .data
#' @export
fig03_coverage_panel <- function(data, width_mm = 180L) {

    resource <- metric <- n <- NULL

    summary <- summarize_mpi_coverage(data)
    if (nrow(summary) == 0L) {
        return(empty_fig03_panel(
            title = 'Coverage',
            subtitle = 'No comparable resources available',
            width_mm = width_mm
        ))
    }

    long <- data.frame(
        resource = rep(summary$resource, 3L),
        metric = rep(c('Interactions', 'Metabolites', 'Proteins'), each = nrow(summary)),
        n = c(summary$n_interactions, summary$n_metabolites, summary$n_proteins)
    )

    ggplot2::ggplot(
        long,
        ggplot2::aes(x = .data$resource, y = .data$n, fill = .data$metric)
    ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::scale_fill_manual(values = c(
            Interactions = palette_lead()[['teal']],
            Metabolites = palette_lead()[['amber']],
            Proteins = palette_lead()[['magenta']]
        )) +
        ggplot2::labs(x = NULL, y = 'Count', fill = NULL, title = 'Coverage') +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(legend.position = 'top')
}


#' Figure 3B -- metabolite-class coverage
#'
#' @param data Harmonized MPI rows.
#' @param width_mm Numeric target width.
#' @param top_n Integer number of classes to keep.
#' @return A ggplot object.
#' @importFrom dplyr count group_by mutate slice_max summarise
#' @importFrom ggplot2 ggplot aes geom_col coord_flip labs scale_fill_identity theme
#' @importFrom rlang .data
#' @export
fig03_metabolite_class_panel <- function(data,
                                         width_mm = 180L,
                                         top_n = 10L) {

    resource <- metabolite_class_label <- n <- fill_hex <- total_n <- NULL

    class_counts <- dplyr::count(
        data,
        .data$resource,
        .data$metabolite_class_label,
        name = 'n'
    )

    if (nrow(class_counts) == 0L) {
        return(empty_fig03_panel(
            title = 'Metabolite classes',
            subtitle = 'No metabolite-class annotations available',
            width_mm = width_mm
        ))
    }

    keep <- dplyr::summarise(
        dplyr::group_by(class_counts, .data$metabolite_class_label),
        total_n = sum(.data$n),
        .groups = 'drop'
    )
    keep <- dplyr::slice_max(keep, .data$total_n, n = top_n, with_ties = FALSE)

    class_counts <- class_counts[class_counts$metabolite_class_label %in% keep$metabolite_class_label, ]

    fills <- setNames(
        palette_n(as.integer(length(unique(class_counts$resource))), unknown = FALSE),
        sort(unique(class_counts$resource))
    )
    register_category_colours('resources', fills)

    ggplot2::ggplot(
        dplyr::mutate(
            class_counts,
            fill_hex = category_colour('resources', .data$resource)
        ),
        ggplot2::aes(
            x = .data$metabolite_class_label,
            y = .data$n,
            fill = .data$fill_hex
        )
    ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::scale_fill_identity() +
        ggplot2::coord_flip() +
        ggplot2::labs(
            x = NULL,
            y = 'Interactions',
            title = 'Metabolite-class breadth'
        ) +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(legend.position = 'none')
}


#' Figure 3C -- protein-class coverage
#'
#' @param data Harmonized MPI rows.
#' @param width_mm Numeric target width.
#' @param top_n Integer number of classes to keep.
#' @return A ggplot object.
#' @importFrom dplyr count group_by mutate slice_max summarise
#' @importFrom ggplot2 ggplot aes geom_col coord_flip labs scale_fill_identity theme
#' @importFrom rlang .data
#' @export
fig03_protein_class_panel <- function(data,
                                      width_mm = 180L,
                                      top_n = 10L) {

    resource <- protein_class_label <- n <- fill_hex <- total_n <- NULL

    class_counts <- dplyr::count(
        data,
        .data$resource,
        .data$protein_class_label,
        name = 'n'
    )

    if (nrow(class_counts) == 0L) {
        return(empty_fig03_panel(
            title = 'Protein classes',
            subtitle = 'No protein-class annotations available',
            width_mm = width_mm
        ))
    }

    keep <- dplyr::summarise(
        dplyr::group_by(class_counts, .data$protein_class_label),
        total_n = sum(.data$n),
        .groups = 'drop'
    )
    keep <- dplyr::slice_max(keep, .data$total_n, n = top_n, with_ties = FALSE)

    class_counts <- class_counts[class_counts$protein_class_label %in% keep$protein_class_label, ]
    fills <- setNames(
        palette_n(as.integer(length(unique(class_counts$resource))), unknown = FALSE),
        sort(unique(class_counts$resource))
    )
    register_category_colours('resources', fills)

    ggplot2::ggplot(
        dplyr::mutate(
            class_counts,
            fill_hex = category_colour('resources', .data$resource)
        ),
        ggplot2::aes(
            x = .data$protein_class_label,
            y = .data$n,
            fill = .data$fill_hex
        )
    ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::scale_fill_identity() +
        ggplot2::coord_flip() +
        ggplot2::labs(
            x = NULL,
            y = 'Interactions',
            title = 'Protein-class breadth'
        ) +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(legend.position = 'none')
}


#' Figure 3D -- evidence / confidence availability
#'
#' @param data Harmonized MPI rows.
#' @param width_mm Numeric target width.
#' @return A ggplot object.
#' @importFrom dplyr bind_rows group_by summarise
#' @importFrom ggplot2 ggplot aes geom_col labs scale_fill_manual theme
#' @importFrom rlang .data
#' @export
fig03_evidence_confidence_panel <- function(data, width_mm = 180L) {

    resource <- dimension <- n <- NULL

    summary <- dplyr::bind_rows(
        dplyr::summarise(dplyr::group_by(data, .data$resource), dimension = 'Multi-source', n = sum(.data$source_count > 1, na.rm = TRUE), .groups = 'drop'),
        dplyr::summarise(dplyr::group_by(data, .data$resource), dimension = 'Citations', n = sum(!is.na(.data$citation_count) & .data$citation_count > 0, na.rm = TRUE), .groups = 'drop'),
        dplyr::summarise(dplyr::group_by(data, .data$resource), dimension = 'Affinity', n = sum(!is.na(.data$affinity_value), na.rm = TRUE), .groups = 'drop'),
        dplyr::summarise(dplyr::group_by(data, .data$resource), dimension = 'Curation mode', n = sum(!is.na(.data$curation_mode) & .data$curation_mode != '', na.rm = TRUE), .groups = 'drop')
    )

    if (nrow(summary) == 0L || all(summary$n == 0L)) {
        return(empty_fig03_panel(
            title = 'Evidence / confidence',
            subtitle = 'Only source-count evidence available so far',
            width_mm = width_mm
        ))
    }

    ggplot2::ggplot(
        summary,
        ggplot2::aes(x = .data$resource, y = .data$n, fill = .data$dimension)
    ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::scale_fill_manual(values = c(
            `Multi-source` = palette_lead()[['teal']],
            Citations = palette_lead()[['amber']],
            Affinity = palette_lead()[['magenta']],
            `Curation mode` = palette_lead()[['green']]
        )) +
        ggplot2::labs(x = NULL, y = 'Interactions', fill = NULL, title = 'Evidence / confidence') +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(legend.position = 'top')
}


#' Figure 3E -- source by relationship type
#'
#' @param data Harmonized MPI rows.
#' @param width_mm Numeric target width.
#' @return A ggplot object.
#' @importFrom dplyr count
#' @importFrom ggplot2 ggplot aes geom_col coord_flip labs scale_fill_manual theme
#' @importFrom rlang .data
#' @export
fig03_source_relationship_panel <- function(data, width_mm = 180L) {

    source <- relation_type <- n <- NULL

    summary <- dplyr::count(data, .data$source, .data$relation_type, name = 'n')
    if (nrow(summary) == 0L) {
        return(empty_fig03_panel(
            title = 'Source × relationship type',
            subtitle = 'No source-composition data available',
            width_mm = width_mm
        ))
    }

    relation_types <- sort(unique(summary$relation_type))
    fills <- setNames(
        palette_n(as.integer(length(relation_types)), unknown = FALSE),
        relation_types
    )
    register_category_colours('interaction_types', fills)

    ggplot2::ggplot(
        summary,
        ggplot2::aes(x = .data$source, y = .data$n, fill = .data$relation_type)
    ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::scale_fill_manual(values = fills) +
        ggplot2::coord_flip() +
        ggplot2::labs(
            x = NULL,
            y = 'Interactions',
            fill = NULL,
            title = 'Source × relationship type'
        ) +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(legend.position = 'top')
}
