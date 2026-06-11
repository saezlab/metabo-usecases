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


# ---------------------------------------------------------------------------
# Figure 4 renderers (T053)
# ---------------------------------------------------------------------------

#' Figure 4 Panel A: species-aware COSMOS old vs. COSMOS+ comparison
#'
#' Grouped-bar plot with three bars per interaction-type group:
#' old-COSMOS (human), COSMOS+ human, COSMOS+ mouse.
#'
#' @param old_pkn_tally Tibble from \code{\link{cosmos_old_pkn}} with
#'     columns \code{interaction_type}, \code{n_edges}.
#' @param cosmos_plus_by_type_species Tibble from
#'     \code{cosmos_plus_data()$by_type_species} with columns
#'     \code{interaction_type}, \code{species}, \code{n_interactions}.
#' @param width_mm Numeric panel width in mm.
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual labs
#'     position_dodge coord_flip
#' @importFrom dplyr mutate bind_rows
#' @importFrom tibble tibble
#' @importFrom rlang .data
#' @export
fig04_cosmos_comparison_panel <- function(
    old_pkn_tally,
    cosmos_plus_by_type_species,
    width_mm = 180L
) {
    interaction_type <- n_interactions <- panel_group <- NULL

    old_rows <- dplyr::mutate(
        old_pkn_tally,
        panel_group    = "Old COSMOS (human)",
        n_interactions = .data$n_edges
    )

    new_rows <- dplyr::mutate(
        cosmos_plus_by_type_species,
        panel_group = dplyr::case_when(
            species == "human" ~ "COSMOS+ (human)",
            species == "mouse" ~ "COSMOS+ (mouse)",
            .default           = paste0("COSMOS+ (", species, ")")
        )
    )

    all_types <- sort(unique(c(
        old_rows$interaction_type,
        new_rows$interaction_type
    )))

    groups <- c("Old COSMOS (human)", "COSMOS+ (human)", "COSMOS+ (mouse)")
    fills  <- setNames(
        palette_n(as.integer(length(groups)), unknown = FALSE),
        groups
    )

    plot_data <- dplyr::bind_rows(
        dplyr::select(
            old_rows,
            interaction_type, n_interactions, panel_group
        ),
        dplyr::select(
            new_rows,
            interaction_type, n_interactions, panel_group
        )
    )
    plot_data$panel_group <- factor(
        plot_data$panel_group,
        levels = groups
    )
    plot_data$interaction_type <- factor(
        plot_data$interaction_type,
        levels = rev(all_types)
    )

    ggplot2::ggplot(
        plot_data,
        ggplot2::aes(
            x    = .data$interaction_type,
            y    = .data$n_interactions,
            fill = .data$panel_group
        )
    ) +
        ggplot2::geom_col(
            position = ggplot2::position_dodge(0.8),
            width    = 0.7
        ) +
        ggplot2::scale_fill_manual(values = fills) +
        ggplot2::coord_flip() +
        ggplot2::labs(
            x     = NULL,
            y     = "Edges",
            fill  = NULL,
            title = "Old COSMOS vs. COSMOS+ by interaction type"
        ) +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(legend.position = "top")
}


#' Figure 4 Panel B: COSMOS+ interactions per compartment
#'
#' @param cosmos_plus_by_compartment Tibble from
#'     \code{cosmos_plus_data()$by_compartment}.
#' @param top_n Integer: keep the top N compartments, collapsing the
#'     rest to "Other" per FR-009c. Default 12.
#' @param width_mm Numeric panel width in mm.
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes geom_col labs scale_fill_manual
#' @importFrom dplyr mutate if_else slice_head summarise
#' @importFrom rlang .data
#' @export
fig04_compartment_panel <- function(
    cosmos_plus_by_compartment,
    top_n    = 12L,
    width_mm = 89L
) {
    compartment <- n_interactions <- NULL

    data <- cosmos_plus_by_compartment

    if (nrow(data) > top_n) {
        top    <- data[seq_len(top_n), ]
        others <- sum(data$n_interactions[seq(top_n + 1L, nrow(data))])
        other_row <- tibble::tibble(
            compartment    = "Other",
            n_interactions = others
        )
        data <- dplyr::bind_rows(top, other_row)
    }

    data$compartment <- factor(
        data$compartment,
        levels = rev(data$compartment)
    )

    unannotated <- identical(data$compartment[[1L]], "unannotated")

    fill_col <- if (unannotated) "#BEBEBE" else palette_n(1L, unknown = FALSE)

    p <- ggplot2::ggplot(
        data,
        ggplot2::aes(
            x = .data$compartment,
            y = .data$n_interactions
        )
    ) +
        ggplot2::geom_col(fill = fill_col[[1L]], width = 0.7) +
        ggplot2::coord_flip() +
        ggplot2::labs(
            x     = NULL,
            y     = "Edges",
            title = "COSMOS+ edges per compartment"
        ) +
        theme_bw_metabo(width_mm = width_mm)

    if (unannotated) {
        p <- p + ggplot2::labs(
            subtitle = "Location annotations not yet available for this build"
        )
    }

    p
}


#' Figure 4 Panel C: entity and interaction count per resource
#'
#' @param cosmos_plus_by_resource Tibble from
#'     \code{cosmos_plus_data()$by_resource}.
#' @param width_mm Numeric panel width in mm.
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual labs
#'     position_stack
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @export
fig04_resource_contribution_panel <- function(
    cosmos_plus_by_resource,
    width_mm = 89L
) {
    resource <- entity_type <- count <- NULL

    data <- cosmos_plus_by_resource
    resources <- data$resource

    fills <- setNames(
        palette_n(as.integer(length(resources)), unknown = FALSE),
        resources
    )
    register_category_colours("resources", fills)

    long <- tidyr::pivot_longer(
        data,
        cols      = c("n_metabolites", "n_proteins"),
        names_to  = "entity_type",
        values_to = "count"
    )
    long$entity_type <- dplyr::case_when(
        long$entity_type == "n_metabolites" ~ "Metabolites",
        long$entity_type == "n_proteins"    ~ "Proteins",
        .default = long$entity_type
    )
    long$resource <- factor(long$resource, levels = rev(resources))

    entity_fills <- c(
        Metabolites = palette_n(1L)[[1L]],
        Proteins    = palette_n(2L)[[2L]]
    )

    ggplot2::ggplot(
        long,
        ggplot2::aes(
            x    = .data$resource,
            y    = .data$count,
            fill = .data$entity_type
        )
    ) +
        ggplot2::geom_col(
            position = ggplot2::position_stack(),
            width    = 0.7
        ) +
        ggplot2::scale_fill_manual(values = entity_fills) +
        ggplot2::coord_flip() +
        ggplot2::labs(
            x     = NULL,
            y     = "Unique entities",
            fill  = NULL,
            title = "COSMOS+ entities per resource"
        ) +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(legend.position = "top")
}
