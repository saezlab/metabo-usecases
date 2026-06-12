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

    # Alias catalysis → metabolic_reactions so COSMOS+ GEM interactions
    # align with the old PKN's metabolic_reactions bar group.
    new_rows <- dplyr::mutate(
        cosmos_plus_by_type_species,
        interaction_type = dplyr::case_match(
            interaction_type,
            "catalysis" ~ "metabolic_reactions",
            .default    = interaction_type
        ),
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


#' Figure 4 Panel B: COSMOS+ interactions per compartment, stacked by type
#'
#' @param cosmos_plus_by_compartment Tibble from
#'     \code{cosmos_plus_data()$by_compartment} with columns
#'     \code{compartment_name}, \code{interaction_type},
#'     \code{n_interactions}.
#' @param top_n Integer: keep the top N compartments by total count,
#'     collapsing the rest to "Other" per FR-009c. Default 12.
#' @param width_mm Numeric panel width in mm.
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes geom_col labs scale_fill_manual
#'     position_stack coord_flip
#' @importFrom dplyr mutate filter group_by summarise bind_rows arrange desc
#' @importFrom tibble tibble
#' @importFrom rlang .data
#' @export
fig04_compartment_panel <- function(
    cosmos_plus_by_compartment,
    top_n    = 12L,
    width_mm = 89L
) {
    compartment_name <- interaction_type <- n_interactions <- total <- NULL

    type_labels <- c(
        catalysis             = "Metabolic reaction (enzyme–metabolite)",
        transport             = "Transport (transporter–metabolite)",
        gene_regulation       = "Gene regulation (TF–target, GRN)",
        signaling             = "Signaling (PPI)",
        allosteric_regulation = "Allosteric regulation (metabolite–enzyme)",
        ligand_receptor       = "Ligand receptor (receptor–metabolite)"
    )

    data <- cosmos_plus_by_compartment

    # Compute per-compartment totals for ordering and top-N collapsing
    totals <- data |>
        dplyr::group_by(compartment_name) |>
        dplyr::summarise(total = sum(n_interactions), .groups = "drop") |>
        dplyr::arrange(dplyr::desc(total))

    unannotated_only <- all(
        data$compartment_name %in% c("Unannotated", NA_character_)
    )

    if (nrow(totals) > top_n) {
        keep <- totals$compartment_name[seq_len(top_n)]
        other_rows <- data |>
            dplyr::filter(!compartment_name %in% keep) |>
            dplyr::group_by(interaction_type) |>
            dplyr::summarise(
                n_interactions = sum(n_interactions), .groups = "drop"
            ) |>
            dplyr::mutate(compartment_name = "Other")
        data <- dplyr::bind_rows(
            data |> dplyr::filter(compartment_name %in% keep),
            other_rows
        )
        other_total <- sum(totals$total[seq(top_n + 1L, nrow(totals))])
        totals <- dplyr::bind_rows(
            totals[seq_len(top_n), ],
            tibble::tibble(compartment_name = "Other", total = other_total)
        )
    }

    comp_levels <- rev(totals$compartment_name)
    data$compartment_name <- factor(data$compartment_name, levels = comp_levels)

    # Apply display labels to interaction_type
    data$interaction_type <- dplyr::case_match(
        data$interaction_type,
        "catalysis"             ~ type_labels[["catalysis"]],
        "transport"             ~ type_labels[["transport"]],
        "gene_regulation"       ~ type_labels[["gene_regulation"]],
        "signaling"             ~ type_labels[["signaling"]],
        "allosteric_regulation" ~ type_labels[["allosteric_regulation"]],
        "ligand_receptor"       ~ type_labels[["ligand_receptor"]],
        .default = data$interaction_type
    )

    type_vals <- unique(data$interaction_type[!is.na(data$interaction_type)])
    fills <- setNames(
        palette_n(as.integer(length(type_vals)), unknown = FALSE),
        type_vals
    )
    if (unannotated_only) fills["Unannotated"] <- "#BEBEBE"

    p <- ggplot2::ggplot(
        data[!is.na(data$interaction_type), ],
        ggplot2::aes(
            x    = .data$compartment_name,
            y    = .data$n_interactions,
            fill = .data$interaction_type
        )
    ) +
        ggplot2::geom_col(
            position = ggplot2::position_stack(),
            width    = 0.7
        ) +
        ggplot2::scale_fill_manual(values = fills) +
        ggplot2::coord_flip() +
        ggplot2::labs(
            x     = NULL,
            y     = "Interactions",
            fill  = NULL,
            title = "COSMOS+ interactions per compartment"
        ) +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(legend.position = "right")

    if (unannotated_only) {
        p <- p + ggplot2::labs(
            subtitle = "Location annotations not yet available for this build"
        )
    }

    p
}


#' Figure 4 Panel C: entity and interaction count per resource
#'
#' Shows the top 15 resources by total interaction count as stacked bars
#' (metabolites + proteins). Resource names are abbreviated via
#' \code{.abbreviate_resource()} before plotting. Remaining resources
#' beyond top 15 are collapsed into a single "Other" bar.
#'
#' @param cosmos_plus_by_resource Tibble from
#'     \code{cosmos_plus_data()$by_resource}.
#' @param top_n Integer: number of individual resources to show.
#'     Default 15.
#' @param width_mm Numeric panel width in mm.
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual labs
#'     position_stack coord_flip
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate bind_rows slice_head
#' @importFrom tibble tibble
#' @importFrom rlang .data
#' @export
fig04_resource_contribution_panel <- function(
    cosmos_plus_by_resource,
    top_n    = 15L,
    width_mm = 89L
) {
    resource <- entity_type <- count <- n_metabolites <- n_proteins <- NULL
    n_interactions <- NULL

    data <- cosmos_plus_by_resource

    # Collapse tail into "Other"
    if (nrow(data) > top_n) {
        top   <- data[seq_len(top_n), ]
        other <- data[seq(top_n + 1L, nrow(data)), ]
        other_row <- tibble::tibble(
            resource       = "Other",
            n_metabolites  = sum(other$n_metabolites),
            n_proteins     = sum(other$n_proteins),
            n_interactions = sum(other$n_interactions)
        )
        data <- dplyr::bind_rows(top, other_row)
    }

    # Apply short abbreviations
    data$resource_label <- .abbreviate_resource(data$resource)

    resources_ordered <- data$resource_label
    long <- tidyr::pivot_longer(
        data,
        cols      = c("n_metabolites", "n_proteins"),
        names_to  = "entity_type",
        values_to = "count"
    )
    long$entity_type <- dplyr::case_match(
        long$entity_type,
        "n_metabolites" ~ "Metabolites",
        "n_proteins"    ~ "Proteins",
        .default        = long$entity_type
    )
    long$resource_label <- factor(
        long$resource_label,
        levels = rev(resources_ordered)
    )

    entity_fills <- c(
        Metabolites = palette_n(1L)[[1L]],
        Proteins    = palette_n(2L)[[2L]]
    )

    ggplot2::ggplot(
        long,
        ggplot2::aes(
            x    = .data$resource_label,
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


#' Figure 4 Panel D: MetaLinksDB 2.0 vs. COSMOS+ by interaction type
#'
#' Compares MetaLinksDB 2.0 (GtP-class-mapped categories from the Panel D
#' SQL) with COSMOS+ (internal interaction types mapped to the same
#' canonical categories). An "Allosteric regulation" row with zero
#' MetaLinksDB count is appended to make the COSMOS+-only scope visible.
#'
#' @param metalinks_counts Tibble with columns \code{interaction_type} and
#'     \code{n_interactions} from \code{pg_query_panel()}, where
#'     \code{interaction_type} already uses the canonical panel-D category
#'     labels (Transport, Ligand receptor, Catalysis, Gene regulation, Other).
#' @param cosmos_plus_by_type_species Tibble from
#'     \code{cosmos_plus_data()$by_type_species}.
#' @param width_mm Numeric panel width in mm.
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual labs
#'     coord_flip position_dodge
#' @importFrom dplyr group_by summarise mutate bind_rows case_match
#' @importFrom tibble tibble
#' @importFrom rlang .data
#' @export
fig04_metalinks_cosmos_panel <- function(
    metalinks_counts,
    cosmos_plus_by_type_species,
    width_mm = 120L
) {
    interaction_type <- n_interactions <- source <- NULL

    # Map COSMOS+ internal types to canonical panel-D categories
    cosmos_agg <- dplyr::group_by(
        cosmos_plus_by_type_species,
        .data$interaction_type
    ) |>
        dplyr::summarise(
            n_interactions = sum(.data$n_interactions),
            .groups        = "drop"
        ) |>
        dplyr::mutate(
            interaction_type = dplyr::case_match(
                interaction_type,
                "catalysis"             ~ "Catalysis",
                "transport"             ~ "Transport",
                "ligand_receptor"       ~ "Ligand receptor",
                "gene_regulation"       ~ "Gene regulation",
                "allosteric_regulation" ~ "Allosteric regulation",
                "signaling"             ~ "Signaling",
                .default = interaction_type
            ),
            source = "COSMOS+"
        )

    metalinks_long <- dplyr::mutate(metalinks_counts, source = "MetaLinksDB 2.0")

    # Add zero rows for COSMOS+-only categories not present in MetaLinksDB
    cosmos_types    <- unique(cosmos_agg$interaction_type)
    metalinks_types <- unique(metalinks_long$interaction_type)
    cosmos_only     <- setdiff(cosmos_types, metalinks_types)

    if (length(cosmos_only) > 0L) {
        zero_rows <- tibble::tibble(
            interaction_type = cosmos_only,
            n_interactions   = 0L,
            source           = "MetaLinksDB 2.0"
        )
        metalinks_long <- dplyr::bind_rows(metalinks_long, zero_rows)
    }

    combined <- dplyr::bind_rows(cosmos_agg, metalinks_long)

    all_types <- sort(unique(combined$interaction_type))
    combined$interaction_type <- factor(
        combined$interaction_type, levels = rev(all_types)
    )
    combined$source <- factor(
        combined$source, levels = c("MetaLinksDB 2.0", "COSMOS+")
    )

    fills <- c(
        "COSMOS+"         = palette_lead()[["teal"]],
        "MetaLinksDB 2.0" = palette_lead()[["amber"]]
    )
    register_category_colours("metalinks_cosmos", fills)

    ggplot2::ggplot(
        combined,
        ggplot2::aes(
            x    = .data$interaction_type,
            y    = .data$n_interactions,
            fill = .data$source
        )
    ) +
        ggplot2::geom_col(
            position = ggplot2::position_dodge(width = 0.8),
            width    = 0.7
        ) +
        ggplot2::scale_fill_manual(values = fills) +
        ggplot2::coord_flip() +
        ggplot2::labs(
            x     = NULL,
            y     = "Interactions",
            fill  = NULL,
            title = "MetaLinksDB 2.0 vs. COSMOS+"
        ) +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(legend.position = "top")
}


# ── Internal helpers ──────────────────────────────────────────────────────────

# Canonical abbreviation rules for COSMOS+ resource name strings.
.resource_abbrev_single <- c(
    "GEM_transporter:Human-GEM"     = "hGEM-T",
    "GEM_transporter:Mouse-GEM"     = "mGEM-T",
    "GEM:Human-GEM"                 = "hGEM",
    "GEM:Mouse-GEM"                 = "mGEM",
    "GEM:Recon3D"                   = "R3D",
    "OmniPath:omnipath,ligrecextra" = "OmniPath-LR",
    "OmniPath:collectri"            = "Collectri",
    "MRCLinksDB"                    = "MRCLinks"
)

.abbreviate_resource <- function(x) {
    vapply(x, function(name) {
        parts  <- strsplit(name, ";", fixed = TRUE)[[1L]]
        abbrevs <- ifelse(
            parts %in% names(.resource_abbrev_single),
            unname(.resource_abbrev_single[parts]),
            parts
        )
        paste(abbrevs, collapse = "+")
    }, character(1L), USE.NAMES = FALSE)
}
