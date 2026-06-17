#' Canonicalise upstream source names for display
#'
#' The MetaLinksDB v2.0 sources arrive in mixed case and acronym
#' conventions (e.g. \code{"chembl"}, \code{"Stitch"},
#' \code{"guidetopharma"}, \code{"Cellinker"}). This helper maps
#' every observed form to a single canonical display label so the
#' figure legend and axis labels look right.
#'
#' Unknown inputs are returned with the first letter upper-cased.
#'
#' @param x Character vector of raw source names.
#' @return Character vector of pretty display labels.
#'
#' @keywords internal
#' @noRd
pretty_source <- function(x) {
    canonical <- c(
        bindingdb     = "BindingDB",
        cellinker     = "CellInker",
        cellphonedb   = "CellPhoneDB",
        chebi         = "ChEBI",
        chembl        = "ChEMBL",
        drugcentral   = "DrugCentral",
        foodb         = "FooDB",
        guidetopharma = "GuideToPharma",
        hmdb          = "HMDB",
        hmr           = "HMR",
        macdb         = "MACDB",
        metatlas      = "MetAtlas",
        mrclinksdb    = "MRClinksDB",
        neuronchat    = "NeuronChat",
        pfocr         = "PFOCR",
        recon         = "Recon",
        recon3d       = "Recon3D",
        rhea          = "Rhea",
        scconnect     = "scConnect",
        stitch        = "STITCH",
        swisslipids   = "SwissLipids",
        tcdb          = "TCDB"
    )

    raw <- as.character(x)
    key <- tolower(raw)
    hit <- canonical[key]

    out <- ifelse(
        is.na(hit),
        sub("^(.)", "\\U\\1", raw, perl = TRUE),
        unname(hit)
    )
    out
}


#' Shorten the longest metabolite-class display labels
#'
#' A few HMDB-derived metabolite classes have names too long to fit
#' the strip-layout's narrow per-panel width. This helper maps those
#' few specifically-named classes to a publication-friendly short
#' form; anything else is returned unchanged.
#'
#' @param x Character vector of raw metabolite-class names.
#' @return Character vector with the long classes shortened.
#'
#' @keywords internal
#' @noRd
pretty_metabolite_class <- function(x) {
    shortmap <- c(
        "Amino acids, peptides, and analogues"      = "Amino acids",
        "Carbohydrates and carbohydrate conjugates" = "Carbohydrates",
        "Fatty acids and conjugates"                = "Fatty acids"
    )
    raw <- as.character(x)
    hit <- shortmap[raw]
    ifelse(is.na(hit), raw, unname(hit))
}


#' Canonicalise relation-type codes for display
#'
#' The v2 MetaLinksDB exposes four relation-type codes:
#' \code{transport}, \code{interaction}, \code{lr} (ligand-receptor),
#' \code{pd} (pharmacodynamic / "potency depends on"). This helper
#' maps each code to a publication-ready label.
#'
#' Unknown inputs are returned as-is (first letter upper-cased).
#'
#' @param x Character vector of raw relation_type codes.
#' @return Character vector of display labels.
#'
#' @keywords internal
#' @noRd
pretty_relation_type <- function(x) {
    canonical <- c(
        transport   = "Transport",
        interaction = "Interaction",
        lr          = "Ligand-receptor",
        pd          = "Drug-target",
        receptor    = "Receptor"
    )

    raw <- as.character(x)
    key <- tolower(raw)
    hit <- canonical[key]
    out <- ifelse(
        is.na(hit),
        sub("^(.)", "\\U\\1", raw, perl = TRUE),
        unname(hit)
    )
    out
}


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
        ggplot2::labs(x = NULL, y = "Count", fill = NULL) +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(
            plot.title       = ggplot2::element_blank(),
            legend.position  = "top",
            axis.text        = ggplot2::element_text(size = 9),
            axis.title       = ggplot2::element_text(size = 11),
            legend.text      = ggplot2::element_text(size = 8),
            legend.title     = ggplot2::element_text(size = 9),
            legend.key.size  = ggplot2::unit(0.3, "cm"),
            axis.text.x      = ggplot2::element_text(
                angle = 30, hjust = 1
            )
        )
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

    # Shorten the few classes whose names overflow the strip
    # layout's narrow per-panel width, then re-aggregate so
    # duplicates (none expected here, but safe) collapse.
    class_counts$metabolite_class_label <- pretty_metabolite_class(
        class_counts$metabolite_class_label
    )
    class_counts <- dplyr::summarise(
        dplyr::group_by(
            class_counts,
            .data$resource,
            .data$metabolite_class_label
        ),
        n = sum(.data$n),
        .groups = "drop"
    )

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
            x = "Metabolite class",
            y = "Interactions"
        ) +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(
            plot.title  = ggplot2::element_blank(),
            legend.position = "none",
            axis.text   = ggplot2::element_text(size = 9),
            axis.title  = ggplot2::element_text(size = 11)
        )
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

    # Normalise the protein-class strings: strip the
    # ":OM:NNNN" / ":MI:NNNN" ontology-id suffix (e.g.
    # "Gpcr:OM:0040" → "Gpcr"), strip surrounding literal quote
    # characters (some upstream sources emit "gpcr" with quotes),
    # replace underscores with spaces, then apply a canonical
    # display map so both "gpcr" and "Gpcr:OM:0040" land on the
    # single "GPCR" label.
    pretty_protein_class <- function(x) {
        out <- as.character(x)
        # Strip OM-style suffix
        out <- sub("\\s*:[A-Z]+:\\d+\\s*$", "", out)
        # Strip leading/trailing literal " or ' characters
        out <- gsub("^['\"]+|['\"]+$", "", out)
        # Underscores → spaces
        out <- gsub("_", " ", out)
        out <- trimws(out)
        renames <- c(
            "gpcr"                     = "GPCR",
            "vgic"                     = "VGIC",
            "lgic"                     = "LGIC",
            "enzyme"                   = "Enzyme",
            "transporter"              = "Transporter",
            "catalytic receptor"       = "Catalytic receptor",
            "nuclear hormone receptor" = "NHR",
            "nhr"                      = "NHR",
            "other protein"            = "Other",
            "other"                    = "Other"
        )
        hits <- match(tolower(out), names(renames))
        has_rename <- !is.na(hits)
        out[has_rename] <- renames[hits[has_rename]]
        out
    }
    class_counts$protein_class_label <- pretty_protein_class(
        class_counts$protein_class_label
    )
    # Re-aggregate after normalisation so duplicate (resource,
    # class) rows (e.g. "Gpcr:OM:0040" + "gpcr" both → "GPCR")
    # collapse into one bar.
    class_counts <- dplyr::summarise(
        dplyr::group_by(
            class_counts, .data$resource, .data$protein_class_label
        ),
        n = sum(.data$n),
        .groups = "drop"
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
        ggplot2::scale_y_continuous(n.breaks = 3L) +
        ggplot2::coord_flip() +
        ggplot2::labs(
            x = "Protein class",
            y = "Interactions"
        ) +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(
            plot.title  = ggplot2::element_blank(),
            legend.position = "none",
            axis.text   = ggplot2::element_text(size = 9),
            axis.title  = ggplot2::element_text(size = 11)
        )
}


#' Figure 4 Panel D — MetaLinksDB 2.0 overview (FR-010d)
#'
#' MetaLinksDB-2.0-only view: for each upstream source the network was
#' built from (chembl, hmdb, cellinker, …), report \code{Interactions},
#' \code{Metabolites}, and \code{Proteins} as a NORMAL (grouped) bar
#' chart — not stacked — per FR-010d. Counts are de-duplicated within
#' each (source, metric) cell so a metabolite or protein contributed
#' by N relations in the same source counts once.
#'
#' @param data Harmonized MPI rows.
#' @param width_mm Numeric target width.
#' @return A ggplot object.
#' @importFrom dplyr filter group_by summarise n_distinct
#' @importFrom ggplot2 ggplot aes geom_col coord_flip labs
#' @importFrom ggplot2 scale_fill_manual theme position_dodge2
#' @importFrom rlang .data
#' @export
fig03_metalinks_overview_panel <- function(data, width_mm = 180L) {

    resource <- source <- metric <- n <- NULL

    metalinks <- dplyr::filter(
        data, .data$resource == 'MetaLinksDB v2.0'
    )

    if (nrow(metalinks) == 0L) {
        return(empty_fig03_panel(
            title    = 'MetaLinksDB 2.0 overview',
            subtitle = 'No MetaLinksDB v2.0 rows available',
            width_mm = width_mm
        ))
    }

    summary <- dplyr::summarise(
        dplyr::group_by(metalinks, .data$source),
        Interactions = dplyr::n_distinct(
            paste(.data$hmdb_id, .data$uniprot_id, sep = '|')
        ),
        Metabolites  = dplyr::n_distinct(.data$hmdb_id),
        Proteins     = dplyr::n_distinct(.data$uniprot_id),
        .groups      = 'drop'
    )

    long <- data.frame(
        source = rep(pretty_source(summary$source), 3L),
        metric = rep(
            c("Interactions", "Metabolites", "Proteins"),
            each = nrow(summary)
        ),
        n      = c(summary$Interactions, summary$Metabolites,
                   summary$Proteins)
    )
    long$metric <- factor(
        long$metric,
        levels = c("Interactions", "Metabolites", "Proteins")
    )

    ggplot2::ggplot(
        long,
        ggplot2::aes(x = .data$source, y = .data$n, fill = .data$metric)
    ) +
        ggplot2::geom_col(
            position = ggplot2::position_dodge2(preserve = "single"),
            width    = 0.8
        ) +
        ggplot2::scale_fill_manual(values = c(
            Interactions = palette_lead()[["teal"]],
            Metabolites  = palette_lead()[["amber"]],
            Proteins     = palette_lead()[["magenta"]]
        )) +
        ggplot2::scale_y_continuous(n.breaks = 3L) +
        ggplot2::coord_flip() +
        ggplot2::labs(
            x    = "Source",
            y    = "Count",
            fill = NULL
        ) +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(
            plot.title       = ggplot2::element_blank(),
            legend.position  = "top",
            axis.text        = ggplot2::element_text(size = 9),
            axis.title       = ggplot2::element_text(size = 11),
            legend.text      = ggplot2::element_text(size = 8),
            legend.title     = ggplot2::element_text(size = 9),
            legend.key.size  = ggplot2::unit(0.3, "cm")
        )
}


#' Figure 4 Panel E — relationship-type composition (FR-010e)
#'
#' Source × relationship-type grouped bars, restricted to the three
#' key categories \code{transport}, \code{receptor}, and
#' \code{interaction} (FR-010e). The set of sources shown is the
#' union of upstream sources the harmonized MPI rows expose across
#' all included resources.
#'
#' @param data Harmonized MPI rows.
#' @param width_mm Numeric target width.
#' @return A ggplot object.
#' @importFrom dplyr filter count
#' @importFrom ggplot2 ggplot aes geom_col coord_flip labs
#' @importFrom ggplot2 scale_fill_manual theme position_dodge2
#' @importFrom rlang .data
#' @export
fig03_relationship_types_panel <- function(data, width_mm = 180L) {

    source <- relation_type <- n <- NULL

    # All four relation-type codes the prod metalinksdb view exposes
    # (lower-case) — was previously narrowed to (transport, receptor,
    # interaction) but `receptor` does not actually appear and `lr` /
    # `pd` are real categories that belong in the panel.
    kept_types <- c("transport", "interaction", "lr", "pd")

    sub <- dplyr::filter(
        data, tolower(.data$relation_type) %in% kept_types
    )

    if (nrow(sub) == 0L) {
        return(empty_fig03_panel(
            title    = "Relationship types",
            subtitle = "No transport / interaction / lr / pd rows",
            width_mm = width_mm
        ))
    }

    summary <- dplyr::count(
        sub, .data$source, .data$relation_type, name = "n"
    )
    summary$source <- pretty_source(summary$source)
    display_levels <- pretty_relation_type(kept_types)
    summary$relation_type <- factor(
        pretty_relation_type(tolower(summary$relation_type)),
        levels = display_levels
    )

    fill_values <- stats::setNames(
        palette_n(as.integer(length(display_levels)), unknown = FALSE),
        display_levels
    )

    ggplot2::ggplot(
        summary,
        ggplot2::aes(
            x    = .data$source,
            y    = .data$n,
            fill = .data$relation_type
        )
    ) +
        ggplot2::geom_col(
            position = ggplot2::position_dodge2(preserve = "single"),
            width    = 0.8
        ) +
        ggplot2::scale_fill_manual(values = fill_values) +
        ggplot2::scale_y_continuous(n.breaks = 3L) +
        ggplot2::coord_flip() +
        ggplot2::labs(
            x    = "Source",
            y    = "Interactions",
            fill = NULL
        ) +
        # 4 relation types with names up to 15 chars don't fit on a
        # single row in the strip-layout's narrow per-panel width;
        # wrap to a 2x2 legend block so the "Pharmacodynamic" entry
        # no longer overflows the right edge of the panel.
        ggplot2::guides(
            fill = ggplot2::guide_legend(nrow = 2L, byrow = TRUE)
        ) +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(
            plot.title       = ggplot2::element_blank(),
            legend.position  = "top",
            axis.text        = ggplot2::element_text(size = 9),
            axis.title       = ggplot2::element_text(size = 11),
            legend.text      = ggplot2::element_text(size = 8),
            legend.title     = ggplot2::element_text(size = 9),
            legend.key.size  = ggplot2::unit(0.3, "cm")
        )
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
#' @importFrom dplyr mutate bind_rows select
#' @importFrom tidyr complete
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
    plot_data <- tidyr::complete(
        plot_data,
        interaction_type, panel_group,
        fill = list(n_interactions = 0)
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
        ggplot2::theme(
            legend.position = "top",
            plot.title      = ggplot2::element_text(size = 14),
            axis.text       = ggplot2::element_text(size = 11),
            axis.title      = ggplot2::element_text(size = 12),
            legend.text     = ggplot2::element_text(size = 10)
        )
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
            x    = "Compartment",
            y    = "Interactions",
            fill = NULL
        ) +
        ggplot2::guides(
            fill = ggplot2::guide_legend(nrow = 3L, ncol = 2L)
        ) +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(
            legend.position  = "top",
            plot.title       = ggplot2::element_blank(),
            axis.text        = ggplot2::element_text(size = 9),
            axis.title       = ggplot2::element_text(size = 12),
            legend.text      = ggplot2::element_text(size = 7),
            legend.title     = ggplot2::element_text(size = 8),
            legend.key.size  = ggplot2::unit(0.3, "cm"),
            legend.spacing.y = ggplot2::unit(0.1, "cm")
        )

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
#' @param position Character: bar arrangement. \code{"stack"} (default)
#'     keeps the existing stacked layout; \code{"dodge"} draws
#'     Metabolites and Proteins as side-by-side bars per resource.
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual labs
#'     position_stack position_dodge coord_flip
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate bind_rows slice_head
#' @importFrom tibble tibble
#' @importFrom rlang .data abort
#' @export
fig04_resource_contribution_panel <- function(
    cosmos_plus_by_resource,
    top_n    = 15L,
    width_mm = 89L,
    position = c("stack", "dodge")
) {

    position <- match.arg(position)
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

    # Sort descending by total entity count; pin "Other" to bottom
    data <- dplyr::bind_rows(
        dplyr::arrange(
            data[data$resource != "Other", ],
            dplyr::desc(n_metabolites + n_proteins)
        ),
        data[data$resource == "Other", ]
    )

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

    # Vivid pair of similar intensity (lead palette) — avoiding the
    # reserved grey #BEBEBE which is the canonical None / Unknown
    # slot (FR-020). Teal vs. magenta gives two saturated, similarly
    # chromatic hues that read well at small panel sizes.
    lead <- palette_lead()
    entity_fills <- c(
        Metabolites = unname(lead[["teal"]]),
        Proteins    = unname(lead[["magenta"]])
    )

    geom_position <- switch(
        position,
        stack = ggplot2::position_stack(),
        dodge = ggplot2::position_dodge(width = 0.75)
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
            position = geom_position,
            width    = 0.7
        ) +
        ggplot2::scale_fill_manual(values = entity_fills) +
        ggplot2::coord_flip() +
        ggplot2::labs(
            x    = "Resource",
            y    = "Unique entities",
            fill = NULL
        ) +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(
            legend.position = "top",
            plot.title      = ggplot2::element_blank(),
            axis.text       = ggplot2::element_text(size = 9),
            axis.title      = ggplot2::element_text(size = 12),
            legend.text     = ggplot2::element_text(size = 7),
            legend.title    = ggplot2::element_text(size = 8),
            legend.key.size = ggplot2::unit(0.3, "cm")
        )
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
#' @importFrom dplyr group_by summarise mutate bind_rows case_match ungroup
#' @importFrom tidyr complete
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

    combined <- dplyr::bind_rows(cosmos_agg, metalinks_long) |>
        tidyr::complete(
            interaction_type, source,
            fill = list(n_interactions = 0L)
        )

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
        ggplot2::theme(
            legend.position = "top",
            plot.title      = ggplot2::element_text(size = 14),
            axis.text       = ggplot2::element_text(size = 11),
            axis.title      = ggplot2::element_text(size = 12),
            legend.text     = ggplot2::element_text(size = 10)
        )
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
