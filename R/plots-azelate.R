#' Azelate upstream-source palette mapping
#'
#' Reproduces the \code{SOURCE_COLORS} dict from
#' \code{omnipath_metabo_case1/azelate_fig_tables/notebooks/azelaic_acid_plots.py}
#' so the R panels render the same source in the same colour as the
#' Python reference figures.
#'
#' @return Named character vector keyed by source name.
#'
#' @keywords internal
#' @noRd
azelate_source_palette <- function() {
    c(
        chembl      = "#31688E",
        drugcentral = "#F28E2B",
        foodb       = "#59A14F",
        hmdb        = "#E15759",
        macdb       = "#8F63B8",
        metatlas    = "#9C755F",
        pfocr       = "#D66AB0",
        rhea        = "#4E79A7",
        stitch      = "#7F7F7F",
        chebi       = "#B07AA1",
        swisslipids = "#76B7B2"
    )
}


#' Azelate sample-type (tissue) palette mapping
#'
#' Reproduces the \code{TISSUE_COLORS} dict from the azelate
#' Python plotting reference.
#'
#' @return Named character vector keyed by sample-type name.
#'
#' @keywords internal
#' @noRd
azelate_sample_type_palette <- function() {
    c(
        Blood  = "#C6403D",
        Serum  = "#F28E2B",
        Plasma = "#8C564B",
        Tissue = "#59A14F",
        Stool  = "#4E79A7",
        Urine  = "#9467BD",
        Sweat  = "#7F7F7F",
        Cells  = "#E377C2"
    )
}


#' Register azelate source + sample-type colour mappings
#'
#' Merges the azelate palettes into the shared category registry
#' (FR-022 / SC-004). Sources land under the existing
#' \code{"resources"} category; sample types land under a new
#' \code{"sample_types"} category. Idempotent — re-registration
#' with the same mapping is a no-op.
#'
#' @return Invisibly a named list with the two registered mappings.
#'
#' @examples
#' \dontrun{
#' register_azelate_colours()
#' category_colour("resources", "chembl")
#' category_colour("sample_types", "Blood")
#' }
#'
#' @export
register_azelate_colours <- function() {

    src_map  <- azelate_source_palette()
    samp_map <- azelate_sample_type_palette()

    register_category_colours("resources", src_map)

    if (!exists("sample_types", envir = .metabo_category_colours)) {
        .metabo_category_colours$sample_types <- character(0)
    }
    register_category_colours("sample_types", samp_map)

    invisible(list(resources = src_map, sample_types = samp_map))
}


#' Figure 6 Panel C -- interaction type composition for resolved Azelate
#'
#' Stacked-bar plot of relation counts per \code{interaction_type},
#' fill-coloured by upstream \code{source}. Reproduces the
#' \code{make_fig3_panel_C_interaction_type_barplot} reference from
#' the azelate Python notebook through the shared style module.
#'
#' Bars are ordered by descending total relation count; the source
#' stack order within each bar follows descending total contribution
#' across the panel so the dominant source sits at the bottom.
#'
#' @param data Tibble from \code{\link{azelate_panel_c_data}}.
#' @param width_mm Numeric: target physical panel width in mm.
#' @param font_scale Numeric: passed to \code{\link{theme_bw_metabo}}.
#' @param legend_scale Numeric or \code{NULL}: passed to
#'     \code{\link{theme_bw_metabo}}.
#'
#' @return A ggplot object.
#'
#' @importFrom dplyr group_by summarise arrange desc mutate
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_identity
#' @importFrom ggplot2 labs theme element_text
#' @importFrom rlang .data
#' @export
azelate_interaction_panel <- function(
    data,
    width_mm = 89L,
    font_scale = 1,
    legend_scale = NULL
) {

    # NSE vs. R CMD check workaround
    interaction_type <- source <- relation_count <- total <- NULL
    fill_hex <- src_total <- NULL

    bar_order <- dplyr::arrange(
        dplyr::summarise(
            dplyr::group_by(data, .data$interaction_type),
            total = sum(.data$relation_count),
            .groups = "drop"
        ),
        dplyr::desc(total)
    )$interaction_type

    src_order <- dplyr::arrange(
        dplyr::summarise(
            dplyr::group_by(data, .data$source),
            src_total = sum(.data$relation_count),
            .groups = "drop"
        ),
        dplyr::desc(src_total)
    )$source

    plot_data <- dplyr::mutate(
        data,
        interaction_type = factor(
            .data$interaction_type, levels = bar_order
        ),
        source           = factor(.data$source, levels = src_order),
        fill_hex         = category_colour(
            "resources", as.character(.data$source)
        )
    )

    ggplot2::ggplot(
        plot_data,
        ggplot2::aes(
            x    = .data$interaction_type,
            y    = .data$relation_count,
            fill = .data$fill_hex
        )
    ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::scale_fill_identity(
            guide  = "legend",
            name   = "Source",
            breaks = unname(category_colour("resources", src_order)),
            labels = src_order
        ) +
        ggplot2::labs(
            title = "Azelate: interaction types by source",
            x     = NULL,
            y     = "Relation count"
        ) +
        theme_bw_metabo(
            width_mm     = width_mm,
            font_scale   = font_scale,
            legend_scale = legend_scale
        ) +
        ggplot2::theme(
            axis.text.x     = ggplot2::element_text(
                angle = 45, hjust = 1
            ),
            legend.position = "right"
        )
}


#' Figure 6 Panel D -- Azelate cancer associations by sample type
#'
#' Stacked-bar plot of evidence counts per \code{disease_type},
#' fill-coloured by \code{sample_type}. Reproduces the
#' \code{make_fig3_panel_D_cancer_stacked} reference from the
#' azelate Python notebook through the shared style module.
#'
#' Bars are ordered by descending total evidence count; the
#' sample-type stack order is by descending total contribution
#' across the panel.
#'
#' @param data Tibble from \code{\link{azelate_panel_d_data}}.
#' @param width_mm Numeric: target physical panel width in mm.
#' @param font_scale Numeric: passed to \code{\link{theme_bw_metabo}}.
#' @param legend_scale Numeric or \code{NULL}: passed to
#'     \code{\link{theme_bw_metabo}}.
#'
#' @return A ggplot object.
#'
#' @importFrom dplyr group_by summarise arrange desc mutate
#' @importFrom stringr str_wrap
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_identity
#' @importFrom ggplot2 scale_x_discrete labs theme element_text
#' @importFrom rlang .data
#' @export
azelate_disease_panel <- function(
    data,
    width_mm = 89L,
    font_scale = 1,
    legend_scale = NULL
) {

    # NSE vs. R CMD check workaround
    disease_type <- sample_type <- evidence_count <- total <- NULL
    fill_hex <- samp_total <- NULL

    bar_order <- dplyr::arrange(
        dplyr::summarise(
            dplyr::group_by(data, .data$disease_type),
            total = sum(.data$evidence_count),
            .groups = "drop"
        ),
        dplyr::desc(total)
    )$disease_type

    samp_order <- dplyr::arrange(
        dplyr::summarise(
            dplyr::group_by(data, .data$sample_type),
            samp_total = sum(.data$evidence_count),
            .groups = "drop"
        ),
        dplyr::desc(samp_total)
    )$sample_type

    plot_data <- dplyr::mutate(
        data,
        disease_type = factor(.data$disease_type, levels = bar_order),
        sample_type  = factor(.data$sample_type, levels = samp_order),
        fill_hex     = category_colour(
            "sample_types", as.character(.data$sample_type)
        )
    )

    ggplot2::ggplot(
        plot_data,
        ggplot2::aes(
            x    = .data$disease_type,
            y    = .data$evidence_count,
            fill = .data$fill_hex
        )
    ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::scale_x_discrete(
            labels = function(x) stringr::str_wrap(x, width = 14L)
        ) +
        ggplot2::scale_fill_identity(
            guide  = "legend",
            name   = "Sample type",
            breaks = unname(category_colour("sample_types", samp_order)),
            labels = samp_order
        ) +
        ggplot2::labs(
            title = "Azelate: cancer associations by sample type",
            x     = NULL,
            y     = "Evidence count"
        ) +
        theme_bw_metabo(
            width_mm     = width_mm,
            font_scale   = font_scale,
            legend_scale = legend_scale
        ) +
        ggplot2::theme(
            axis.text.x     = ggplot2::element_text(
                angle = 45, hjust = 1
            ),
            legend.position = "right"
        )
}
