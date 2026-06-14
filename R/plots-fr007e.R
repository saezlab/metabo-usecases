#' Render the FR-007e structural specificity × chemical category panel
#'
#' Vertical bar chart faceted by chemical category, with the six
#' \code{structural_specificity} levels on the x axis and entity
#' counts on a log10 y axis (the dynamic range across cells spans
#' 6 orders of magnitude, so a linear y collapses everything onto
#' the largest bin).
#'
#' Per FR-007e the \code{no_structure} bucket is rendered as-counted
#' — the panel caption notes the dev4 build state (~44 % of
#' chemicals are still pre-T020 structure-less hashes).
#'
#' @param data Tibble from \code{\link{fr007e_specificity_by_category}}.
#' @param width_mm Numeric: target physical width.
#'
#' @return A ggplot.
#'
#' @importFrom ggplot2 ggplot aes geom_col facet_wrap labs theme
#' @importFrom ggplot2 element_text element_blank scale_y_log10
#' @importFrom ggplot2 scale_fill_manual scale_x_discrete expansion
#' @importFrom rlang .data
#' @export
plot_fr007e <- function(data, width_mm = 180L) {

    spec_levels <- c(
        "stereospecific",
        "cis_trans_only",
        "constitution_only",
        "variable_constitution",
        "unknown_constitution",
        "no_structure"
    )
    spec_labels <- c(
        "Stereospecific",
        "Cis/trans",
        "Constitution",
        "Variable",
        "Unknown",
        "No structure"
    )

    # Category ordering: chemical_class first, then metabolic_domain,
    # roughly by size (so the largest panels are at the top-left of
    # the 2-row layout).
    cat_levels <- c(
        "drugs", "metabolites", "lipids", "food compounds",
        "amino-acid metabolism", "nucleic-acid metabolism",
        "carbohydrates"
    )
    # Short facet titles — each strip cell is only ~20 mm wide
    # in the composite slot (180 mm / 3 cols / 3 facets across),
    # so the original full-length names (e.g. "Amino-acid
    # metabolism") overflow and get clipped.
    cat_labels <- c(
        "drugs"                   = "Drugs",
        "metabolites"             = "Metabolites",
        "lipids"                  = "Lipids",
        "food compounds"          = "Foods",
        "amino-acid metabolism"   = "Amino acids",
        "nucleic-acid metabolism" = "Nucleic acids",
        "carbohydrates"           = "Carbohydrates"
    )
    keep <- intersect(cat_levels, unique(as.character(data$category)))
    data$category <- factor(
        data$category,
        levels = keep,
        labels = cat_labels[keep]
    )
    data$specificity <- factor(data$specificity, levels = spec_levels)

    # Replace n = 0 with NA so log10 doesn't choke; geom_col skips NA.
    data$n <- ifelse(data$n == 0, NA_real_, as.numeric(data$n))

    # Single fill colour — the x-axis already encodes specificity, so a
    # categorical fill scale would be redundant.
    bar_fill <- palette_lead()[[1L]]

    ggplot2::ggplot(
        data,
        ggplot2::aes(
            x = .data$specificity,
            y = .data$n
        )
    ) +
        ggplot2::geom_col(fill = bar_fill) +
        ggplot2::facet_wrap(~ .data$category, nrow = 2L,
                            scales = "free_y") +
        ggplot2::scale_y_log10(
            labels = scales::label_number(
                scale_cut = scales::cut_short_scale()
            )
        ) +
        ggplot2::scale_x_discrete(labels = spec_labels) +
        ggplot2::labs(
            x = "Structural specificity",
            y = "Entities (log)"
        ) +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(
            legend.position = "none",
            plot.margin     = ggplot2::margin(2, 2, 2, 2),
            # Composite slot is ~60 mm wide × 50 mm tall; every text
            # element drops to the FR-017a 6 pt floor (or just
            # above) so the plot area dominates and the rotated
            # y-axis title stays within the panel height.
            axis.text.x     = ggplot2::element_text(
                angle = 35, hjust = 1, size = 6
            ),
            axis.text.y     = ggplot2::element_text(size = 6),
            axis.title.x    = ggplot2::element_text(
                size = 7, margin = ggplot2::margin(t = 1)
            ),
            axis.title.y    = ggplot2::element_text(
                size = 7, margin = ggplot2::margin(r = 1)
            ),
            strip.text      = ggplot2::element_text(
                face = "bold", size = 7,
                margin = ggplot2::margin(t = 0.5, b = 0.5)
            ),
            strip.background = ggplot2::element_blank(),
            panel.spacing.x  = grid::unit(1, "mm"),
            panel.spacing.y  = grid::unit(1, "mm")
        )
}
