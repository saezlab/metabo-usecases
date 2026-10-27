#' Render the FR-007d entity × interaction-type matrix panel
#'
#' Bar plot with vertical facets by participant entity type; bars
#' show interaction class (Signaling / Transport / Other) counts.
#' Spec asked for Proteins / Metabolites / Lipids / Drugs; the
#' \code{participant_type} facet has the broader entity-type set
#' from the cycle-001 build (Chemical, Gene, Food, Cv Term, Protein,
#' Reaction, Physical Entity, Pathway, Complex, Organism, Mirna, …).
#' We use the top-N participant types by relation count.
#'
#' Log10 y because relation counts span 5+ orders of magnitude.
#'
#' @param data Tibble from
#'     \code{\link{entity_by_interaction_type}}.
#' @param width_mm Numeric.
#'
#' @return A ggplot.
#'
#' @importFrom ggplot2 ggplot aes geom_col facet_wrap labs theme
#' @importFrom ggplot2 element_text scale_y_log10
#' @importFrom rlang .data
#' @export
plot_entity_interaction_matrix <- function(data, width_mm = 180L) {

    # Order facets by total relation count (preserved by the SQL).
    type_order <- unique(as.character(data$participant_type))
    data$participant_type <- factor(
        data$participant_type, levels = type_order
    )
    # Strip the trailing :MI:/:OM: codes and apply per-token
    # corrections — "Cv Term" should read "CV term" in display.
    pretty_label <- function(x) {
        x <- sub("\\s*:[A-Z]+:\\d+\\s*$", "", x)
        x <- sub("^Cv Term$", "CV term", x)
        x
    }
    short_labels <- stats::setNames(pretty_label(type_order), type_order)

    # Interaction-class order: most concrete → most generic.
    class_order <- c(
        "Signaling", "Transport", "Interaction", "Reaction",
        "Association", "Membership", "Other"
    )
    data$interaction_class <- factor(
        data$interaction_class,
        levels = intersect(class_order, unique(data$interaction_class))
    )

    # Replace n=0 with NA so log10 doesn't choke; geom_col skips NA.
    data$n <- ifelse(data$n == 0, NA_real_, as.numeric(data$n))

    bar_fill <- palette_lead()[[1L]]

    ggplot2::ggplot(
        data,
        ggplot2::aes(x = .data$interaction_class, y = .data$n)
    ) +
        ggplot2::geom_col(fill = bar_fill) +
        # 2 rows × 4 cols so the 8 participant-entity facets sit
        # in a compact rectangle (Session 2026-06-14 layout).
        ggplot2::facet_wrap(
            ~ .data$participant_type, nrow = 2L,
            labeller = ggplot2::as_labeller(short_labels),
            scales = "free_y"
        ) +
        ggplot2::scale_y_log10(
            labels = scales::label_number(
                scale_cut = scales::cut_short_scale()
            )
        ) +
        ggplot2::labs(x = "Interaction type", y = "Interactions (log)") +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(
            legend.position = "none",
            plot.margin     = ggplot2::margin(2, 2, 2, 2),
            # Composite slot is ~120 mm wide × 50 mm tall (2/3 of
            # the bottom row). Sizes match Panel C (specificity) so the
            # two panels read as a coordinated barplot row.
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
