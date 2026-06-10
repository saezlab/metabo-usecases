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
#'     \code{\link{fr007d_entity_x_interaction}}.
#' @param width_mm Numeric.
#'
#' @return A ggplot.
#'
#' @importFrom ggplot2 ggplot aes geom_col facet_wrap labs theme
#' @importFrom ggplot2 element_text scale_y_log10
#' @importFrom rlang .data
#' @export
plot_fr007d_matrix <- function(data, width_mm = 180L) {

    # Order facets by total relation count (preserved by the SQL).
    type_order <- unique(as.character(data$participant_type))
    data$participant_type <- factor(
        data$participant_type, levels = type_order
    )
    # Shorten the type labels for the strip (strip out :MI:/:OM:
    # codes) to keep facet headings compact.
    short_labels <- stats::setNames(
        sub("\\s*:[A-Z]+:\\d+\\s*$", "", type_order),
        type_order
    )

    # Interaction-class order: Signaling, Transport, Other (Other last
    # so the longest bar / dominant class isn't always first).
    class_order <- c("Signaling", "Transport", "Other")
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
        ggplot2::labs(x = NULL, y = "Interactions (log scale)") +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(
            legend.position = "none",
            axis.text.x     = ggplot2::element_text(
                angle = 30, hjust = 1, size = 13
            ),
            axis.text.y     = ggplot2::element_text(size = 13),
            axis.title      = ggplot2::element_text(size = 15),
            strip.text      = ggplot2::element_text(face = "bold", size = 14)
        )
}
