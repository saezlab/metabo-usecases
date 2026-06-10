#' Render the FR-007c resource-overlap networks panel
#'
#' Composes one network plot per available content kind
#' (\code{Molecular entities}, \code{Interactions}) into a single
#' patchwork. Per FR-007c the LAYOUT is computed once on the union
#' graph and reused for every sub-network — the same resource sits
#' at the same screen position in every panel so the eye can
#' compare overlaps directly.
#'
#' Edge thickness scales with the log of the overlap count (counts
#' span 4-5 orders of magnitude). Edges below the threshold
#' \code{min_overlap} (default 100) are pruned to leave only the
#' strongest overlaps visible.
#'
#' @param data Tibble from \code{\link{fr007c_overlap}}.
#' @param min_overlap Integer: threshold below which edges are
#'     pruned. Default 100.
#' @param width_mm Numeric.
#'
#' @return A patchwork composite.
#'
#' @importFrom igraph graph_from_data_frame layout_with_kk V degree
#' @importFrom ggraph ggraph geom_edge_link geom_node_point
#' @importFrom ggraph geom_node_text scale_edge_width
#' @importFrom ggplot2 labs theme element_text element_blank
#' @importFrom ggplot2 aes
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom rlang .data
#' @export
plot_fr007c_networks <- function(data,
                                 min_overlap = 100L,
                                 width_mm    = 320L) {

    kind_labels <- c(entity = "Molecular entities",
                     relation = "Interactions")

    data <- data[data$overlap >= as.integer(min_overlap), , drop = FALSE]

    # Union graph for the shared layout.
    union_edges <- unique(rbind(
        data.frame(from = data$label_a, to = data$label_b),
        data.frame(from = data$label_b, to = data$label_a)
    ))
    union_edges <- union_edges[union_edges$from < union_edges$to, ,
                               drop = FALSE]
    union_g <- igraph::graph_from_data_frame(
        union_edges, directed = FALSE
    )
    set.seed(2026L)
    union_layout <- igraph::layout_with_kk(union_g)
    layout_df <- data.frame(
        name = igraph::V(union_g)$name,
        x    = union_layout[, 1L],
        y    = union_layout[, 2L]
    )

    # Per-kind plot.
    make_kind_plot <- function(kind) {
        kdata <- data[data$content_kind == kind, , drop = FALSE]
        if (nrow(kdata) == 0L) return(NULL)

        edges <- data.frame(
            from   = kdata$label_a,
            to     = kdata$label_b,
            weight = kdata$overlap
        )
        nodes <- layout_df

        g <- igraph::graph_from_data_frame(
            edges, vertices = nodes, directed = FALSE
        )

        ggraph::ggraph(g, layout = "manual",
                       x = layout_df$x, y = layout_df$y) +
            ggraph::geom_edge_link(
                ggplot2::aes(width = .data$weight),
                alpha = 0.4, colour = palette_lead()[[1L]]
            ) +
            ggraph::geom_node_point(
                size = 2.5, colour = palette_lead()[[1L]]
            ) +
            ggraph::geom_node_text(
                ggplot2::aes(label = .data$name),
                size = 3.5, repel = TRUE,
                max.overlaps = 30L
            ) +
            ggraph::scale_edge_width(
                range = c(0.2, 2.5),
                trans = "log10",
                guide = "none"
            ) +
            ggplot2::labs(title = kind_labels[[kind]]) +
            ggplot2::theme(
                panel.background = ggplot2::element_blank(),
                plot.background  = ggplot2::element_blank(),
                axis.text        = ggplot2::element_blank(),
                axis.title       = ggplot2::element_blank(),
                axis.ticks       = ggplot2::element_blank(),
                plot.title       = ggplot2::element_text(
                    face = "bold", size = 15, hjust = 0.5
                )
            )
    }

    plots <- list(
        make_kind_plot("entity"),
        make_kind_plot("relation")
    )
    plots <- Filter(Negate(is.null), plots)

    patchwork::wrap_plots(plots, nrow = 1L) +
        patchwork::plot_annotation(
            caption = sprintf(
                "Shared kamada-kawai layout; edges with overlap >= %d shown.",
                as.integer(min_overlap)
            )
        )
}
