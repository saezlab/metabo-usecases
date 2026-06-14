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
#' @param node_sizes Tibble from \code{\link{fr007c_node_sizes}} with
#'     per-resource entity / interaction counts used to scale node
#'     size. Optional — when \code{NULL} the function queries it.
#'
#' @importFrom igraph graph_from_data_frame V degree
#' @importFrom ggraph ggraph geom_edge_link geom_node_point
#' @importFrom ggraph geom_node_text scale_edge_width
#' @importFrom ggplot2 labs theme element_text element_blank aes
#' @importFrom ggplot2 scale_size guides guide_legend
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom rlang .data
#' @export
plot_fr007c_networks <- function(data,
                                 min_overlap = 100L,
                                 width_mm    = 320L,
                                 node_sizes  = NULL) {

    kind_labels <- c(entity = "Molecular entities",
                     relation = "Interactions")

    data <- data[data$overlap >= as.integer(min_overlap), , drop = FALSE]

    if (is.null(node_sizes)) {
        node_sizes <- fr007c_node_sizes("fig02-overview")
    }

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

    # qgraph Fruchterman-Reingold with the area / repulse parameters
    # from the user's snippet — nodes distribute more evenly inside
    # the bounding disc than layout_with_kk's tight-core / sparse-rim
    # default.
    set.seed(2026L)
    el <- igraph::as_edgelist(union_g, names = FALSE)
    vc <- igraph::vcount(union_g)
    union_layout <- qgraph::qgraph.layout.fruchtermanreingold(
        el,
        vcount      = vc,
        area        = vc ^ 2.3,
        repulse.rad = vc ^ 2.1,
        niter       = 3000L
    )
    layout_df <- data.frame(
        name = igraph::V(union_g)$name,
        x    = union_layout[, 1L],
        y    = union_layout[, 2L]
    )

    # Per-kind plot.
    make_kind_plot <- function(kind, show_singletons = TRUE) {
        kdata <- data[data$content_kind == kind, , drop = FALSE]
        if (nrow(kdata) == 0L) return(NULL)

        # Session 2026-06-14: the entity panel is too dense to show a
        # pattern, so drop the thinnest edges. Keep edges in the top
        # 30 % by overlap weight. The interactions / literature panels
        # have fewer edges already and keep their full set.
        if (identical(kind, "entity") && nrow(kdata) > 0L) {
            cutoff <- stats::quantile(
                kdata$overlap, probs = 0.70, na.rm = TRUE
            )
            kdata <- kdata[kdata$overlap >= cutoff, , drop = FALSE]
        }

        edges <- data.frame(
            from   = kdata$label_a,
            to     = kdata$label_b,
            weight = kdata$overlap
        )

        # Node-size column for this content kind, joined onto the
        # shared layout. Resources without size data (rare) get the
        # min size so they still show up.
        kind_sizes <- node_sizes[node_sizes$content_kind == kind, ,
                                 drop = FALSE]
        nodes <- merge(
            layout_df,
            data.frame(name = kind_sizes$label, size_n = kind_sizes$n),
            by   = "name",
            all.x = TRUE,
            sort  = FALSE
        )
        # Re-order to layout_df order so the layout (x, y) lines up
        # with the node table for ggraph.
        nodes <- nodes[match(layout_df$name, nodes$name), , drop = FALSE]
        # Clamp NA / zero sizes so log scale stays defined.
        nodes$size_n <- pmax(
            ifelse(is.na(nodes$size_n), 1, as.numeric(nodes$size_n)),
            1
        )

        g <- igraph::graph_from_data_frame(
            edges, vertices = nodes, directed = FALSE
        )

        if (!show_singletons) {
            # Drop nodes whose degree is 0 in this kind's edge list.
            # The layout coordinates stay anchored to the union graph
            # so cross-panel comparability is preserved for the
            # remaining nodes.
            deg <- igraph::degree(g)
            keep_nodes <- names(deg)[deg > 0L]
            g <- igraph::induced_subgraph(g, vids = keep_nodes)
            layout_subset <- layout_df[
                match(keep_nodes, layout_df$name), , drop = FALSE
            ]
        } else {
            layout_subset <- layout_df[
                match(igraph::V(g)$name, layout_df$name), ,
                drop = FALSE
            ]
        }

        # Node colour — bright green pairs well with the dark-blue
        # edge colour (the prior palette's first entry).
        node_colour <- "#3FA34D"
        edge_colour <- palette_lead()[[1L]]

        ggraph::ggraph(g, layout = "manual",
                       x = layout_subset$x, y = layout_subset$y) +
            ggraph::geom_edge_link(
                ggplot2::aes(width = .data$weight),
                alpha  = 0.45,
                colour = edge_colour
            ) +
            ggraph::geom_node_point(
                ggplot2::aes(size = .data$size_n),
                colour = node_colour,
                alpha  = 0.95
            ) +
            ggraph::geom_node_text(
                ggplot2::aes(label = .data$name),
                size = 2.2, repel = TRUE,
                max.overlaps = 30L,
                box.padding   = 0.1,
                point.padding = 0.1
            ) +
            ggraph::scale_edge_width(
                range = c(0.2, 2.2),
                trans = "log10",
                name  = "Resource\noverlap"
            ) +
            ggplot2::scale_size(
                range = c(0.8, 8),
                trans = "log10",
                labels = scales::label_number(
                    scale_cut = scales::cut_short_scale()
                ),
                name   = if (identical(kind, "entity")) "Entities"
                         else "Interactions"
            ) +
            ggplot2::guides(
                size = ggplot2::guide_legend(
                    order = 1L,
                    keywidth  = grid::unit(2, "mm"),
                    keyheight = grid::unit(2, "mm")
                ),
                edge_width = ggplot2::guide_legend(
                    order = 2L,
                    keywidth  = grid::unit(4, "mm"),
                    keyheight = grid::unit(1.5, "mm")
                )
            ) +
            ggplot2::labs(title = kind_labels[[kind]]) +
            ggplot2::theme(
                panel.background = ggplot2::element_blank(),
                plot.background  = ggplot2::element_blank(),
                axis.text        = ggplot2::element_blank(),
                axis.title       = ggplot2::element_blank(),
                axis.ticks       = ggplot2::element_blank(),
                plot.margin      = ggplot2::margin(2, 2, 2, 2),
                plot.title       = ggplot2::element_text(
                    face = "bold", size = 9, hjust = 0.5,
                    margin = ggplot2::margin(b = 1)
                ),
                legend.position  = "right",
                legend.justification = c(0, 0.5),
                legend.box.margin    = ggplot2::margin(l = 1),
                legend.text      = ggplot2::element_text(size = 6),
                legend.title     = ggplot2::element_text(
                    size = 7, face = "bold"
                ),
                legend.key.size  = grid::unit(2, "mm"),
                legend.spacing.y = grid::unit(0.3, "mm")
            )
    }

    plots <- list(
        make_kind_plot("entity",   show_singletons = TRUE),
        # Session 2026-06-14 (R7 review): both networks share the
        # union-graph node set so the layout is *identical*. Nodes
        # without interaction edges show up as isolated dots at the
        # same screen position they occupy in the Molecular-entities
        # panel — the reader can visually pair every resource
        # between the two networks without doing a node-by-node
        # cross-reference.
        make_kind_plot("relation", show_singletons = TRUE)
    )
    plots <- Filter(Negate(is.null), plots)

    # No inner plot_annotation here — when this patchwork is
    # composed into Figure 2 by figures/fig02-overview/build.R,
    # the outer plot_annotation(tag_levels = "A") owns the figure
    # caption + tag styling. The methodological detail (shared
    # qgraph FR layout, min_overlap = 100) lives in caption.tex
    # instead.
    patchwork::wrap_plots(plots, nrow = 1L)
}
