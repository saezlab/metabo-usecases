#' FR-007c resource-overlap data layer
#'
#' Returns the pre-computed pairwise resource overlaps from
#' \code{resource_overlap_summary} (cycle-001 derived table, one row
#' per (source_a, source_b, content_kind)), with the resource ids
#' resolved to display labels via \code{\link{resources_label_map}}.
#'
#' Content kinds currently available in the cycle-001 build:
#' \code{entity} (used as "Molecular entities" in the FR-007c
#' renderer) and \code{relation} ("Interactions"). The
#' \code{literature} content kind specified by FR-007c is deferred
#' until the build emits a literature-overlap shape.
#'
#' @param panel_id Character: panel identifier.
#'
#' @return Tibble with \code{source_a}, \code{source_b},
#'     \code{content_kind}, \code{overlap}, \code{label_a},
#'     \code{label_b}.
#'
#' @importFrom DBI dbGetQuery
#' @export
fr007c_overlap <- function(panel_id = "fig01-overview") {

    sql <- "
        SELECT da.name AS source_a,
               db.name AS source_b,
               ros.content_kind,
               ros.overlap::bigint AS overlap
        FROM   resource_overlap_summary ros
        JOIN   data_source da ON da.source_id = ros.source_a_id
        JOIN   data_source db ON db.source_id = ros.source_b_id
        ORDER  BY ros.content_kind, ros.overlap DESC
    "
    res <- pg_query_panel(panel_id, sql)

    labels <- resources_label_map(panel_id)
    res$label_a <- ifelse(
        is.na(labels[res$source_a]),
        res$source_a, unname(labels[res$source_a])
    )
    res$label_b <- ifelse(
        is.na(labels[res$source_b]),
        res$source_b, unname(labels[res$source_b])
    )
    res
}
