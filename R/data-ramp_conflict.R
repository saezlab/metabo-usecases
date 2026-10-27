#' RaMP InChIKey-conflict counts per conflict reason (FR-007f)
#'
#' Reads the cycle-001 \code{metabo_ramp_inchikey_conflict} table on
#' \code{dev4} (the structural-specificity / RaMP-conflict deployment
#' per FR-030). The table classifies each RaMP id that maps to
#' multiple distinct InChIKeys by the reason for the discrepancy:
#' \code{stereo} (stereochemistry only), \code{specificity}
#' (different specificity levels), \code{tautomer} (tautomeric
#' variants), \code{similar} (similar but not equivalent structures),
#' \code{unrelated} (no chemical similarity). The classification is
#' produced in-database via the RDKit Postgres cartridge — no
#' in-pipeline RDKit invocation here.
#'
#' \code{pg_query_panel} routes the query to \code{dev4}
#' automatically via the FR-030 registry entry
#' (\code{database-content} + \code{facet = "ramp_conflict"}).
#'
#' @param panel_id Character: panel identifier. Default
#'     \code{"database-content"}.
#'
#' @return A tibble with \code{conflict_reason} and \code{n},
#'     ordered by \code{n} descending.
#'
#' @examples
#' \dontrun{
#' counts <- ramp_conflict_counts()
#' attr(counts, "deployment")   # "dev4"
#' }
#'
#' @importFrom DBI dbGetQuery
#' @export
ramp_conflict_counts <- function(panel_id = "database-content") {

    sql <- "
        SELECT conflict_reason,
               COUNT(*)::bigint AS n
        FROM   metabo_ramp_inchikey_conflict
        GROUP  BY conflict_reason
        ORDER  BY n DESC, conflict_reason
    "
    pg_query_panel(panel_id, sql, facet = "ramp_conflict")
}
