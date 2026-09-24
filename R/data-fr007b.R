#' FR-007b coverage profile data layer
#'
#' For each variant, returns the cumulative coverage curve: for each
#' \code{n_resources = N}, the number of items supported by AT LEAST
#' \code{N} resources. Implemented as a single CTE per variant
#' aggregating \code{entity_source_count} (cycle-001 pre-computed
#' per-entity source-count table) over the appropriate entity filter.
#'
#' Variants implemented at the data layer:
#' \describe{
#'   \item{\code{entities}}{All entities in \code{entity_source_count}.
#'     Whole-DB view.}
#'   \item{\code{molecular_entities}}{Filtered to the FR-007a
#'     "molecular" entity types — Chemical, Protein, Gene, Mirna,
#'     Complex, Physical Entity.}
#'   \item{\code{structures}}{Filtered to entities that have a
#'     non-\code{no_structure} bucket in
#'     \code{metabo_entity_structural_specificity}; routes to
#'     \code{dev4}.}
#' }
#'
#' Not yet implemented (would benefit from derived
#' \code{identifier_source_count} / \code{relation_source_count}
#' tables; see
#' \code{saezverse/human/plans/omnipath-improvements-2026-06-identifier-source-count.md}):
#' \code{identifiers}, \code{literature}, \code{interactions}.
#'
#' @param variant Character: one of \code{"entities"},
#'     \code{"molecular_entities"}, \code{"structures"}, or
#'     \code{"all"} to return a combined long-format tibble of every
#'     supported variant.
#' @param panel_id Character: panel identifier (default
#'     \code{"database-content"}).
#'
#' @return A tibble with \code{n_resources}, \code{n_items_ge_N}, and
#'     (for \code{variant = "all"}) a \code{variant} column.
#'
#' @examples
#' \dontrun{
#' fr007b_coverage("entities")
#' fr007b_coverage("all")
#' }
#'
#' @importFrom DBI dbGetQuery
#' @importFrom dplyr bind_rows
#' @export
fr007b_coverage <- function(
    variant = c("entities", "molecular_entities", "structures", "all"),
    panel_id = "database-content"
) {
    variant <- match.arg(variant)

    if (variant == "all") {
        out <- list()
        for (v in c("entities", "molecular_entities", "structures")) {
            rows <- fr007b_coverage(v, panel_id)
            rows$variant <- v
            out[[v]] <- rows
        }
        return(dplyr::bind_rows(out))
    }

    # Each variant: an entity-id source query that selects
    # entity_ids in entity_source_count matching the variant's
    # filter, plus the source_count window.
    src <- switch(
        variant,
        # Whole-DB view. The cycle-001 entity_source_count table only
        # populates resolved Chemical / Gene / Complex / Pathway /
        # Reaction rows (the long tail of Cv Terms, Mirna, Protein,
        # Physical Entity, Tissue, Phenotype etc. has no row in
        # entity_source_count). To honour the "all entities" semantic
        # the curve includes those tail entities at source_count = 1
        # — they exist in the DB so their cumulative-coverage curve
        # contribution starts at the 1-source bucket.
        entities = "
            WITH all_entity_source AS (
                SELECT entity_id, source_count FROM entity_source_count
                UNION ALL
                SELECT e.entity_id, 1 AS source_count
                FROM   entity e
                LEFT   JOIN entity_source_count esc USING (entity_id)
                WHERE  esc.entity_id IS NULL
            )
            SELECT source_count, COUNT(*)::bigint AS n
            FROM   all_entity_source
            GROUP  BY source_count
        ",
        molecular_entities = "
            SELECT esc.source_count, COUNT(*)::bigint AS n
            FROM   entity_source_count esc
            JOIN   entity e USING (entity_id)
            JOIN   vocab_entity_type vet
                   ON vet.entity_type_id = e.entity_type_id
            WHERE  vet.name IN (
                'Chemical:OM:0037',
                'Protein:MI:0326',
                'Gene:MI:0250',
                'Mirna:OM:0038',
                'Complex:MI:0314',
                'Physical Entity:OM:0016'
            )
            GROUP  BY esc.source_count
        ",
        structures = NULL  # special: routes to dev4
    )

    if (variant == "structures") {
        # dev4-only — metabo_entity_structural_specificity exists there.
        sql <- "
            WITH src AS (
                SELECT esc.source_count, COUNT(*)::bigint AS n
                FROM   entity_source_count esc
                JOIN   metabo_entity_structural_specificity mess
                       USING (entity_id)
                WHERE  mess.structural_specificity_id != (
                    SELECT structural_specificity_id
                    FROM   metabo_vocab_structural_specificity
                    WHERE  name = 'no_structure'
                )
                GROUP  BY esc.source_count
            )
            SELECT source_count AS n_resources,
                   SUM(n) OVER (ORDER BY source_count DESC)::bigint
                       AS n_items_ge_n
            FROM   src
            ORDER  BY source_count
        "
        return(pg_query_panel(panel_id, sql, facet = "panel_e"))
    }

    sql <- sprintf("
        WITH src AS (%s)
        SELECT source_count AS n_resources,
               SUM(n) OVER (ORDER BY source_count DESC)::bigint
                   AS n_items_ge_n
        FROM   src
        ORDER  BY source_count
    ", src)
    pg_query_panel(panel_id, sql)
}
