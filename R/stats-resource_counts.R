#' SQL templates per FR-043f section
#'
#' Each section returns distinct \code{data_source.name} values
#' contributing >= 1 row to the section's underlying join, plus the
#' row count per source. Kept as a single dispatch table so the
#' orchestrator can hash the verbatim SQL into the digest sidecar.
#'
#' Post-2026-06-14 dev5 integrated-build promotion every section runs
#' against the single \code{dev5} deployment; the per-panel split
#' established at the 2026-06-09 handover is retired.
#'
#' Schema assumptions match the integrated build's tables:
#' \itemize{
#'   \item \code{entity_evidence(entity_id, source_id)}
#'   \item \code{relation(relation_id, subject_id, object_id,
#'         predicate_id, relation_category_id)} +
#'         \code{relation_evidence_relation(relation_id, source_id)}
#'   \item \code{vocab_relation_category(name)}
#'   \item \code{metabo_entity_structural_specificity(entity_id, …)}
#'   \item \code{entity_evidence_annotation(entity_id, source_id, …)}
#'   \item \code{data_source(source_id, name)}
#' }
#'
#' @return Named list keyed by section id (\code{1L..5L}).
#'
#' @keywords internal
#' @noRd
section_resource_sql <- function() {

    list(
        # Section 1 — Entities
        `1` = paste(
            "SELECT ds.name AS resource_name,",
            "       COUNT(DISTINCT ee.entity_id) AS n_rows",
            "  FROM entity_evidence ee",
            "  JOIN data_source ds USING (source_id)",
            " GROUP BY ds.name",
            " ORDER BY n_rows DESC",
            sep = "\n"
        ),
        # Section 2 — Metabolite-Protein Interactions
        `2` = paste(
            "WITH chem AS (",
            "    SELECT entity_id FROM entity e",
            "      JOIN vocab_entity_type vet USING (entity_type_id)",
            "     WHERE vet.name = 'Chemical:OM:0037'",
            "), pog AS (",
            "    SELECT entity_id FROM entity e",
            "      JOIN vocab_entity_type vet USING (entity_type_id)",
            "     WHERE vet.name IN ('Gene:MI:0250','Protein:MI:0326')",
            "), mpi AS (",
            "    SELECT r.relation_id",
            "      FROM relation r",
            "      JOIN vocab_relation_category vrc",
            "        ON vrc.relation_category_id = r.relation_category_id",
            "     WHERE vrc.name = 'interaction'",
            "       AND ( (r.subject_id IN (SELECT entity_id FROM chem)",
            "              AND r.object_id  IN (SELECT entity_id FROM pog))",
            "          OR (r.object_id  IN (SELECT entity_id FROM chem)",
            "              AND r.subject_id IN (SELECT entity_id FROM pog)) )",
            ")",
            "SELECT ds.name AS resource_name,",
            "       COUNT(DISTINCT mpi.relation_id) AS n_rows",
            "  FROM mpi",
            "  JOIN relation_evidence_relation rer USING (relation_id)",
            "  JOIN data_source ds USING (source_id)",
            " GROUP BY ds.name",
            " ORDER BY n_rows DESC",
            sep = "\n"
        ),
        # Section 3 — Interactions (includes MPI per FR-043c)
        `3` = paste(
            "SELECT ds.name AS resource_name,",
            "       COUNT(DISTINCT r.relation_id) AS n_rows",
            "  FROM relation r",
            "  JOIN vocab_relation_category vrc",
            "    ON vrc.relation_category_id = r.relation_category_id",
            "  JOIN relation_evidence_relation rer USING (relation_id)",
            "  JOIN data_source ds USING (source_id)",
            " WHERE vrc.name = 'interaction'",
            " GROUP BY ds.name",
            " ORDER BY n_rows DESC",
            sep = "\n"
        ),
        # Section 4 — Structures
        `4` = paste(
            "SELECT ds.name AS resource_name,",
            "       COUNT(DISTINCT mess.entity_id) AS n_rows",
            "  FROM metabo_entity_structural_specificity mess",
            "  JOIN entity_evidence ee USING (entity_id)",
            "  JOIN data_source ds USING (source_id)",
            " GROUP BY ds.name",
            " ORDER BY n_rows DESC",
            sep = "\n"
        ),
        # Section 5 — Annotation
        `5` = paste(
            "SELECT ds.name AS resource_name,",
            "       COUNT(*) AS n_rows",
            "  FROM entity_evidence_annotation eea",
            "  JOIN data_source ds USING (source_id)",
            " GROUP BY ds.name",
            " ORDER BY n_rows DESC",
            sep = "\n"
        )
    )
}


#' Per-section contributing-resource list (FR-043f)
#'
#' Returns the distinct \code{data_source.name} values contributing
#' >= 1 row to the named section's underlying join, plus the row
#' count per source. The full list is written into
#' \code{stats.json -> sections[i].resources} (not just the count) so
#' the author can cross-check which resources back each headline
#' number on the Panel A diagram.
#'
#' Runs via \code{\link{pg_query_panel}} so the deployment routing
#' goes through the FR-030 matrix.
#'
#' @param section_id Integer (1..5).
#' @param panel_id Character: panel id used for deployment routing.
#'     Defaults to \code{"fig01-overview"}.
#'
#' @return A tibble with columns \code{resource_name} (character) and
#'     \code{n_rows} (integer). The tibble also carries
#'     \code{attr(., "sql")} and \code{attr(., "result_hash")} +
#'     \code{attr(., "deployment")} from \code{pg_query_panel}.
#'
#' @importFrom logger log_info
#' @importFrom rlang abort
#' @export
section_contributing_resources <- function(
    section_id,
    panel_id = "fig01-overview"
) {

    key <- as.character(as.integer(section_id))

    sql_map <- section_resource_sql()
    sql <- sql_map[[key]]
    if (is.null(sql)) {
        rlang::abort(sprintf(
            "section_contributing_resources: section_id=%s is not 1..5",
            key
        ))
    }

    facet <- paste0("panel_a_stats_section_", key, "_resources")

    rows <- pg_query_panel(panel_id, sql, facet = facet)

    logger::log_info(
        "section_contributing_resources(section={key}) → ",
        "{nrow(rows)} contributing data sources"
    )

    rows
}
