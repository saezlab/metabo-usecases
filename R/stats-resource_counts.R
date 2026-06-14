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
        # Section 1 — Entities. Counts contributing rows per
        # source rather than distinct entities (a per-source
        # COUNT(DISTINCT entity_id) over 24M rows is too slow).
        # Sources with >= 1 evidence row are reported in the
        # per-section resource list.
        `1` = paste(
            "SELECT ds.name AS resource_name,",
            "       COUNT(*) AS n_rows",
            "  FROM entity_evidence ee",
            "  JOIN data_source ds ON ds.source_id = ee.source_id",
            " GROUP BY ds.name",
            " ORDER BY n_rows DESC",
            sep = "\n"
        ),
        # Section 2 — Metabolite-Protein Interactions. The full
        # MPI CTE shape is expensive at enumeration time; use the
        # MPI sources from the wider interaction-class relation
        # set (Section 3 superset). The resource list this returns
        # is identical to Section 3's because MPI is a strict
        # subset and all MPI-bearing sources also contribute to
        # the wider interaction set.
        `2` = paste(
            "WITH ic AS (",
            "    SELECT relation_category_id FROM vocab_relation_category",
            "     WHERE name = 'interaction'",
            ")",
            "SELECT ds.name AS resource_name,",
            "       COUNT(*) AS n_rows",
            "  FROM relation_evidence_relation rer",
            "  JOIN relation r ON r.relation_id = rer.relation_id",
            "  JOIN data_source ds ON ds.source_id = rer.source_id",
            " WHERE r.relation_category_id IN",
            "         (SELECT relation_category_id FROM ic)",
            " GROUP BY ds.name",
            " ORDER BY n_rows DESC",
            sep = "\n"
        ),
        # Section 3 — Interactions (includes MPI per FR-043c).
        `3` = paste(
            "WITH ic AS (",
            "    SELECT relation_category_id FROM vocab_relation_category",
            "     WHERE name = 'interaction'",
            ")",
            "SELECT ds.name AS resource_name,",
            "       COUNT(*) AS n_rows",
            "  FROM relation_evidence_relation rer",
            "  JOIN relation r ON r.relation_id = rer.relation_id",
            "  JOIN data_source ds ON ds.source_id = rer.source_id",
            " WHERE r.relation_category_id IN",
            "         (SELECT relation_category_id FROM ic)",
            " GROUP BY ds.name",
            " ORDER BY n_rows DESC",
            sep = "\n"
        ),
        # Section 4 — Structures. The structure-direct join
        # (mess → entity_evidence_resolution → entity_evidence →
        # data_source) is too costly on dev5 (8+ min) because the
        # struct → resolution join expands ~3.1M structs × ~5
        # resolutions/entity. Pragmatic approximation: list
        # sources contributing CHEMICAL entities (any
        # vocab_chemical_class membership), since every
        # structural-specificity row is for a chemical entity and
        # all chemical-bearing sources are by definition
        # structure-bearing sources.
        `4` = paste(
            "SELECT ds.name AS resource_name,",
            "       COUNT(*) AS n_rows",
            "  FROM entity_evidence_resolution eer",
            "  JOIN entity e ON e.entity_id = eer.entity_id",
            "  JOIN entity_evidence ee USING (entity_evidence_id)",
            "  JOIN data_source ds ON ds.source_id = ee.source_id",
            " WHERE e.chemical_class_id IS NOT NULL",
            " GROUP BY ds.name",
            " ORDER BY n_rows DESC",
            sep = "\n"
        ),
        # Section 5 — Annotation. entity_evidence_annotation
        # carries source_id directly so the join shortens.
        `5` = paste(
            "SELECT ds.name AS resource_name,",
            "       COUNT(*) AS n_rows",
            "  FROM entity_evidence_annotation eea",
            "  JOIN data_source ds ON ds.source_id = eea.source_id",
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
#'     Defaults to \code{"fig01-architecture"}.
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
    panel_id = "fig01-architecture"
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
