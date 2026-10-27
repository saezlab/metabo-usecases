#' Number of distinct molecular entities per upstream resource
#'
#' Figure 1 Panel B. Reads the per-resource entity_count from the
#' cycle-001 \code{resources} table (pre-computed by the derive
#' phase), so this is a single sequential scan instead of a
#' \code{COUNT(DISTINCT entity_id)} join across tens of millions of
#' rows.
#'
#' @param panel_id Character: panel identifier (default
#'     \code{"database-content"}). Passed to
#'     \code{\link{pg_query_panel}} so the deployment routing is
#'     recorded in the result's \code{"deployment"} attribute
#'     (which surfaces in the sidecar).
#'
#' @return A tibble with \code{resource} and \code{n_entities},
#'     ordered by \code{n_entities} descending.
#'
#' @importFrom DBI dbGetQuery
#' @export
entities_by_resource <- function(panel_id = "database-content") {

    sql <- "
        SELECT resource_id AS resource,
               entity_count AS n_entities
        FROM   resources
        WHERE  entity_count > 0
        ORDER  BY entity_count DESC, resource_id
    "
    pg_query_panel(panel_id, sql)
}


#' Number of interactions per upstream resource
#'
#' Figure 1 Panel C. Reads the pre-computed
#' \code{resources.interaction_count}.
#'
#' @inheritParams entities_by_resource
#'
#' @return A tibble with \code{resource} and \code{n_relations},
#'     ordered by \code{n_relations} descending.
#'
#' @importFrom DBI dbGetQuery
#' @export
interactions_by_resource <- function(panel_id = "database-content") {

    sql <- "
        SELECT resource_id AS resource,
               interaction_count AS n_relations
        FROM   resources
        WHERE  interaction_count > 0
        ORDER  BY interaction_count DESC, resource_id
    "
    pg_query_panel(panel_id, sql)
}


#' Number of interactions per interaction-class vocabulary
#'
#' Figure 1 Panel D. Joins \code{vocab_relation_predicate} to
#' \code{vocab_interaction_class} so the rows surface the coarse
#' interaction class (\code{Signaling} / \code{Transport} /
#' \code{Other}) per the 2026-06-09 handover: the cycle-001
#' predicate vocabulary populates only those three classes; finer
#' categories (TF-target, allosteric, ligand-receptor, drug-target)
#' collapse into \code{Other}. The query also returns the per-class
#' \emph{predicate} list so a renderer can surface predicate-level
#' detail without a second query.
#'
#' @inheritParams entities_by_resource
#'
#' @return A tibble with \code{interaction_class}, \code{predicate}
#'     (text array), and \code{n} per class, ordered by \code{n}
#'     descending.
#'
#' @importFrom DBI dbGetQuery
#' @export
interactions_by_type <- function(panel_id = "database-content") {

    sql <- "
        SELECT
            COALESCE(vic.name, 'Other') AS interaction_class,
            ARRAY_AGG(DISTINCT vrp.name ORDER BY vrp.name) AS predicates,
            SUM(frb.relation_count)::bigint AS n
        FROM   facet_relation_bitmap frb
        JOIN   vocab_relation_predicate vrp
               ON vrp.name = frb.facet_value
        LEFT   JOIN vocab_interaction_class vic
               ON vrp.interaction_class_id = vic.interaction_class_id
        WHERE  frb.facet_name = 'predicate'
        GROUP  BY COALESCE(vic.name, 'Other')
        ORDER  BY n DESC, interaction_class
    "
    pg_query_panel(panel_id, sql)
}


#' Number of associations per resource
#'
#' Figure 1 Panel E (FR-007a Associations facet). Reads the
#' pre-computed \code{resources.association_count}. Replaces the
#' pre-handover \code{annotation_classes_by_resource()} which joined
#' the unpopulated \code{entity_annotation_relation} table and so
#' returned zero rows on cycle-001 builds.
#'
#' @inheritParams entities_by_resource
#'
#' @return A tibble with \code{resource} and \code{n_associations},
#'     ordered by \code{n_associations} descending.
#'
#' @importFrom DBI dbGetQuery
#' @export
associations_by_resource <- function(panel_id = "database-content") {

    sql <- "
        SELECT resource_id AS resource,
               association_count AS n_associations
        FROM   resources
        WHERE  association_count > 0
        ORDER  BY association_count DESC, resource_id
    "
    pg_query_panel(panel_id, sql)
}


#' Number of ontology terms per ontology prefix
#'
#' Figure 1 Panel F. Replaces the deprecated \code{ontology_terms}
#' table read with a scan of the cycle-001 entity-relations
#' \code{entity_ontology_term} table. The top-N ontologies are
#' returned as-is; everything below the cap collapses into a single
#' \code{Misc} bucket so the panel honours the FR-009c distinct-colour
#' cap (10–12 categories).
#'
#' @inheritParams entities_by_resource
#' @param top_n Integer: how many ontologies to keep before the
#'     \code{Misc} collapse. Default \code{10L} (≤ FR-009c cap).
#'
#' @return A tibble with \code{ontology} and \code{n_terms},
#'     ordered by \code{n_terms} descending. At most
#'     \code{top_n + 1} rows (the trailing \code{Misc} appears only
#'     when the long tail exists).
#'
#' @importFrom DBI dbGetQuery
#' @export
ontology_terms_by_ontology <- function(panel_id = "database-content",
                                       top_n = 10L) {

    sql <- sprintf("
        WITH ranked AS (
            SELECT
                COALESCE(ontology_prefix, 'unspecified') AS ontology,
                COUNT(*) AS n_terms,
                ROW_NUMBER() OVER (ORDER BY COUNT(*) DESC) AS rk
            FROM   entity_ontology_term
            GROUP  BY ontology_prefix
        )
        SELECT
            CASE WHEN rk <= %d THEN ontology ELSE 'Misc' END AS ontology,
            SUM(n_terms)::bigint AS n_terms
        FROM   ranked
        GROUP  BY CASE WHEN rk <= %d THEN ontology ELSE 'Misc' END
        ORDER  BY n_terms DESC, ontology
    ", as.integer(top_n), as.integer(top_n))
    pg_query_panel(panel_id, sql)
}
