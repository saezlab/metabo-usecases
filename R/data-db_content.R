#' Number of distinct molecular entities per upstream resource
#'
#' Figure 1 Panel B. Joins \code{data_source} to
#' \code{entity_evidence_resolution} (which links source_id → entity_id
#' once an evidence row is resolved to a canonical entity) and counts
#' distinct entities per resource.
#'
#' Source schema reference:
#' \code{omnipath-build/omnipath_build/db/schema.py} —
#' \code{data_source}, \code{entity_evidence_resolution}.
#'
#' @param con A DBI connection from \code{\link{pg_connect}}.
#' @return A tibble with \code{resource} and \code{n_entities},
#'     ordered by \code{n_entities} descending.
#'
#' @importFrom DBI dbGetQuery
#' @importFrom tibble as_tibble
#' @export
entities_by_resource <- function(con) {

    sql <- "
        SELECT ds.name AS resource,
               COUNT(DISTINCT eer.entity_id) AS n_entities
        FROM data_source ds
        JOIN entity_evidence_resolution eer USING (source_id)
        WHERE eer.entity_id IS NOT NULL
        GROUP BY ds.name
        ORDER BY n_entities DESC, ds.name
    "
    pg_query(con, sql)
}


#' Number of interactions per upstream resource
#'
#' Figure 1 Panel C. The \code{relation} table is a deduplicated graph
#' with no direct \code{source_id}; \code{relation_evidence_relation}
#' is the bridge that ties each (source, relation_evidence) to a
#' \code{relation_id}. Distinct relation_ids per source give the
#' interaction count.
#'
#' @param con A DBI connection.
#' @return A tibble with \code{resource} and \code{n_relations}.
#'
#' @importFrom DBI dbGetQuery
#' @export
interactions_by_resource <- function(con) {

    sql <- "
        SELECT ds.name AS resource,
               COUNT(DISTINCT rer.relation_id) AS n_relations
        FROM data_source ds
        JOIN relation_evidence_relation rer USING (source_id)
        GROUP BY ds.name
        ORDER BY n_relations DESC, ds.name
    "
    pg_query(con, sql)
}


#' Number of interactions per predicate (interaction type)
#'
#' Figure 1 Panel D. Uses \code{vocab_relation_predicate} for the
#' human-readable predicate name.
#'
#' @param con A DBI connection.
#' @return A tibble with \code{interaction_type} and \code{n}.
#'
#' @importFrom DBI dbGetQuery
#' @export
interactions_by_type <- function(con) {

    sql <- "
        SELECT vp.name AS interaction_type,
               COUNT(*) AS n
        FROM relation r
        JOIN vocab_relation_predicate vp
            ON r.predicate_id = vp.relation_predicate_id
        GROUP BY vp.name
        ORDER BY n DESC, vp.name
    "
    pg_query(con, sql)
}


#' Number of distinct annotation terms per resource
#'
#' Figure 1 Panel E. \code{entity_annotation_relation} ties each
#' (source, evidence, annotation) to a \code{relation_id}; the
#' \code{annotation} table stores term/value/unit per annotation_key.
#' Distinct terms-per-resource gives the "annotation classes" count.
#'
#' @param con A DBI connection.
#' @return A tibble with \code{resource} and \code{n_classes}.
#'
#' @importFrom DBI dbGetQuery
#' @export
annotation_classes_by_resource <- function(con) {

    sql <- "
        SELECT ds.name AS resource,
               COUNT(DISTINCT a.term) AS n_classes
        FROM data_source ds
        JOIN entity_annotation_relation ear USING (source_id)
        JOIN annotation a USING (annotation_key)
        GROUP BY ds.name
        ORDER BY n_classes DESC, ds.name
    "
    pg_query(con, sql)
}


#' Number of ontology terms per ontology prefix
#'
#' Figure 1 Panel F. \code{ontology_terms.ontology_prefix} is the
#' canonical short name (e.g. \code{"GO"}, \code{"CHEBI"}).
#' \code{NULL} prefixes collapse into \code{"unspecified"}.
#'
#' @param con A DBI connection.
#' @return A tibble with \code{ontology} and \code{n_terms}.
#'
#' @importFrom DBI dbGetQuery
#' @export
ontology_terms_by_ontology <- function(con) {

    sql <- "
        SELECT COALESCE(ontology_prefix, 'unspecified') AS ontology,
               COUNT(*) AS n_terms
        FROM ontology_terms
        GROUP BY ontology_prefix
        ORDER BY n_terms DESC, ontology_prefix
    "
    pg_query(con, sql)
}
