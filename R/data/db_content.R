#' Number of molecular entities per upstream resource
#'
#' Figure 1 Panel B. One row per resource, joined via
#' \code{entity_evidence -> data_source}.
#'
#' @param con A DBI connection from \code{\link{pg_connect}}.
#' @return A tibble with columns \code{resource} and \code{n_entities},
#'     ordered by \code{n_entities} descending.
#'
#' @importFrom DBI dbGetQuery
#' @importFrom tibble as_tibble
#' @export
entities_by_resource <- function(con) {

    sql <- "
        SELECT ds.name AS resource,
               COUNT(DISTINCT ee.entity_id) AS n_entities
        FROM data_source ds
        JOIN entity_evidence ee USING (source_id)
        GROUP BY ds.name
        ORDER BY n_entities DESC, ds.name
    "
    pg_query(con, sql)
}


#' Number of interactions per upstream resource
#'
#' Figure 1 Panel C.
#'
#' @param con A DBI connection.
#' @return A tibble with columns \code{resource} and \code{n_relations}.
#'
#' @importFrom DBI dbGetQuery
#' @export
interactions_by_resource <- function(con) {

    sql <- "
        SELECT ds.name AS resource,
               COUNT(*) AS n_relations
        FROM data_source ds
        JOIN relation r USING (source_id)
        GROUP BY ds.name
        ORDER BY n_relations DESC, ds.name
    "
    pg_query(con, sql)
}


#' Number of interactions per predicate (interaction type)
#'
#' Figure 1 Panel D. Uses the canonical \code{predicate} column on
#' \code{relation}.
#'
#' @param con A DBI connection.
#' @return A tibble with columns \code{interaction_type} and \code{n}.
#'
#' @importFrom DBI dbGetQuery
#' @export
interactions_by_type <- function(con) {

    sql <- "
        SELECT predicate AS interaction_type,
               COUNT(*)  AS n
        FROM relation
        GROUP BY predicate
        ORDER BY n DESC, predicate
    "
    pg_query(con, sql)
}


#' Number of annotation classes per resource
#'
#' Figure 1 Panel E. Counts distinct \code{annotation_id} contributed
#' by each resource via \code{entity_annotation_relation}.
#'
#' @param con A DBI connection.
#' @return A tibble with columns \code{resource} and \code{n_classes}.
#'
#' @importFrom DBI dbGetQuery
#' @export
annotation_classes_by_resource <- function(con) {

    sql <- "
        SELECT ds.name AS resource,
               COUNT(DISTINCT ea.annotation_id) AS n_classes
        FROM data_source ds
        JOIN entity_annotation_relation ea USING (source_id)
        GROUP BY ds.name
        ORDER BY n_classes DESC, ds.name
    "
    pg_query(con, sql)
}


#' Number of ontology terms per ontology
#'
#' Figure 1 Panel F.
#'
#' @param con A DBI connection.
#' @return A tibble with columns \code{ontology} and \code{n_terms}.
#'
#' @importFrom DBI dbGetQuery
#' @export
ontology_terms_by_ontology <- function(con) {

    sql <- "
        SELECT ontology, COUNT(*) AS n_terms
        FROM ontology_terms
        GROUP BY ontology
        ORDER BY n_terms DESC, ontology
    "
    pg_query(con, sql)
}
