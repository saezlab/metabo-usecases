#' FR-007d entity × interaction-type matrix data layer
#'
#' Cross-tab of relation \code{participant_type} ×
#' \code{vocab_interaction_class}, computed via bitmap intersection
#' between \code{facet_relation_bitmap} facets
#' (\code{participant_type} and \code{predicate} → class). Returns
#' counts for the top-\code{n_types} participant types by relation
#' count.
#'
#' Performance: ~150 ms via roaringbitmap algebra; no row scans.
#'
#' @param panel_id Character: panel identifier.
#' @param n_types Integer: how many top participant types to keep
#'     (default \code{8L}).
#'
#' @return Tibble with \code{participant_type},
#'     \code{interaction_class}, \code{n}.
#'
#' @importFrom DBI dbGetQuery
#' @export
fr007d_entity_x_interaction <- function(
    panel_id = "fig01-overview",
    n_types  = 8L
) {
    sql <- sprintf("
        WITH
        pt_all AS (
            SELECT facet_value AS participant_type, relation_bitmap AS bm,
                   relation_count
            FROM   facet_relation_bitmap
            WHERE  facet_name = 'participant_type'
        ),
        pt AS (
            SELECT * FROM pt_all
            ORDER  BY relation_count DESC
            LIMIT  %d
        ),
        pred_bm AS (
            SELECT facet_value AS predicate, relation_bitmap AS bm
            FROM   facet_relation_bitmap WHERE facet_name = 'predicate'
        ),
        pred_class AS (
            SELECT vrp.name AS predicate,
                   COALESCE(vic.name, 'Other') AS class
            FROM   vocab_relation_predicate vrp
            LEFT   JOIN vocab_interaction_class vic
                   ON vic.interaction_class_id = vrp.interaction_class_id
        ),
        class_bm AS (
            SELECT pc.class, rb_or_agg(pb.bm) AS bm
            FROM   pred_bm pb JOIN pred_class pc USING (predicate)
            GROUP  BY pc.class
        )
        SELECT pt.participant_type,
               cb.class AS interaction_class,
               rb_and_cardinality(pt.bm, cb.bm)::bigint AS n
        FROM   pt CROSS JOIN class_bm cb
        ORDER  BY pt.relation_count DESC, cb.class
    ", as.integer(n_types))

    pg_query_panel(panel_id, sql)
}
