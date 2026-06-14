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
    panel_id = "fig02-overview",
    n_types  = 8L
) {
    # The cycle-001 build's vocab_relation_predicate.interaction_class_id
    # only assigns predicates to Signaling, Transport and Other (the
    # finer classes — Allosteric, Ligand-receptor, TF-target,
    # Drug-target, Orthosteric, Maturation — have no predicate
    # mappings yet). To get the 5-7 categories the manuscript needs
    # without re-using a vocabulary that doesn't exist, we re-bucket
    # the predicates manually into broader functional categories:
    # Signaling, Transport, Interaction, Reaction, Association,
    # Membership, Other. These map directly to the "what kind of
    # relation is this" question the entity × interaction matrix
    # answers — finer than the cycle-001 class map but still derived
    # from named predicates.
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
            SELECT predicate,
                   CASE
                     WHEN predicate IN (
                       'controls', 'regulates',
                       'positively_regulates', 'negatively_regulates'
                     )                                  THEN 'Signaling'
                     WHEN predicate = 'transports'      THEN 'Transport'
                     WHEN predicate = 'interacts_with'  THEN 'Interaction'
                     WHEN predicate = 'has_participant' THEN 'Reaction'
                     WHEN predicate = 'associated_with' THEN 'Association'
                     WHEN predicate = 'has_member'      THEN 'Membership'
                     ELSE 'Other'
                   END AS class
            FROM pred_bm
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
