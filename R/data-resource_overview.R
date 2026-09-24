#' FR-007a faceted resource overview — data layer
#'
#' Six per-facet query functions plus a combiner. Each per-facet
#' function returns a long-format tibble with the columns
#' \code{resource}, \code{bar_type}, \code{category}, \code{n}. The
#' combiner concatenates them, adding a leading \code{facet} column
#' so the renderer can drive \code{ggplot2::facet_wrap()} directly.
#'
#' Bar types per the spec (FR-007a, Session 2026-06-02):
#' \describe{
#'   \item{\code{shared_unique}}{Top bar — two categories,
#'     \code{shared} and \code{unique}. Per-resource:
#'     \code{unique} = items only this resource contributes;
#'     \code{shared} = items also contributed by ≥ 1 other resource.
#'     For the \code{Total} row: \code{unique} = items contributed
#'     by exactly one resource overall, \code{shared} = ≥ 2.}
#'   \item{\code{major_class}}{Second bar — facet-specific class
#'     buckets (see per-function docs).}
#' }
#'
#' @name resource_overview_data_layer
NULL


#' SQL projection template for an FR-007a facet
#'
#' Returns the trailing 4-way \code{UNION ALL} that re-shapes a
#' per-resource and a total-by-class CTE into the long-format
#' \code{(resource, bar_type, category, n)} tibble shape. Hides the
#' duplication that previously led to a bug where the major-class
#' projection forgot to aggregate across the \code{shared_unique}
#' axis.
#'
#' Callers MUST define two CTEs upstream:
#' \describe{
#'   \item{\code{per_resource}}{Columns: \code{resource},
#'     \code{<class_col>}, \code{shared_unique}, \code{n}.}
#'   \item{\code{total_by_class}}{Columns: \code{<class_col>},
#'     \code{shared_unique}, \code{n}.}
#' }
#'
#' @param class_col Character: the major-class column name.
#'
#' @return Character: SQL fragment ready to follow a comma-separated
#'     pair of CTEs.
#'
#' @keywords internal
#' @noRd
resource_overview_projection <- function(class_col) {
    sprintf("
        SELECT 'Total' AS resource, 'shared_unique' AS bar_type,
               shared_unique AS category, SUM(n)::bigint AS n
        FROM   total_by_class GROUP BY shared_unique
        UNION ALL
        SELECT 'Total' AS resource, 'major_class' AS bar_type,
               %s AS category, SUM(n)::bigint AS n
        FROM   total_by_class GROUP BY %s
        UNION ALL
        SELECT resource, 'shared_unique' AS bar_type,
               shared_unique AS category, SUM(n)::bigint AS n
        FROM   per_resource GROUP BY resource, shared_unique
        UNION ALL
        SELECT resource, 'major_class' AS bar_type,
               %s AS category, SUM(n)::bigint AS n
        FROM   per_resource GROUP BY resource, %s
    ", class_col, class_col, class_col, class_col)
}


#' FR-007a Entities facet — per-resource + Total (bitmap path)
#'
#' Reads \code{facet_entity_bitmap} (cycle-001 pre-computed roaring
#' bitmaps) for the \code{source}, \code{entity_type} and
#' \code{chemical_class} facets, derives the FR-007a Entities major-
#' class bitmaps (proteins/genes/RNA, drugs, metabolites, lipids,
#' food compounds, other chemicals, other) via bitmap algebra, and
#' computes per-resource × class shared/unique counts via
#' \code{rb_and_cardinality} / \code{rb_andnot_cardinality}.
#'
#' Performance: ~150 ms end-to-end vs ~7.7 s for the equivalent
#' \code{entity_source_count}-scanning row-based query.
#' Numerically the bitmap totals match
#' \code{resources.entity_count} (the pre-computed source of truth)
#' rather than the legacy \code{entity_source_count} unnest path
#' which silently drops entities filtered out of the derive phase.
#'
#' Routed to \code{dev3} via the default
#' \code{\link{pg_query_panel}} dispatch.
#'
#' @param panel_id Character: panel identifier.
#'
#' @return Long-format tibble per \code{\link{resource_overview_data_layer}}.
#'
#' @importFrom DBI dbGetQuery
#' @export
resource_overview_entities <- function(panel_id = "database-content") {

    sql <- "
        WITH
        src AS (
            SELECT facet_value AS resource, entity_bitmap AS bm
            FROM   facet_entity_bitmap WHERE facet_name = 'source'
        ),
        etype AS (
            SELECT facet_value, entity_bitmap AS bm
            FROM   facet_entity_bitmap WHERE facet_name = 'entity_type'
        ),
        cclass AS (
            SELECT facet_value, entity_bitmap AS bm
            FROM   facet_entity_bitmap WHERE facet_name = 'chemical_class'
        ),
        all_entities AS (SELECT rb_or_agg(bm) AS bm FROM src),
        others AS (
            SELECT s.resource,
                   (SELECT rb_or_agg(s2.bm)
                    FROM src s2 WHERE s2.resource != s.resource) AS bm
            FROM src s
        ),
        exactly_one_per_src AS (
            SELECT s.resource, rb_andnot(s.bm, o.bm) AS bm
            FROM src s JOIN others o USING (resource)
        ),
        exactly_one AS (
            SELECT rb_or_agg(bm) AS bm FROM exactly_one_per_src
        ),
        multi_source AS (
            SELECT rb_andnot((SELECT bm FROM all_entities),
                             (SELECT bm FROM exactly_one)) AS bm
        ),
        m_proteins AS (
            SELECT rb_or_agg(bm) AS bm FROM etype
            WHERE facet_value IN ('Gene:MI:0250', 'Mirna:OM:0038',
                                  'Protein:MI:0326')
        ),
        m_drugs       AS (SELECT bm FROM cclass WHERE facet_value = 'drug'),
        m_metabolites AS (SELECT bm FROM cclass WHERE facet_value = 'metabolite'),
        m_lipids      AS (SELECT bm FROM cclass WHERE facet_value = 'lipid'),
        m_food        AS (SELECT bm FROM cclass WHERE facet_value = 'food'),
        m_chem        AS (
            SELECT bm FROM etype WHERE facet_value = 'Chemical:OM:0037'
        ),
        classified_chem AS (
            SELECT rb_or((SELECT bm FROM m_drugs),
                         rb_or((SELECT bm FROM m_metabolites),
                               rb_or((SELECT bm FROM m_lipids),
                                     (SELECT bm FROM m_food)))) AS bm
        ),
        m_other_chem AS (
            SELECT rb_andnot((SELECT bm FROM m_chem),
                             (SELECT bm FROM classified_chem)) AS bm
        ),
        all_classified AS (
            SELECT rb_or((SELECT bm FROM m_proteins),
                         rb_or((SELECT bm FROM classified_chem),
                               (SELECT bm FROM m_other_chem))) AS bm
        ),
        m_other AS (
            SELECT rb_andnot((SELECT bm FROM all_entities),
                             (SELECT bm FROM all_classified)) AS bm
        ),
        classes(name, bm) AS (
            VALUES
                ('proteins/genes/RNA', (SELECT bm FROM m_proteins)),
                ('drugs',              (SELECT bm FROM m_drugs)),
                ('metabolites',        (SELECT bm FROM m_metabolites)),
                ('lipids',             (SELECT bm FROM m_lipids)),
                ('food compounds',     (SELECT bm FROM m_food)),
                ('other chemicals',    (SELECT bm FROM m_other_chem)),
                ('other',              (SELECT bm FROM m_other))
        ),
        per_res_cell AS (
            SELECT s.resource,
                   c.name AS major_class,
                   rb_andnot_cardinality(rb_and(s.bm, c.bm), o.bm)::bigint
                       AS unique_n,
                   rb_and_cardinality(rb_and(s.bm, c.bm), o.bm)::bigint
                       AS shared_n
            FROM src s
            CROSS JOIN classes c
            JOIN others o USING (resource)
        ),
        total_cell AS (
            SELECT c.name AS major_class,
                   rb_andnot_cardinality(c.bm,
                       (SELECT bm FROM multi_source))::bigint AS unique_n,
                   rb_and_cardinality(c.bm,
                       (SELECT bm FROM multi_source))::bigint AS shared_n
            FROM classes c
        )
        -- Total row, shared_unique bar
        SELECT 'Total' AS resource, 'shared_unique' AS bar_type,
               'unique' AS category, SUM(unique_n)::bigint AS n
        FROM total_cell
        UNION ALL
        SELECT 'Total', 'shared_unique', 'shared', SUM(shared_n)::bigint
        FROM total_cell
        UNION ALL
        -- Total row, major_class bar
        SELECT 'Total', 'major_class', major_class,
               (unique_n + shared_n)::bigint
        FROM total_cell
        UNION ALL
        -- Per-resource, shared_unique bar
        SELECT resource, 'shared_unique', 'unique',
               SUM(unique_n)::bigint
        FROM per_res_cell GROUP BY resource
        UNION ALL
        SELECT resource, 'shared_unique', 'shared',
               SUM(shared_n)::bigint
        FROM per_res_cell GROUP BY resource
        UNION ALL
        -- Per-resource, major_class bar
        SELECT resource, 'major_class', major_class,
               (unique_n + shared_n)::bigint
        FROM per_res_cell
        WHERE unique_n + shared_n > 0
    "
    pg_query_panel(panel_id, sql)
}


#' FR-007a Interactions facet — per-resource + Total (bitmap path)
#'
#' Reads \code{facet_relation_bitmap} for the \code{source} and
#' \code{predicate} facets, joins predicate → interaction class
#' (Signaling / Transport / Other), and computes per-resource ×
#' class shared/unique counts via roaringbitmap algebra.
#'
#' Performance: ~150 ms vs ~44 s for the equivalent
#' \code{relation_evidence_relation}-scanning row-based query.
#'
#' @inheritParams resource_overview_entities
#'
#' @return Long-format tibble per \code{\link{resource_overview_data_layer}}.
#'
#' @importFrom DBI dbGetQuery
#' @export
resource_overview_interactions <- function(panel_id = "database-content") {

    sql <- "
        WITH
        src AS (
            SELECT facet_value AS resource, relation_bitmap AS bm
            FROM   facet_relation_bitmap WHERE facet_name = 'source'
        ),
        pred_bm AS (
            SELECT facet_value AS predicate, relation_bitmap AS bm
            FROM   facet_relation_bitmap WHERE facet_name = 'predicate'
        ),
        -- The cycle-001 build's vocab_relation_predicate
        -- .interaction_class_id only buckets predicates into
        -- Signaling / Transport / Other. To get the 7-class
        -- vocabulary the manuscript figure uses (Signaling,
        -- Transport, Interaction, Reaction, Association,
        -- Membership, Other), re-bucket predicates by name using
        -- the same CASE statement as entity_by_interaction_type.
        pred_class AS (
            SELECT facet_value AS predicate,
                   CASE
                     WHEN facet_value IN (
                       'controls', 'regulates',
                       'positively_regulates', 'negatively_regulates'
                     )                                  THEN 'Signaling'
                     WHEN facet_value = 'transports'    THEN 'Transport'
                     WHEN facet_value = 'interacts_with' THEN 'Interaction'
                     WHEN facet_value = 'has_participant' THEN 'Reaction'
                     WHEN facet_value = 'associated_with' THEN 'Association'
                     WHEN facet_value = 'has_member'    THEN 'Membership'
                     ELSE 'Other'
                   END AS class
            FROM   facet_relation_bitmap
            WHERE  facet_name = 'predicate'
        ),
        class_bm AS (
            SELECT pc.class, rb_or_agg(pb.bm) AS bm
            FROM   pred_bm pb JOIN pred_class pc USING (predicate)
            GROUP  BY pc.class
        ),
        all_rels AS (SELECT rb_or_agg(bm) AS bm FROM src),
        others AS (
            SELECT s.resource,
                   (SELECT rb_or_agg(s2.bm)
                    FROM src s2 WHERE s2.resource != s.resource) AS bm
            FROM src s
        ),
        exactly_one_per_src AS (
            SELECT s.resource, rb_andnot(s.bm, o.bm) AS bm
            FROM src s JOIN others o USING (resource)
        ),
        exactly_one AS (
            SELECT rb_or_agg(bm) AS bm FROM exactly_one_per_src
        ),
        multi_source AS (
            SELECT rb_andnot((SELECT bm FROM all_rels),
                             (SELECT bm FROM exactly_one)) AS bm
        ),
        per_res_cell AS (
            SELECT s.resource,
                   c.class AS interaction_class,
                   rb_andnot_cardinality(rb_and(s.bm, c.bm), o.bm)::bigint
                       AS unique_n,
                   rb_and_cardinality(rb_and(s.bm, c.bm), o.bm)::bigint
                       AS shared_n
            FROM src s
            CROSS JOIN class_bm c
            JOIN others o USING (resource)
        ),
        total_cell AS (
            SELECT c.class AS interaction_class,
                   rb_andnot_cardinality(c.bm,
                       (SELECT bm FROM multi_source))::bigint AS unique_n,
                   rb_and_cardinality(c.bm,
                       (SELECT bm FROM multi_source))::bigint AS shared_n
            FROM class_bm c
        )
        SELECT 'Total' AS resource, 'shared_unique' AS bar_type,
               'unique' AS category, SUM(unique_n)::bigint AS n
        FROM total_cell
        UNION ALL
        SELECT 'Total', 'shared_unique', 'shared', SUM(shared_n)::bigint
        FROM total_cell
        UNION ALL
        SELECT 'Total', 'major_class', interaction_class,
               (unique_n + shared_n)::bigint
        FROM total_cell
        UNION ALL
        SELECT resource, 'shared_unique', 'unique',
               SUM(unique_n)::bigint
        FROM per_res_cell GROUP BY resource
        UNION ALL
        SELECT resource, 'shared_unique', 'shared',
               SUM(shared_n)::bigint
        FROM per_res_cell GROUP BY resource
        UNION ALL
        SELECT resource, 'major_class', interaction_class,
               (unique_n + shared_n)::bigint
        FROM per_res_cell
        WHERE unique_n + shared_n > 0
    "
    pg_query_panel(panel_id, sql)
}


#' FR-007a Associations facet — per-resource + Total
#'
#' Reads \code{entity_ontology_term} and buckets by
#' \code{ontology_prefix} into broad association classes
#' (\code{disease}, \code{phenotype}, \code{function (GO)},
#' \code{anatomy}, \code{chemical taxonomy}, \code{enzymes /
#' pathways}, \code{other}). The \code{shared_unique} split uses
#' \code{cardinality(sources)} on the table's \code{text[]} sources
#' column directly — no source-id join needed.
#'
#' @inheritParams resource_overview_entities
#'
#' @return Long-format tibble per \code{\link{resource_overview_data_layer}}.
#'
#' @importFrom DBI dbGetQuery
#' @export
resource_overview_associations <- function(panel_id = "database-content") {

    sql <- sprintf("
        WITH classified AS (
            SELECT
                eot.term_entity_id,
                eot.sources,
                CASE
                    WHEN eot.ontology_prefix = 'mondo' THEN 'disease'
                    WHEN eot.ontology_prefix = 'hp' THEN 'phenotype'
                    WHEN eot.ontology_prefix IN ('go', 'gene_ontology')
                        THEN 'function (GO)'
                    WHEN eot.ontology_prefix = 'uberon' THEN 'anatomy'
                    WHEN eot.ontology_prefix IN ('chebi', 'chemontid',
                                                 'slm')
                        THEN 'chemical taxonomy'
                    WHEN eot.ontology_prefix IN ('ec',
                                                 'enzyme_classification',
                                                 'reactome_pathways')
                        THEN 'enzymes / pathways'
                    ELSE 'other'
                END AS major_class
            FROM   entity_ontology_term eot
        ),
        per_resource AS (
            SELECT
                src AS resource,
                c.major_class,
                CASE WHEN cardinality(c.sources) = 1
                     THEN 'unique' ELSE 'shared' END AS shared_unique,
                COUNT(*)::bigint AS n
            FROM   classified c
            CROSS  JOIN LATERAL unnest(c.sources) AS src
            GROUP  BY src, c.major_class, shared_unique
        ),
        total_by_class AS (
            SELECT
                c.major_class,
                CASE WHEN cardinality(c.sources) = 1
                     THEN 'unique' ELSE 'shared' END AS shared_unique,
                COUNT(*)::bigint AS n
            FROM   classified c
            GROUP  BY c.major_class, shared_unique
        )
        %s
    ", resource_overview_projection("major_class"))

    pg_query_panel(panel_id, sql)
}


#' FR-007a Identifiers facet — major-ID-types variant
#'
#' Per-resource identifier counts split by shared/unique and by major
#' identifier type (top-N from \code{vocab_identifier_type}) with a
#' \code{Misc} long-tail collapse.
#'
#' @inheritParams resource_overview_entities
#' @param top_n Integer: how many ID types to keep before the
#'     \code{Misc} collapse. Default \code{8L}.
#'
#' @return Long-format tibble per \code{\link{resource_overview_data_layer}}.
#'
#' @importFrom DBI dbGetQuery
#' @export
resource_overview_identifiers_major <- function(panel_id = "database-content",
                                     top_n = 8L) {

    sql <- sprintf("
        WITH id_sources AS (
            SELECT
                eei.identifier_id,
                ARRAY_AGG(DISTINCT eei.source_id) AS source_list,
                COUNT(DISTINCT eei.source_id) AS source_count
            FROM   entity_evidence_identifier eei
            GROUP  BY eei.identifier_id
        ),
        ranked_types AS (
            SELECT vit.name AS id_type, COUNT(*) AS n,
                   ROW_NUMBER() OVER (ORDER BY COUNT(*) DESC) AS rk
            FROM   identifier_evidence ie
            JOIN   vocab_identifier_type vit
                   ON vit.identifier_type_id = ie.identifier_type_id
            GROUP  BY vit.name
        ),
        id_class AS (
            SELECT DISTINCT
                ie.identifier_id,
                CASE WHEN rt.rk <= %d
                     THEN rt.id_type ELSE 'Misc' END AS id_class
            FROM   identifier_evidence ie
            JOIN   vocab_identifier_type vit
                   ON vit.identifier_type_id = ie.identifier_type_id
            JOIN   ranked_types rt ON rt.id_type = vit.name
        ),
        per_resource AS (
            SELECT
                ds.name AS resource,
                ic.id_class,
                CASE WHEN ids.source_count = 1
                     THEN 'unique' ELSE 'shared' END AS shared_unique,
                COUNT(*)::bigint AS n
            FROM   id_sources ids
            CROSS  JOIN LATERAL unnest(ids.source_list) AS src
            JOIN   data_source ds ON ds.source_id = src
            JOIN   id_class ic USING (identifier_id)
            GROUP  BY ds.name, ic.id_class, shared_unique
        ),
        total_by_class AS (
            SELECT
                ic.id_class,
                CASE WHEN ids.source_count = 1
                     THEN 'unique' ELSE 'shared' END AS shared_unique,
                COUNT(*)::bigint AS n
            FROM   id_sources ids
            JOIN   id_class ic USING (identifier_id)
            GROUP  BY ic.id_class, shared_unique
        )
        %s
    ", as.integer(top_n), resource_overview_projection("id_class"))

    pg_query_panel(panel_id, sql)
}


#' FR-007a Identifiers facet — authoritative-vs-foreign variant
#'
#' Each identifier is either the entity's \emph{authoritative} ID
#' (matches \code{entity.canonical_identifier_type_id}) or a
#' \emph{foreign} cross-reference. Per-resource counts of each.
#'
#' @inheritParams resource_overview_entities
#'
#' @return Long-format tibble per \code{\link{resource_overview_data_layer}}.
#'
#' @importFrom DBI dbGetQuery
#' @export
resource_overview_identifiers_authoritative <- function(panel_id = "database-content") {

    sql <- sprintf("
        WITH id_sources AS (
            SELECT
                eei.identifier_id,
                ARRAY_AGG(DISTINCT eei.source_id) AS source_list,
                COUNT(DISTINCT eei.source_id) AS source_count
            FROM   entity_evidence_identifier eei
            GROUP  BY eei.identifier_id
        ),
        id_class AS (
            SELECT DISTINCT
                ie.identifier_id,
                CASE WHEN e.canonical_identifier_type_id = ie.identifier_type_id
                     THEN 'authoritative' ELSE 'foreign' END AS id_class
            FROM   identifier_evidence ie
            JOIN   entity_evidence_identifier eei USING (identifier_id)
            JOIN   entity_evidence_resolution eer
                   ON eer.entity_evidence_id = eei.entity_evidence_id
            JOIN   entity e ON e.entity_id = eer.entity_id
        ),
        per_resource AS (
            SELECT
                ds.name AS resource,
                ic.id_class,
                CASE WHEN ids.source_count = 1
                     THEN 'unique' ELSE 'shared' END AS shared_unique,
                COUNT(*)::bigint AS n
            FROM   id_sources ids
            CROSS  JOIN LATERAL unnest(ids.source_list) AS src
            JOIN   data_source ds ON ds.source_id = src
            JOIN   id_class ic USING (identifier_id)
            GROUP  BY ds.name, ic.id_class, shared_unique
        ),
        total_by_class AS (
            SELECT
                ic.id_class,
                CASE WHEN ids.source_count = 1
                     THEN 'unique' ELSE 'shared' END AS shared_unique,
                COUNT(*)::bigint AS n
            FROM   id_sources ids
            JOIN   id_class ic USING (identifier_id)
            GROUP  BY ic.id_class, shared_unique
        )
        %s
    ", resource_overview_projection("id_class"))

    pg_query_panel(panel_id, sql)
}


#' FR-007a Structures facet — per-resource + Total (dev4, bitmap path)
#'
#' Reads \code{facet_entity_bitmap} for the \code{source} and
#' \code{structural_specificity} facets on dev4 (both facets are
#' pre-computed by the cycle-001 derive phase; the latter is
#' \code{dev4}-only per FR-030). Specificity levels:
#' \code{stereospecific} / \code{cis_trans_only} /
#' \code{constitution_only} / \code{variable_constitution} /
#' \code{unknown_constitution} / \code{no_structure}. The
#' \code{no_structure} bucket is rendered as-counted (dev4 is
#' pre-T020 chemical fallback so ~44 % of chemicals are still
#' structure-less hashes).
#'
#' Performance: ~150 ms vs ~8.4 s for the equivalent
#' \code{metabo_entity_structural_specificity}-scanning row query.
#'
#' @inheritParams resource_overview_entities
#'
#' @return Long-format tibble per \code{\link{resource_overview_data_layer}}.
#'
#' @importFrom DBI dbGetQuery
#' @export
resource_overview_structures <- function(panel_id = "database-content") {

    sql <- "
        WITH
        src AS (
            SELECT facet_value AS resource, entity_bitmap AS bm
            FROM   facet_entity_bitmap WHERE facet_name = 'source'
        ),
        spec AS (
            SELECT facet_value AS class, entity_bitmap AS bm
            FROM   facet_entity_bitmap
            WHERE  facet_name = 'structural_specificity'
        ),
        all_entities AS (SELECT rb_or_agg(bm) AS bm FROM src),
        others AS (
            SELECT s.resource,
                   (SELECT rb_or_agg(s2.bm)
                    FROM src s2 WHERE s2.resource != s.resource) AS bm
            FROM src s
        ),
        exactly_one_per_src AS (
            SELECT s.resource, rb_andnot(s.bm, o.bm) AS bm
            FROM src s JOIN others o USING (resource)
        ),
        exactly_one AS (
            SELECT rb_or_agg(bm) AS bm FROM exactly_one_per_src
        ),
        multi_source AS (
            SELECT rb_andnot((SELECT bm FROM all_entities),
                             (SELECT bm FROM exactly_one)) AS bm
        ),
        per_res_cell AS (
            SELECT s.resource, sp.class AS major_class,
                   rb_andnot_cardinality(rb_and(s.bm, sp.bm), o.bm)::bigint
                       AS unique_n,
                   rb_and_cardinality(rb_and(s.bm, sp.bm), o.bm)::bigint
                       AS shared_n
            FROM src s
            CROSS JOIN spec sp
            JOIN others o USING (resource)
        ),
        total_cell AS (
            SELECT sp.class AS major_class,
                   rb_andnot_cardinality(sp.bm,
                       (SELECT bm FROM multi_source))::bigint AS unique_n,
                   rb_and_cardinality(sp.bm,
                       (SELECT bm FROM multi_source))::bigint AS shared_n
            FROM spec sp
        )
        SELECT 'Total' AS resource, 'shared_unique' AS bar_type,
               'unique' AS category, SUM(unique_n)::bigint AS n
        FROM total_cell
        UNION ALL
        SELECT 'Total', 'shared_unique', 'shared', SUM(shared_n)::bigint
        FROM total_cell
        UNION ALL
        SELECT 'Total', 'major_class', major_class,
               (unique_n + shared_n)::bigint
        FROM total_cell
        UNION ALL
        SELECT resource, 'shared_unique', 'unique',
               SUM(unique_n)::bigint
        FROM per_res_cell GROUP BY resource
        UNION ALL
        SELECT resource, 'shared_unique', 'shared',
               SUM(shared_n)::bigint
        FROM per_res_cell GROUP BY resource
        UNION ALL
        SELECT resource, 'major_class', major_class,
               (unique_n + shared_n)::bigint
        FROM per_res_cell
        WHERE unique_n + shared_n > 0
    "
    pg_query_panel(panel_id, sql, facet = "structures")
}


#' FR-007a Literature references facet — per-resource + Total
#'
#' Reads \code{relation_evidence_annotation} joined to
#' \code{annotation} where the \code{term} matches a literature
#' identifier vocabulary (\code{Pubmed:MI:0446} /
#' \code{Pubmed Central:MI:1042} / \code{Doi:MI:0574}). Major class
#' is the literature-id vocabulary. Shared/unique computed on
#' \code{annotation_key}.
#'
#' @inheritParams resource_overview_entities
#'
#' @return Long-format tibble per \code{\link{resource_overview_data_layer}}.
#'
#' @importFrom DBI dbGetQuery
#' @export
resource_overview_literature <- function(panel_id = "database-content") {

    sql <- sprintf("
        WITH lit_ann AS (
            SELECT a.annotation_key,
                   CASE
                       WHEN a.term LIKE 'Pubmed:%%' THEN 'PubMed'
                       WHEN a.term LIKE 'Pubmed Central:%%' THEN 'PMC'
                       WHEN a.term LIKE 'Doi:%%' THEN 'DOI'
                       ELSE NULL
                   END AS major_class
            FROM   annotation a
            WHERE  a.term LIKE 'Pubmed:%%'
                OR a.term LIKE 'Pubmed Central:%%'
                OR a.term LIKE 'Doi:%%'
        ),
        lit_sources AS (
            SELECT
                rea.annotation_key,
                ARRAY_AGG(DISTINCT rea.source_id) AS source_list,
                COUNT(DISTINCT rea.source_id) AS source_count
            FROM   relation_evidence_annotation rea
            JOIN   lit_ann USING (annotation_key)
            GROUP  BY rea.annotation_key
        ),
        per_resource AS (
            SELECT
                ds.name AS resource,
                la.major_class,
                CASE WHEN ls.source_count = 1
                     THEN 'unique' ELSE 'shared' END AS shared_unique,
                COUNT(*)::bigint AS n
            FROM   lit_sources ls
            CROSS  JOIN LATERAL unnest(ls.source_list) AS src
            JOIN   data_source ds ON ds.source_id = src
            JOIN   lit_ann la USING (annotation_key)
            GROUP  BY ds.name, la.major_class, shared_unique
        ),
        total_by_class AS (
            SELECT
                la.major_class,
                CASE WHEN ls.source_count = 1
                     THEN 'unique' ELSE 'shared' END AS shared_unique,
                COUNT(*)::bigint AS n
            FROM   lit_sources ls
            JOIN   lit_ann la USING (annotation_key)
            GROUP  BY la.major_class, shared_unique
        )
        %s
    ", resource_overview_projection("major_class"))

    pg_query_panel(panel_id, sql)
}


#' Resource short-label lookup
#'
#' Reads \code{resources.resource_id} + \code{resources.resource_short}
#' (with \code{resource_id} as a fallback) so figure labels show
#' "ChEMBL", "RaMP", "SwissLipids" rather than the lowercase
#' \code{chembl} / \code{ramp} / \code{swisslipids} slugs.
#'
#' @inheritParams resource_overview_entities
#'
#' @return Named character vector mapping \code{resource_id} →
#'     display label.
#'
#' @importFrom DBI dbGetQuery
#' @export
resources_label_map <- function(panel_id = "database-content") {

    sql <- "
        SELECT resource_id,
               COALESCE(NULLIF(resource_short, ''), resource_id) AS label
        FROM   resources
    "
    rows <- pg_query_panel(panel_id, sql)
    stats::setNames(rows$label, rows$resource_id)
}


#' Combine FR-007a facets into a single long-format tibble
#'
#' Runs the per-facet queries and concatenates the results, adding a
#' leading \code{facet} column so the renderer can drive
#' \code{ggplot2::facet_wrap()}.
#'
#' @param panel_id Character: panel identifier (forwarded to each
#'     per-facet query).
#' @param identifiers Character: which Identifiers variant to
#'     include. \code{"major"} (default, top-N ID types + Misc) or
#'     \code{"authoritative"} (authoritative-vs-foreign).
#'
#' @return A tibble with columns \code{facet}, \code{resource},
#'     \code{bar_type}, \code{category}, \code{n}.
#'
#' @importFrom dplyr bind_rows
#' @importFrom logger log_info
#' @export
resource_overview_overview <- function(panel_id = "database-content",
                            identifiers = c("major", "authoritative")) {

    identifiers <- match.arg(identifiers)

    facet_fn <- list(
        Entities       = resource_overview_entities,
        Associations   = resource_overview_associations,
        Interactions   = resource_overview_interactions,
        Identifiers    = if (identifiers == "major") resource_overview_identifiers_major
                         else resource_overview_identifiers_authoritative,
        Structures     = resource_overview_structures,
        Literature     = resource_overview_literature
    )

    out_list <- list()
    for (facet_name in names(facet_fn)) {
        logger::log_info("FR-007a {facet_name} facet")
        rows <- facet_fn[[facet_name]](panel_id)
        rows$facet <- facet_name
        out_list[[facet_name]] <- rows
    }
    combined <- dplyr::bind_rows(out_list)

    # Map resource_id → resource_short label (Total stays "Total").
    labels <- resources_label_map(panel_id)
    combined$resource_label <- ifelse(
        combined$resource == "Total",
        "Total",
        unname(labels[combined$resource]) %|na|% combined$resource
    )

    # Magnitude band — based on each resource's shared/unique total
    # in the Entities facet. Used by the renderer to split the panel
    # into row bands so a 1-resource-with-1M-entities tail doesn't
    # squash the rest onto an invisible scale.
    #
    # Critical: RPostgres returns BIGINT as bit64::integer64. cut()
    # and max() on integer64 silently misclassify (the bit pattern
    # gets reinterpreted as ~1e-317). Coerce to double FIRST.
    entities_total <- combined[
        combined$facet == "Entities" &
            combined$bar_type == "shared_unique",
        c("resource", "n"), drop = FALSE
    ]
    entities_total$n <- as.numeric(entities_total$n)
    entities_total <- stats::aggregate(
        n ~ resource, data = entities_total, FUN = sum
    )

    band_of <- function(value) {
        cut(
            as.numeric(value),
            breaks = c(-Inf, 1e4, 1e5, Inf),
            labels = c(
                "small (< 10K)",
                "medium (10K – 100K)",
                "large (>= 100K)"
            ),
            right = FALSE
        )
    }
    band_lookup <- stats::setNames(
        as.character(band_of(entities_total$n)),
        entities_total$resource
    )
    # Total row always sits in the "large" band so it heads the
    # top-magnitude row.
    band_lookup["Total"] <- "large (>= 100K)"
    combined$magnitude_band <- factor(
        ifelse(
            is.na(band_lookup[combined$resource]),
            "small (< 10K)",
            band_lookup[combined$resource]
        ),
        levels = c(
            "large (>= 100K)",
            "medium (10K – 100K)",
            "small (< 10K)"
        )
    )

    combined
}


#' Local replacement-for-NA helper used by \code{resource_overview_overview}
#'
#' @keywords internal
#' @noRd
`%|na|%` <- function(x, y) ifelse(is.na(x), y, x)
