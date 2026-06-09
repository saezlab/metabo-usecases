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
#' @name fr007a_data_layer
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
fr007a_projection <- function(class_col) {
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


#' FR-007a Entities facet — per-resource + Total
#'
#' Reads \code{entity_source_count} (cycle-001 pre-computed),
#' \code{data_source}, \code{vocab_chemical_class},
#' \code{vocab_entity_type}; buckets entities into the FR-007a
#' Entities major-class list (proteins/genes/RNA, drugs,
#' metabolites, lipids, food compounds, other chemicals, other).
#' Routed to \code{dev3} via the default
#' \code{\link{pg_query_panel}} dispatch.
#'
#' @param panel_id Character: panel identifier.
#'
#' @return Long-format tibble per \code{\link{fr007a_data_layer}}.
#'
#' @importFrom DBI dbGetQuery
#' @export
fr007a_entities <- function(panel_id = "fig01-overview") {

    sql <- sprintf("
        WITH entity_class AS (
            SELECT
                e.entity_id,
                CASE
                    WHEN vet.name IN ('Gene:MI:0250', 'Mirna:OM:0038',
                                      'Protein:MI:0326')
                        THEN 'proteins/genes/RNA'
                    WHEN vcc.name = 'drug'       THEN 'drugs'
                    WHEN vcc.name = 'metabolite' THEN 'metabolites'
                    WHEN vcc.name = 'lipid'      THEN 'lipids'
                    WHEN vcc.name = 'food'       THEN 'food compounds'
                    WHEN vet.name = 'Chemical:OM:0037'
                         AND e.chemical_class_id IS NULL
                        THEN 'other chemicals'
                    ELSE 'other'
                END AS major_class
            FROM   entity e
            JOIN   vocab_entity_type vet
                   ON vet.entity_type_id = e.entity_type_id
            LEFT   JOIN vocab_chemical_class vcc
                   ON vcc.chemical_class_id = e.chemical_class_id
        ),
        per_resource AS (
            SELECT
                ds.name AS resource,
                ec.major_class,
                CASE WHEN esc.source_count = 1
                     THEN 'unique' ELSE 'shared' END AS shared_unique,
                COUNT(*)::bigint AS n
            FROM   entity_source_count esc
            CROSS  JOIN LATERAL unnest(esc.source_list) AS src
            JOIN   data_source ds ON ds.source_id = src
            JOIN   entity_class ec USING (entity_id)
            GROUP  BY ds.name, ec.major_class, shared_unique
        ),
        total_by_class AS (
            SELECT
                ec.major_class,
                CASE WHEN esc.source_count = 1
                     THEN 'unique' ELSE 'shared' END AS shared_unique,
                COUNT(*)::bigint AS n
            FROM   entity_source_count esc
            JOIN   entity_class ec USING (entity_id)
            GROUP  BY ec.major_class, shared_unique
        )
        %s
    ", fr007a_projection("major_class"))

    pg_query_panel(panel_id, sql)
}


#' FR-007a Interactions facet — per-resource + Total
#'
#' Per-resource interaction counts split by shared/unique (via
#' \code{relation_evidence_relation.source_id}) and by coarse
#' interaction class (\code{vocab_interaction_class}: Signaling /
#' Transport / Other). Cycle-001 predicate vocabulary populates only
#' those three classes today; finer breakdowns collapse to
#' \code{Other} (handover note).
#'
#' @inheritParams fr007a_entities
#'
#' @return Long-format tibble per \code{\link{fr007a_data_layer}}.
#'
#' @importFrom DBI dbGetQuery
#' @export
fr007a_interactions <- function(panel_id = "fig01-overview") {

    sql <- sprintf("
        WITH rel_class AS (
            SELECT
                r.relation_id,
                COALESCE(vic.name, 'Other') AS interaction_class
            FROM   relation r
            JOIN   vocab_relation_predicate vrp
                   ON vrp.relation_predicate_id = r.predicate_id
            LEFT   JOIN vocab_interaction_class vic
                   ON vic.interaction_class_id = vrp.interaction_class_id
        ),
        rel_sources AS (
            SELECT
                rer.relation_id,
                ARRAY_AGG(DISTINCT rer.source_id) AS source_list,
                COUNT(DISTINCT rer.source_id) AS source_count
            FROM   relation_evidence_relation rer
            GROUP  BY rer.relation_id
        ),
        per_resource AS (
            SELECT
                ds.name AS resource,
                rc.interaction_class,
                CASE WHEN rs.source_count = 1
                     THEN 'unique' ELSE 'shared' END AS shared_unique,
                COUNT(*)::bigint AS n
            FROM   rel_sources rs
            CROSS  JOIN LATERAL unnest(rs.source_list) AS src
            JOIN   data_source ds ON ds.source_id = src
            JOIN   rel_class rc USING (relation_id)
            GROUP  BY ds.name, rc.interaction_class, shared_unique
        ),
        total_by_class AS (
            SELECT
                rc.interaction_class,
                CASE WHEN rs.source_count = 1
                     THEN 'unique' ELSE 'shared' END AS shared_unique,
                COUNT(*)::bigint AS n
            FROM   rel_sources rs
            JOIN   rel_class rc USING (relation_id)
            GROUP  BY rc.interaction_class, shared_unique
        )
        %s
    ", fr007a_projection("interaction_class"))

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
#' @inheritParams fr007a_entities
#'
#' @return Long-format tibble per \code{\link{fr007a_data_layer}}.
#'
#' @importFrom DBI dbGetQuery
#' @export
fr007a_associations <- function(panel_id = "fig01-overview") {

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
    ", fr007a_projection("major_class"))

    pg_query_panel(panel_id, sql)
}


#' FR-007a Identifiers facet — major-ID-types variant
#'
#' Per-resource identifier counts split by shared/unique and by major
#' identifier type (top-N from \code{vocab_identifier_type}) with a
#' \code{Misc} long-tail collapse.
#'
#' @inheritParams fr007a_entities
#' @param top_n Integer: how many ID types to keep before the
#'     \code{Misc} collapse. Default \code{8L}.
#'
#' @return Long-format tibble per \code{\link{fr007a_data_layer}}.
#'
#' @importFrom DBI dbGetQuery
#' @export
fr007a_identifiers_major <- function(panel_id = "fig01-overview",
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
    ", as.integer(top_n), fr007a_projection("id_class"))

    pg_query_panel(panel_id, sql)
}


#' FR-007a Identifiers facet — authoritative-vs-foreign variant
#'
#' Each identifier is either the entity's \emph{authoritative} ID
#' (matches \code{entity.canonical_identifier_type_id}) or a
#' \emph{foreign} cross-reference. Per-resource counts of each.
#'
#' @inheritParams fr007a_entities
#'
#' @return Long-format tibble per \code{\link{fr007a_data_layer}}.
#'
#' @importFrom DBI dbGetQuery
#' @export
fr007a_identifiers_authoritative <- function(panel_id = "fig01-overview") {

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
    ", fr007a_projection("id_class"))

    pg_query_panel(panel_id, sql)
}


#' FR-007a Structures facet — per-resource + Total (dev4)
#'
#' Reads \code{metabo_entity_structural_specificity} +
#' \code{metabo_vocab_structural_specificity} on dev4 (the tables
#' are dev4-only per FR-030); shared/unique computed from
#' \code{entity_source_count} on the same deployment. Major-class
#' bar is the six specificity levels themselves
#' (\code{stereospecific} / \code{cis_trans_only} /
#' \code{constitution_only} / \code{variable_constitution} /
#' \code{unknown_constitution} / \code{no_structure}). The
#' \code{no_structure} bucket is rendered as-counted (dev4 is
#' pre-T020 chemical fallback so ~44 % of chemicals are still
#' structure-less hashes).
#'
#' @inheritParams fr007a_entities
#'
#' @return Long-format tibble per \code{\link{fr007a_data_layer}}.
#'
#' @importFrom DBI dbGetQuery
#' @export
fr007a_structures <- function(panel_id = "fig01-overview") {

    sql <- sprintf("
        WITH structure_class AS (
            SELECT mess.entity_id, mvss.name AS major_class
            FROM   metabo_entity_structural_specificity mess
            JOIN   metabo_vocab_structural_specificity mvss
                   ON mvss.structural_specificity_id
                      = mess.structural_specificity_id
        ),
        per_resource AS (
            SELECT
                ds.name AS resource,
                sc.major_class,
                CASE WHEN esc.source_count = 1
                     THEN 'unique' ELSE 'shared' END AS shared_unique,
                COUNT(*)::bigint AS n
            FROM   structure_class sc
            JOIN   entity_source_count esc USING (entity_id)
            CROSS  JOIN LATERAL unnest(esc.source_list) AS src
            JOIN   data_source ds ON ds.source_id = src
            GROUP  BY ds.name, sc.major_class, shared_unique
        ),
        total_by_class AS (
            SELECT
                sc.major_class,
                CASE WHEN esc.source_count = 1
                     THEN 'unique' ELSE 'shared' END AS shared_unique,
                COUNT(*)::bigint AS n
            FROM   structure_class sc
            JOIN   entity_source_count esc USING (entity_id)
            GROUP  BY sc.major_class, shared_unique
        )
        %s
    ", fr007a_projection("major_class"))

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
#' @inheritParams fr007a_entities
#'
#' @return Long-format tibble per \code{\link{fr007a_data_layer}}.
#'
#' @importFrom DBI dbGetQuery
#' @export
fr007a_literature <- function(panel_id = "fig01-overview") {

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
    ", fr007a_projection("major_class"))

    pg_query_panel(panel_id, sql)
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
fr007a_overview <- function(panel_id = "fig01-overview",
                            identifiers = c("major", "authoritative")) {

    identifiers <- match.arg(identifiers)

    facet_fn <- list(
        Entities       = fr007a_entities,
        Associations   = fr007a_associations,
        Interactions   = fr007a_interactions,
        Identifiers    = if (identifiers == "major") fr007a_identifiers_major
                         else fr007a_identifiers_authoritative,
        Structures     = fr007a_structures,
        Literature     = fr007a_literature
    )

    out_list <- list()
    for (facet_name in names(facet_fn)) {
        logger::log_info("FR-007a {facet_name} facet")
        rows <- facet_fn[[facet_name]](panel_id)
        rows$facet <- facet_name
        out_list[[facet_name]] <- rows
    }

    dplyr::bind_rows(out_list)
}
