#' FR-007e specificity × chemical category data layer
#'
#' Cross-tab of \code{structural_specificity} × chemical category
#' bitmaps on dev4 via roaringbitmap intersection. Categories are
#' a mix of \code{chemical_class} (drugs, lipids, food compounds)
#' and \code{metabolic_domain} (amino-acid metabolism, nucleic-acid
#' metabolism, carbohydrates) per the FR-007e spec — both facets
#' are pre-computed on \code{dev4}'s \code{facet_entity_bitmap}.
#'
#' Per-cell counts via \code{rb_and_cardinality}: ~25 ms for the
#' 6 × 6 = 36-cell cross-tab.
#'
#' Routes to \code{dev4} via \code{pg_query_panel("database-content",
#' facet = "panel_e")}.
#'
#' @param panel_id Character: panel identifier. Default
#'     \code{"database-content"}.
#'
#' @return Tibble with columns \code{category}, \code{specificity},
#'     \code{n}.
#'
#' @examples
#' \dontrun{
#' fr007e_specificity_by_category()
#' }
#'
#' @importFrom DBI dbGetQuery
#' @export
fr007e_specificity_by_category <- function(
    panel_id = "database-content"
) {
    sql <- "
        WITH
        spec AS (
            SELECT facet_value AS specificity, entity_bitmap AS bm
            FROM   facet_entity_bitmap
            WHERE  facet_name = 'structural_specificity'
        ),
        cclass AS (
            SELECT facet_value AS class, entity_bitmap AS bm
            FROM   facet_entity_bitmap WHERE facet_name = 'chemical_class'
        ),
        mdom AS (
            SELECT facet_value AS domain, entity_bitmap AS bm
            FROM   facet_entity_bitmap
            WHERE  facet_name = 'metabolic_domain'
        ),
        categories(category, bm) AS (
            SELECT 'lipids',                  bm FROM cclass WHERE class = 'lipid'
            UNION ALL SELECT 'drugs',         bm FROM cclass WHERE class = 'drug'
            UNION ALL SELECT 'food compounds', bm FROM cclass WHERE class = 'food'
            UNION ALL SELECT 'amino-acid metabolism',
                bm FROM mdom WHERE domain = 'amino_acid'
            UNION ALL SELECT 'nucleic-acid metabolism',
                bm FROM mdom WHERE domain = 'nucleotide'
            UNION ALL SELECT 'carbohydrates',
                bm FROM mdom WHERE domain = 'carbohydrate'
        )
        SELECT c.category, s.specificity,
               rb_and_cardinality(c.bm, s.bm)::bigint AS n
        FROM   categories c
        CROSS  JOIN spec s
        ORDER  BY c.category, s.specificity
    "
    pg_query_panel(panel_id, sql, facet = "panel_e")
}
