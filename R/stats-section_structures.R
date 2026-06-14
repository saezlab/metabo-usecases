#' Canonical six structural-specificity levels
#'
#' The names MUST match the values stored in
#' \code{metabo_entity_structural_specificity.specificity_level} on
#' the integrated dev5 build (post-2026-06-14 dev5 integrated-build
#' promotion). The metric-name suffix is the snake_cased level.
#'
#' @return Character vector.
#'
#' @keywords internal
#' @noRd
structural_specificity_levels <- function() {
    c(
        "stereospecific",
        "cis_trans_only",
        "constitution_only",
        "variable_constitution",
        "unknown_constitution",
        "no_structure"
    )
}


#' Section 4 — Structures (FR-043d)
#'
#' Reports total structure records plus the six per-specificity-level
#' counts plus distinct constitutional skeletons (InChIKey first
#' block) plus distinct full InChIKey values plus the structures
#' with complete stereo (subset of the stereospecific level) plus the
#' RaMP-conflict count. The InChIKey column name is read from
#' \code{runtime.structures.inchikey_column} so a build-side rename
#' is a one-line config change.
#'
#' Routes to \code{dev5} via \code{\link{pg_query_panel}}
#' (post-2026-06-14 dev5 integrated-build promotion — the integrated
#' build carries \code{metabo_entity_structural_specificity} and
#' \code{metabo_ramp_inchikey_conflict}, which were
#' \code{dev4}-only prior to the promotion).
#'
#' @param panel_id Character.
#' @param runtime Named list: from \code{\link{digest_runtime}};
#'     \code{NULL} → load fresh.
#'
#' @return Tibble of DigestMetric rows.
#'
#' @importFrom logger log_info
#' @importFrom tibble tibble
#' @export
section_structures <- function(
    panel_id = "fig01-architecture",
    runtime = NULL
) {

    rt <- runtime %||% digest_runtime()
    inchikey_col <- rt$structures$inchikey_column %||% "standard_inchikey"
    facet <- "panel_a_stats_structures"

    levels_sql <- structures_levels_sql()
    levels_rows <- pg_query_panel(panel_id, levels_sql, facet = facet)

    inchikey_sql <- structures_inchikey_sql(inchikey_col)
    inchikey_rows <- pg_query_panel(panel_id, inchikey_sql, facet = facet)

    ramp_sql <- ramp_conflict_sql()
    ramp_rows <- pg_query_panel(panel_id, ramp_sql, facet = facet)

    levels <- structural_specificity_levels()
    # Map the integer level_id 1..6 back to the canonical name.
    levels_lookup <- as.list(
        setNames(
            as.integer(levels_rows$n_structures),
            levels[as.integer(levels_rows$level_id)]
        )
    )
    n_total <- as.integer(sum(unlist(levels_lookup)))

    level_metrics <- purrr::map(levels, function(level) {
        value <- as.integer(levels_lookup[[level]] %||% 0L)
        tibble::tibble(
            section_id       = 4L,
            metric_name      = paste0("n_structures_", level),
            value            = value,
            definition_label = sprintf(
                "metabo_entity_structural_specificity.specificity_level:%s",
                level
            ),
            state            = if (value > 0L) "populated" else "empty",
            deployment       = attr(levels_rows, "deployment"),
            sql_hash         = sql_sha256(levels_sql)
        )
    })

    metrics <- dplyr::bind_rows(
        tibble::tibble(
            section_id       = 4L,
            metric_name      = "n_structures_total",
            value            = n_total,
            definition_label = "metabo_entity_structural_specificity:total",
            state            = if (n_total > 0L) "populated" else "empty",
            deployment       = attr(levels_rows, "deployment"),
            sql_hash         = sql_sha256(levels_sql)
        ),
        dplyr::bind_rows(level_metrics),
        tibble::tibble(
            section_id       = 4L,
            metric_name      = "n_constitutional_skeletons",
            value            = as.integer(
                inchikey_rows$n_skeletons %||% 0L
            ),
            definition_label = sprintf(
                "inchikey:first_block(%s)",
                inchikey_col
            ),
            state            = if (
                as.integer(inchikey_rows$n_skeletons %||% 0L) > 0L
            ) {
                "populated"
            } else {
                "empty"
            },
            deployment       = attr(inchikey_rows, "deployment"),
            sql_hash         = sql_sha256(inchikey_sql)
        ),
        tibble::tibble(
            section_id       = 4L,
            metric_name      = "n_full_inchikey",
            value            = as.integer(
                inchikey_rows$n_full_inchikey %||% 0L
            ),
            definition_label = sprintf(
                "inchikey:distinct(%s)",
                inchikey_col
            ),
            state            = if (
                as.integer(inchikey_rows$n_full_inchikey %||% 0L) > 0L
            ) {
                "populated"
            } else {
                "empty"
            },
            deployment       = attr(inchikey_rows, "deployment"),
            sql_hash         = sql_sha256(inchikey_sql)
        ),
        tibble::tibble(
            section_id       = 4L,
            metric_name      = "n_ramp_conflicts",
            value            = as.integer(
                ramp_rows$n_ramp_conflicts %||% 0L
            ),
            definition_label = "metabo_ramp_inchikey_conflict:total",
            state            = if (
                as.integer(ramp_rows$n_ramp_conflicts %||% 0L) > 0L
            ) {
                "populated"
            } else {
                "empty"
            },
            deployment       = attr(ramp_rows, "deployment"),
            sql_hash         = sql_sha256(ramp_sql)
        )
    )

    queries <- list(
        levels = list(
            metric_names = c(
                "n_structures_total",
                paste0("n_structures_", levels)
            ),
            sql          = levels_sql,
            sql_hash     = sql_sha256(levels_sql),
            deployment   = attr(levels_rows, "deployment"),
            result_hash  = attr(levels_rows, "result_hash"),
            row_count    = nrow(levels_rows)
        ),
        inchikey = list(
            metric_names = c(
                "n_constitutional_skeletons",
                "n_full_inchikey"
            ),
            sql          = inchikey_sql,
            sql_hash     = sql_sha256(inchikey_sql),
            deployment   = attr(inchikey_rows, "deployment"),
            result_hash  = attr(inchikey_rows, "result_hash"),
            row_count    = nrow(inchikey_rows)
        ),
        ramp = list(
            metric_names = "n_ramp_conflicts",
            sql          = ramp_sql,
            sql_hash     = sql_sha256(ramp_sql),
            deployment   = attr(ramp_rows, "deployment"),
            result_hash  = attr(ramp_rows, "result_hash"),
            row_count    = nrow(ramp_rows)
        )
    )

    attr(metrics, "queries") <- queries

    logger::log_info(
        "section_structures: n_total={n_total}, ",
        "skeletons={as.integer(inchikey_rows$n_skeletons %||% 0L)}, ",
        "full_inchikey={as.integer(inchikey_rows$n_full_inchikey %||% 0L)}, ",
        "ramp_conflicts={as.integer(ramp_rows$n_ramp_conflicts %||% 0L)}"
    )

    metrics
}


#' SQL — per-specificity-level structure counts
#'
#' \code{metabo_entity_structural_specificity} stores the level as
#' an integer \code{structural_specificity_id}; on dev5 there is no
#' vocab table mapping IDs to names so the digest hard-codes the
#' canonical six-level mapping (FR-043d /
#' \code{structural_specificity_levels()}). IDs 1..6 map to
#' \code{stereospecific, cis_trans_only, constitution_only,
#' variable_constitution, unknown_constitution, no_structure}
#' respectively per the omnipath-metabo build's enum order.
#'
#' @return Character.
#'
#' @keywords internal
#' @noRd
structures_levels_sql <- function() {
    paste(
        "SELECT structural_specificity_id AS level_id,",
        "       COUNT(*)::bigint AS n_structures",
        "  FROM metabo_entity_structural_specificity",
        " GROUP BY structural_specificity_id",
        sep = "\n"
    )
}


#' SQL — distinct constitutional skeletons + distinct full InChIKey
#'
#' Splits the InChIKey on its `-` separator (the standard
#' first-block / second-block layout: 14-char skeleton hash,
#' 10-char protonation+stereo hash, 1-char version). The skeleton
#' count counts distinct 14-character prefixes.
#'
#' @param inchikey_col Character: column name carrying the InChIKey
#'     in \code{metabo_entity_structural_specificity}.
#'
#' @return Character.
#'
#' @keywords internal
#' @noRd
structures_inchikey_sql <- function(inchikey_col) {
    paste(
        sprintf(
            paste0(
                "SELECT",
                "    COUNT(DISTINCT SUBSTRING(%s FROM 1 FOR 14))::bigint",
                "        AS n_skeletons,",
                "    COUNT(DISTINCT %s)::bigint AS n_full_inchikey"
            ),
            inchikey_col,
            inchikey_col
        ),
        "  FROM metabo_entity_structural_specificity",
        sprintf(" WHERE %s IS NOT NULL", inchikey_col),
        sep = "\n"
    )
}


#' SQL — RaMP-conflict count (distinct RaMP IDs)
#'
#' \code{metabo_ramp_inchikey_conflict} on dev5 has one row per
#' conflict pair (\code{ramp_id, inchikey_a, inchikey_b,
#' conflict_reason}); the digest reports DISTINCT RaMP IDs since
#' that is the user-facing count ("how many RaMP ids carry a
#' conflict").
#'
#' @return Character.
#'
#' @keywords internal
#' @noRd
ramp_conflict_sql <- function() {
    paste(
        "SELECT COUNT(DISTINCT ramp_id)::bigint AS n_ramp_conflicts",
        "  FROM metabo_ramp_inchikey_conflict",
        sep = "\n"
    )
}
