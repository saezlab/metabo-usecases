#' RaMP InChIKey-conflict summary for the FR-015 Methods table
#'
#' Per-reason aggregation of \code{metabo_ramp_inchikey_conflict} on
#' \code{dev4} (the structural-specificity / RaMP-conflict build per
#' FR-030). Each row of the source table is a conflict pair
#' \code{(ramp_id, inchikey_a, inchikey_b, conflict_reason)} produced
#' in-database by the cycle-001 RDKit Postgres cartridge; this
#' renderer reports per-reason distinct RaMP-id and conflict-pair
#' counts plus two example RaMP ids.
#'
#' \code{pg_query_panel} routes the query to \code{dev4} automatically
#' via the \code{ramp-comparison} registry default.
#'
#' @param panel_id Character: panel identifier (default
#'     \code{"ramp-comparison"}).
#' @param example_count Integer: how many example RaMP ids to include
#'     per reason. Default \code{2L}.
#'
#' @return A tibble with columns \code{conflict_reason},
#'     \code{n_ramp_ids}, \code{n_conflict_pairs}, \code{share_pct},
#'     \code{examples}. Carries the \code{"deployment"} attribute set
#'     by \code{\link{pg_query_panel}}.
#'
#' @examples
#' \dontrun{
#' summary <- tbl_ramp_comparison_summary()
#' attr(summary, "deployment")   # "dev4"
#' }
#'
#' @importFrom DBI dbGetQuery
#' @export
tbl_ramp_comparison_summary <- function(
    panel_id      = "ramp-comparison",
    example_count = 2L
) {

    sql <- sprintf("
        WITH per_reason AS (
            SELECT
                conflict_reason,
                COUNT(*)::bigint AS n_conflict_pairs,
                COUNT(DISTINCT ramp_id)::bigint AS n_ramp_ids,
                (
                    SELECT array_agg(r ORDER BY r)
                    FROM (
                        SELECT DISTINCT ramp_id AS r
                        FROM   metabo_ramp_inchikey_conflict mc
                        WHERE  mc.conflict_reason = m.conflict_reason
                        ORDER  BY ramp_id
                        LIMIT  %d
                    ) s
                ) AS examples
            FROM metabo_ramp_inchikey_conflict m
            GROUP BY conflict_reason
        ),
        grand_totals AS (
            SELECT
                SUM(n_ramp_ids)::bigint AS total_ramp_ids
            FROM per_reason
        )
        SELECT
            pr.conflict_reason,
            pr.n_ramp_ids,
            pr.n_conflict_pairs,
            ROUND(
                100.0 * pr.n_ramp_ids / gt.total_ramp_ids,
                1
            )::numeric AS share_pct,
            array_to_string(pr.examples, ', ') AS examples
        FROM per_reason pr
        CROSS JOIN grand_totals gt
        ORDER BY pr.n_conflict_pairs DESC, pr.conflict_reason
    ", as.integer(example_count))

    rows <- pg_query_panel(panel_id, sql)

    # Replace the cycle-001 RDKit classifier's terse one-word reason
    # codes with the manuscript-facing labels.
    label_map <- tbl_ramp_comparison_reason_labels()
    rows$conflict_reason <- ifelse(
        is.na(label_map[rows$conflict_reason]),
        rows$conflict_reason,
        unname(label_map[rows$conflict_reason])
    )
    rows
}


#' RaMP conflict-reason code → manuscript label map
#'
#' The cycle-001 RDKit cartridge tags each conflict pair with a terse
#' one-word reason code (\code{stereo}, \code{specificity},
#' \code{tautomer}, \code{similar}, \code{unrelated}). These are
#' opaque outside the build context; the Methods table needs to spell
#' them out for the reader.
#'
#' @return Named character vector mapping the raw reason code to the
#'     manuscript label.
#'
#' @keywords internal
#' @export
tbl_ramp_comparison_reason_labels <- function() {
    c(
        stereo      = "Stereochemistry only",
        specificity = "Different specificity levels",
        tautomer    = "Tautomeric variants",
        similar     = "Similar structures",
        unrelated   = "Unrelated structures"
    )
}


#' Render the FR-015 RaMP comparison table as a styled \code{gt} object
#'
#' One row per conflict reason plus a grand-total row produced via
#' \code{gt::grand_summary_rows}. The example RaMP ids are rendered
#' in monospace via \code{gt::cell_text(font = "monospace")} so the
#' identifiers stay readable; the share column is formatted as a
#' percentage to one decimal.
#'
#' @param summary_tibble Tibble from
#'     \code{\link{tbl_ramp_comparison_summary}}.
#'
#' @return A \code{gt_tbl} styled per FR-016 / FR-017a.
#'
#' @examples
#' \dontrun{
#' tbl_ramp_comparison_gt(tbl_ramp_comparison_summary())
#' }
#'
#' @importFrom gt gt cols_label fmt_number fmt_percent tab_header
#' @importFrom gt tab_style cells_body cell_text grand_summary_rows
#' @importFrom rlang check_installed
#' @export
tbl_ramp_comparison_gt <- function(summary_tibble) {

    rlang::check_installed("gt")

    summary_tibble %>%
        gt::gt(rowname_col = "conflict_reason") %>%
        gt::cols_label(
            n_ramp_ids       = "Distinct RaMP ids",
            n_conflict_pairs = "Conflict pairs",
            share_pct        = "Share (%)",
            examples         = "Example RaMP ids"
        ) %>%
        gt::fmt_number(
            columns  = c("n_ramp_ids", "n_conflict_pairs"),
            decimals = 0,
            sep_mark = ","
        ) %>%
        gt::fmt_number(
            columns  = "share_pct",
            decimals = 1
        ) %>%
        gt::tab_style(
            style    = gt::cell_text(font = "monospace"),
            locations = gt::cells_body(columns = "examples")
        ) %>%
        gt::grand_summary_rows(
            columns = c("n_ramp_ids", "n_conflict_pairs"),
            fns     = list(Total = ~ sum(.x, na.rm = TRUE)),
            fmt     = list(~ gt::fmt_number(., decimals = 0, sep_mark = ",")),
            missing_text = ""
        ) %>%
        gt::tab_style(
            style    = list(
                gt::cell_text(weight = "bold"),
                gt::cell_borders(
                    sides = "top",
                    color = "black",
                    weight = gt::px(1)
                )
            ),
            locations = gt::cells_grand_summary()
        ) %>%
        gt_metabo_style()
}
