#' Per-resource identifier-type counts for the FR-014 Methods table
#'
#' Queries \code{identifier_evidence} joined to
#' \code{vocab_identifier_type} and \code{entity_evidence_identifier}
#' on \code{dev3} (default for \code{tab01-id-resolving}). The top-N
#' identifier types by global frequency are kept as named columns;
#' the long tail collapses into a single \code{Misc} bucket so the
#' table stays inside the page width.
#'
#' Performance note: the underlying join is the same shape as the
#' FR-007a Identifiers facet of Figure 1, which has no
#' \code{facet_identifier_bitmap} yet (see
#' \code{saezverse/human/plans/omnipath-improvements-2026-06-identifier-source-count.md}
#' — proposes an \code{identifier_source_count} derived table).
#' Expect ~60-90 s on dev3 against the cycle-001 build until that
#' derived table lands.
#'
#' @param panel_id Character: panel identifier (default
#'     \code{"tab01-id-resolving"} — routes to \code{dev3}).
#' @param top_n Integer: number of identifier types to keep before
#'     collapsing the tail into \code{Misc}. Default \code{10L}.
#'
#' @return A long-format tibble with columns \code{resource},
#'     \code{id_type}, \code{n_identifiers}. Carries the
#'     \code{"deployment"} attribute set by \code{\link{pg_query_panel}}.
#'
#' @examples
#' \dontrun{
#' long <- tbl_id_resolving_counts()
#' attr(long, "deployment")   # "dev3"
#' }
#'
#' @importFrom DBI dbGetQuery
#' @export
tbl_id_resolving_counts <- function(panel_id = "tab01-id-resolving",
                                    top_n = 10L) {

    sql <- sprintf("
        WITH ranked_types AS (
            SELECT vit.name AS id_type, COUNT(*) AS n,
                   ROW_NUMBER() OVER (ORDER BY COUNT(*) DESC) AS rk
            FROM   identifier_evidence ie
            JOIN   vocab_identifier_type vit
                   ON vit.identifier_type_id = ie.identifier_type_id
            GROUP  BY vit.name
        ),
        id_classified AS (
            SELECT DISTINCT
                ie.identifier_id,
                CASE WHEN rt.rk <= %d
                     THEN rt.id_type ELSE 'Misc' END AS id_type
            FROM   identifier_evidence ie
            JOIN   vocab_identifier_type vit
                   ON vit.identifier_type_id = ie.identifier_type_id
            JOIN   ranked_types rt ON rt.id_type = vit.name
        )
        SELECT
            ds.name AS resource,
            ic.id_type,
            COUNT(DISTINCT ic.identifier_id)::bigint AS n_identifiers
        FROM   id_classified ic
        JOIN   entity_evidence_identifier eei USING (identifier_id)
        JOIN   data_source ds ON ds.source_id = eei.source_id
        GROUP  BY ds.name, ic.id_type
        ORDER  BY ds.name, ic.id_type
    ", as.integer(top_n))

    pg_query_panel(panel_id, sql)
}


#' Pivot the long id-resolving tibble to a wide resource × id-type matrix
#'
#' Sorts the id-type columns by total identifier count descending so
#' the most populated columns sit at the left of the table; the
#' \code{Misc} column (if present) is always pinned to the right end
#' just before the row \code{Total}. Resources are sorted by their row
#' total, descending. Empty cells are filled with \code{0L}.
#'
#' @param long_tibble Tibble from \code{\link{tbl_id_resolving_counts}}.
#'
#' @return A wide-format tibble with one row per resource and one
#'     column per identifier type, plus a final \code{Total} column.
#'     The first column is \code{resource}.
#'
#' @examples
#' \dontrun{
#' wide <- tbl_id_resolving_wide(tbl_id_resolving_counts())
#' }
#'
#' @importFrom dplyr arrange desc group_by summarise ungroup mutate select
#' @importFrom dplyr left_join across all_of
#' @importFrom tidyr pivot_wider replace_na
#' @importFrom stringr str_replace_all
#' @importFrom rlang .data
#' @export
tbl_id_resolving_wide <- function(long_tibble) {

    if (nrow(long_tibble) == 0L) {
        return(tibble::tibble(resource = character(), Total = integer()))
    }

    # Strip the trailing ":NS:NNNN" CV suffix (e.g.
    # "Standard Inchi Key:MI:1101" → "Standard Inchi Key") so the
    # column headers stay readable in the Methods table without
    # losing the SQL-level traceability — the raw vocab name is
    # still in the sidecar via the recorded query result hash.
    long_tibble <- dplyr::mutate(
        long_tibble,
        id_type = stringr::str_replace_all(
            .data$id_type, ":[A-Za-z]{2}:\\d+$", ""
        )
    )

    type_totals <- long_tibble %>%
        dplyr::group_by(.data$id_type) %>%
        dplyr::summarise(
            type_total = sum(.data$n_identifiers),
            .groups    = "drop"
        )

    misc <- type_totals$id_type == "Misc"
    type_order <- c(
        type_totals$id_type[!misc][
            order(type_totals$type_total[!misc], decreasing = TRUE)
        ],
        type_totals$id_type[misc]
    )

    wide <- long_tibble %>%
        tidyr::pivot_wider(
            id_cols     = "resource",
            names_from  = "id_type",
            values_from = "n_identifiers",
            values_fill = 0L
        ) %>%
        dplyr::mutate(
            Total = rowSums(
                dplyr::across(dplyr::all_of(type_order)),
                na.rm = TRUE
            )
        ) %>%
        dplyr::arrange(dplyr::desc(.data$Total)) %>%
        dplyr::select(
            "resource",
            dplyr::all_of(type_order),
            "Total"
        )

    wide
}


#' Render the FR-014 id-resolving table as a styled \code{gt} object
#'
#' Applies \code{\link{gt_metabo_style}} for FR-016 / FR-017a
#' typography, formats counts with thousand-separators via
#' \code{gt::fmt_number}, bolds the row total column, and emits a
#' grand-total summary row below the body via
#' \code{gt::grand_summary_rows}.
#'
#' @param wide_tibble Wide-format tibble from
#'     \code{\link{tbl_id_resolving_wide}}.
#'
#' @return A \code{gt_tbl} styled per FR-016 / FR-017a.
#'
#' @examples
#' \dontrun{
#' wide <- tbl_id_resolving_wide(tbl_id_resolving_counts())
#' tbl_id_resolving_gt(wide)
#' }
#'
#' @importFrom gt gt cols_label fmt_number tab_header tab_style
#' @importFrom gt cells_body cells_column_labels cells_grand_summary
#' @importFrom gt cell_text grand_summary_rows
#' @importFrom rlang check_installed
#' @export
tbl_id_resolving_gt <- function(wide_tibble) {

    rlang::check_installed("gt")

    id_cols <- setdiff(names(wide_tibble), "resource")

    wide_tibble %>%
        gt::gt(rowname_col = "resource") %>%
        gt::tab_header(
            title = "Identifier resolving across the integrated resources"
        ) %>%
        gt::fmt_number(
            columns  = dplyr::all_of(id_cols),
            decimals = 0,
            sep_mark = ","
        ) %>%
        gt::tab_style(
            style    = gt::cell_text(weight = "bold"),
            locations = gt::cells_body(columns = "Total")
        ) %>%
        gt::grand_summary_rows(
            columns = dplyr::all_of(id_cols),
            fns     = list(Total = ~ sum(.x, na.rm = TRUE)),
            fmt     = list(~ gt::fmt_number(., decimals = 0, sep_mark = ","))
        ) %>%
        gt_metabo_style()
}
