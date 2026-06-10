#' FR-015a row → bitmap-facet mapping for the resource × record-type table
#'
#' Each row of the resource × record-type checkmark table maps to a
#' (bitmap source, facet name, facet values, deployment) tuple that
#' defines the set of items the row covers. Entries with
#' \code{kind = "blank"} render as empty rows in the table — their
#' underlying vocabulary has not yet landed in the cycle-001 build (see
#' the 2026-06-09 handover notes); per the spec they MUST NOT be
#' synthesised from secondary tables (FR-007d "Other" honesty
#' constraint).
#'
#' Returns a list of named lists with the following fields:
#' \describe{
#'   \item{label}{Character: row label as it appears in the rendered
#'     table (first column).}
#'   \item{kind}{One of \code{"entity"}, \code{"relation"},
#'     \code{"blank"}.}
#'   \item{facet}{Character or NULL: the \code{facet_name} key into
#'     \code{facet_*_bitmap}.}
#'   \item{values}{Character vector or NULL: the
#'     \code{facet_value}(s) to union together for the row's scope
#'     bitmap.}
#'   \item{deployment_facet}{Character or NULL: the registry-facet
#'     name to pass to \code{\link{pg_query_panel}} so dev3 vs dev4
#'     routing follows the registry. \code{NULL} means the panel's
#'     default deployment (dev3 for tab03-record-coverage).}
#'   \item{note}{Character: short rationale when \code{kind == "blank"}.}
#' }
#'
#' @return A list of row definitions in the spec-declared order.
#'
#' @keywords internal
#' @export
record_coverage_row_definitions <- function() {
    list(
        list(label = "Proteins", kind = "entity",
             facet = "entity_type",
             values = c("Gene:MI:0250", "Protein:MI:0326")),
        list(label = "Metabolites", kind = "entity",
             facet = "chemical_class", values = "metabolite"),
        list(label = "Lipids", kind = "entity",
             facet = "chemical_class", values = "lipid"),
        list(label = "Drugs", kind = "entity",
             facet = "chemical_class", values = "drug"),
        list(label = "Food compounds", kind = "entity",
             facet = "chemical_class", values = "food"),
        list(label = "Xenobiotics", kind = "blank",
             note = "chemical_class.xenobiotic deferred (empty class)"),
        list(label = "Diseases", kind = "entity",
             facet = "ontology_id", values = "mondo"),
        list(label = "Phenotypes", kind = "entity",
             facet = "ontology_id", values = "hpo"),
        list(label = "Molecular classes", kind = "entity",
             facet = "entity_type", values = "Cv Term:OM:0012"),
        list(label = "Metabolism (GEMs)", kind = "blank",
             note = "GEM not yet exposed as a distinct facet"),
        list(label = "Allosteric regulation", kind = "blank",
             note = "predicate vocab collapses into Other; see FR-007d"),
        list(label = "Transport", kind = "relation",
             facet = "predicate", values = "transports"),
        list(label = "Ligand-receptor", kind = "blank",
             note = "predicate vocab collapses into Other; see FR-007d"),
        list(label = "Drug-target", kind = "blank",
             note = "predicate vocab collapses into Other; see FR-007d"),
        list(label = "TF-target", kind = "blank",
             note = "predicate vocab collapses into Other; see FR-007d"),
        list(label = "Signaling", kind = "relation",
             facet = "predicate",
             values = c("controls", "positively_regulates",
                        "negatively_regulates")),
        list(label = "Subcellular localization", kind = "blank",
             note = "subcellular component not yet a distinct facet"),
        list(label = "Structures", kind = "entity",
             facet = "structural_specificity",
             values = c("stereospecific", "cis_trans_only",
                        "constitution_only", "variable_constitution",
                        "unknown_constitution"),
             deployment_facet = "structures"),
        list(label = "Literature references", kind = "blank",
             note = "facet_evidence_bitmap not yet derived")
    )
}


#' Per-resource record-coverage counts for one FR-015a row
#'
#' Issues a bitmap-intersection query for one row's scope:
#' \code{rb_and_cardinality(<source bitmap>, <scope bitmap>)} per
#' resource. \code{kind = "entity"} reads \code{facet_entity_bitmap},
#' \code{kind = "relation"} reads \code{facet_relation_bitmap}.
#' Routing follows \code{panel_deployment("tab03-record-coverage",
#' facet = deployment_facet)}: rows with
#' \code{deployment_facet = "structures"} route to dev4, every other
#' row to dev3.
#'
#' @param row_def One element of
#'     \code{\link{record_coverage_row_definitions}} with
#'     \code{kind \%in\% c("entity", "relation")}.
#' @param panel_id Character: panel identifier (default
#'     \code{"tab03-record-coverage"}).
#'
#' @return A tibble with columns \code{resource}, \code{n} (the count
#'     for that resource on that row). Carries the \code{"deployment"}
#'     attribute set by \code{\link{pg_query_panel}}.
#'
#' @keywords internal
#' @export
record_coverage_row_counts <- function(
    row_def,
    panel_id = "tab03-record-coverage"
) {

    bitmap_table <- if (identical(row_def$kind, "entity")) {
        list(table = "facet_entity_bitmap",
             bitmap_col = "entity_bitmap")
    } else {
        list(table = "facet_relation_bitmap",
             bitmap_col = "relation_bitmap")
    }

    values_sql <- paste0(
        "ARRAY[",
        paste(sprintf("'%s'", row_def$values), collapse = ", "),
        "]::text[]"
    )

    sql <- sprintf("
        WITH scope AS (
            SELECT rb_or_agg(%s) AS bm
            FROM   %s
            WHERE  facet_name = '%s'
              AND  facet_value = ANY(%s)
        )
        SELECT
            s.facet_value AS resource,
            rb_and_cardinality(s.%s, scope.bm)::bigint AS n
        FROM   %s s
        CROSS  JOIN scope
        WHERE  s.facet_name = 'source'
        ORDER  BY s.facet_value
    ",
        bitmap_table$bitmap_col,
        bitmap_table$table,
        row_def$facet,
        values_sql,
        bitmap_table$bitmap_col,
        bitmap_table$table
    )

    pg_query_panel(panel_id, sql, facet = row_def$deployment_facet)
}


#' Long tibble of (row, resource, n) for the FR-015a checkmark table
#'
#' Walks \code{\link{record_coverage_row_definitions}}, issues a
#' bitmap-intersection query per non-blank row, and concatenates the
#' results. Blank rows are emitted with \code{n = 0L} for every
#' resource that appears in the entity-bitmap source facet so the
#' wide pivot has consistent shape.
#'
#' @param panel_id Character: panel identifier (default
#'     \code{"tab03-record-coverage"}).
#'
#' @return A long-format tibble with columns \code{row_label},
#'     \code{row_kind}, \code{resource}, \code{n},
#'     \code{deployment}, \code{note} (when row_kind == "blank").
#'     The tibble carries one query-record attribute per non-blank
#'     row, exposed via \code{\link{record_coverage_queries}}.
#'
#' @examples
#' \dontrun{
#' long <- record_coverage_long()
#' }
#'
#' @importFrom dplyr bind_rows mutate
#' @importFrom tibble tibble
#' @export
record_coverage_long <- function(
    panel_id = "tab03-record-coverage"
) {

    rows <- record_coverage_row_definitions()

    # Resource set comes from the dev3 entity-bitmap source facet so
    # blank rows have the same resource columns as populated rows.
    resources_dev3 <- pg_query_panel(
        panel_id,
        "
        SELECT facet_value AS resource
        FROM   facet_entity_bitmap
        WHERE  facet_name = 'source'
        ORDER  BY facet_value
        "
    )

    queries <- list(query_record(resources_dev3))
    all_long <- list()

    for (row_def in rows) {
        if (identical(row_def$kind, "blank")) {
            cells <- tibble::tibble(
                row_label = row_def$label,
                row_kind  = "blank",
                resource  = resources_dev3$resource,
                n         = 0L,
                note      = row_def$note %||% ""
            )
            all_long <- c(all_long, list(cells))
            next
        }

        cell_rows <- record_coverage_row_counts(row_def, panel_id)
        queries <- c(queries, list(query_record(cell_rows)))

        cells <- tibble::tibble(
            row_label = row_def$label,
            row_kind  = row_def$kind,
            resource  = cell_rows$resource,
            n         = as.integer(cell_rows$n),
            note      = ""
        )
        all_long <- c(all_long, list(cells))
    }

    out <- dplyr::bind_rows(all_long)
    attr(out, "queries") <- queries
    out
}


#' Per-row query records accumulated by \code{\link{record_coverage_long}}
#'
#' Convenience accessor for the per-query records the long-format
#' producer accumulates as an attribute. Each record carries the SQL
#' text, the row count, the result hash, and the deployment label
#' (\code{dev3} for the entity-bitmap / relation-bitmap rows,
#' \code{dev4} for the Structures row).
#'
#' @param long_tibble Output of \code{\link{record_coverage_long}}.
#'
#' @return A list of query records (one per non-blank row).
#'
#' @keywords internal
#' @export
record_coverage_queries <- function(long_tibble) {
    attr(long_tibble, "queries") %||% list()
}


#' Wide-format (rows × resources) checkmark tibble for the CSV
#'
#' Pivots \code{\link{record_coverage_long}}'s output to the
#' resource × record-type matrix shape, with cells filled in by
#' whether \code{n \\u2265 threshold}. The first columns are the row
#' label and the row kind (or "blank"); the remaining columns are
#' resources in alphabetical order.
#'
#' @param long_tibble Output of \code{\link{record_coverage_long}}.
#' @param threshold Integer: minimum cell count for a checkmark. Default
#'     \code{1L}. Carried into the artifact sidecar so the caption and
#'     the CSV agree (FR-015a 2026-06-09 note).
#'
#' @return A wide-format tibble with the cell values as \code{"X"}
#'     when checked and empty string otherwise (CSV-friendly).
#'
#' @importFrom dplyr mutate select
#' @importFrom tidyr pivot_wider
#' @importFrom rlang .data
#' @export
record_coverage_wide <- function(long_tibble, threshold = 1L) {

    long_tibble %>%
        dplyr::mutate(
            cell = ifelse(.data$n >= threshold, "X", "")
        ) %>%
        dplyr::select("row_label", "resource", "cell") %>%
        tidyr::pivot_wider(
            names_from  = "resource",
            values_from = "cell",
            values_fill = ""
        )
}


#' LaTeX body for the FR-015a resource × record-type checkmark table
#'
#' Hand-built \code{tabular} with \code{\\rotatebox{45}{...}} column
#' labels — gt's LaTeX backend escapes backslashes inside column
#' labels, so a custom render is the simpler path. Cells are
#' \code{\\checkmark} for true and blank otherwise (FR-015a). The
#' first column is the record-type label; subsequent columns are
#' resources in the order they appear in the wide tibble.
#'
#' @param wide_tibble Output of \code{\link{record_coverage_wide}}.
#'
#' @return Character scalar: a LaTeX body renderable inside the
#'     standalone wrapper from \code{\link{tables_save_latex_pdf_csv}}.
#'
#' @keywords internal
#' @export
record_coverage_latex <- function(wide_tibble) {

    resources <- setdiff(names(wide_tibble), "row_label")
    n_res <- length(resources)

    col_spec <- sprintf("l|%s", paste(rep("c", n_res), collapse = ""))

    header_cells <- vapply(resources, function(r) {
        sprintf(
            "\\rotatebox{60}{\\sffamily\\fontsize{8pt}{10pt}\\selectfont %s}",
            latex_escape_text(r)
        )
    }, character(1L))
    header_row <- paste(c("", header_cells), collapse = " & ")

    body_rows <- vapply(seq_len(nrow(wide_tibble)), function(i) {
        row <- wide_tibble[i, ]
        cells <- vapply(resources, function(r) {
            if (identical(row[[r]], "X")) "\\checkmark" else ""
        }, character(1L))
        sprintf(
            "%s & %s \\\\",
            sprintf(
                "\\sffamily\\fontsize{8pt}{10pt}\\selectfont %s",
                latex_escape_text(row$row_label)
            ),
            paste(cells, collapse = " & ")
        )
    }, character(1L))

    paste(c(
        "\\begingroup",
        "\\setlength{\\tabcolsep}{3pt}",
        sprintf("\\begin{tabular}{%s}", col_spec),
        "\\toprule",
        paste(header_row, "\\\\[1ex]"),
        "\\midrule",
        body_rows,
        "\\bottomrule",
        "\\end{tabular}",
        "\\endgroup"
    ), collapse = "\n")
}


#' Conservatively escape table-cell text for LaTeX
#'
#' Escapes the small set of LaTeX specials we expect in resource
#' identifiers and row labels (underscore, ampersand, percent, hash,
#' dollar). Leaves other characters as-is; the wrapper preamble
#' loads \code{fontspec}-friendly settings so UTF-8 passes through.
#'
#' @param text Character vector.
#'
#' @return Character vector with specials escaped.
#'
#' @keywords internal
#' @noRd
latex_escape_text <- function(text) {
    text <- gsub("\\\\", "\\\\textbackslash{}", text)
    text <- gsub("&",  "\\\\&", text)
    text <- gsub("%",  "\\\\%", text)
    text <- gsub("\\$", "\\\\$", text)
    text <- gsub("#",  "\\\\#", text)
    text <- gsub("_",  "\\\\_", text)
    text
}
