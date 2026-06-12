#' Per-resource identifier-type counts for the FR-014 Methods table
#'
#' Queries \code{identifier_evidence} joined to
#' \code{vocab_identifier_type} and \code{entity_evidence_identifier}
#' on \code{dev3} (default for \code{tab01-id-resolving}). The
#' identifier types are bucketed into a hand-picked set of
#' manuscript-relevant categories (chemical structure, chemical names,
#' resource-specific chemical ids, gene/protein ids) rather than the
#' raw top-N — the original top-N output mixed populated chemical ids
#' with `Synonym` / `Name` columns of unclear distinction and left no
#' room for the gene-side ids the manuscript needs. The long tail
#' collapses into a single \code{Misc} bucket.
#'
#' Performance note: the underlying join is the same shape as the
#' FR-007a Identifiers facet of Figure 1, which has no
#' \code{facet_identifier_bitmap} yet (see
#' \code{saezverse/human/plans/omnipath-improvements-2026-06-identifier-source-count.md}
#' — proposes an \code{identifier_source_count} derived table).
#' Expect ~30 minutes on dev3 against the cycle-001 build until that
#' derived table lands.
#'
#' @param panel_id Character: panel identifier (default
#'     \code{"tab01-id-resolving"} — routes to \code{dev3}).
#'
#' @return A long-format tibble with columns \code{resource},
#'     \code{id_type} (the hand-picked bucket label),
#'     \code{n_identifiers}. Carries the \code{"deployment"}
#'     attribute set by \code{\link{pg_query_panel}}.
#'
#' @examples
#' \dontrun{
#' long <- tbl_id_resolving_counts()
#' attr(long, "deployment")   # "dev3"
#' }
#'
#' @importFrom DBI dbGetQuery
#' @export
tbl_id_resolving_counts <- function(panel_id = "tab01-id-resolving") {

    buckets <- tbl_id_resolving_buckets()
    case_when <- vapply(seq_along(buckets), function(i) {
        b <- buckets[[i]]
        sql_vals <- paste(sprintf("'%s'", b$values), collapse = ", ")
        sprintf("WHEN vit.name IN (%s) THEN '%s'", sql_vals, b$label)
    }, character(1L))
    case_sql <- paste(c(
        "CASE",
        paste0("    ", case_when),
        "    ELSE 'Misc'",
        "END AS id_type"
    ), collapse = "\n")

    excluded <- paste(
        sprintf("'%s'", tbl_id_resolving_excluded_types()),
        collapse = ", "
    )

    sql <- sprintf("
        WITH id_classified AS (
            SELECT DISTINCT
                ie.identifier_id,
                %s
            FROM   identifier_evidence ie
            JOIN   vocab_identifier_type vit
                   ON vit.identifier_type_id = ie.identifier_type_id
            WHERE  vit.name NOT IN (%s)
        )
        SELECT
            ds.name AS resource,
            ic.id_type,
            COUNT(DISTINCT ic.identifier_id)::bigint AS n_identifiers
        FROM   id_classified ic
        JOIN   entity_evidence_identifier eei USING (identifier_id)
        JOIN   data_source ds ON ds.source_id = eei.source_id
        WHERE  ds.name <> 'omnipath_ontology'
        GROUP  BY ds.name, ic.id_type
        ORDER  BY ds.name, ic.id_type
    ", case_sql, excluded)

    pg_query_panel(panel_id, sql)
}


#' Manuscript-relevant identifier-type buckets for the FR-014 table
#'
#' Each entry maps a manuscript-facing column label to the set of raw
#' \code{vocab_identifier_type.name} values that contribute to it.
#' Buckets are stable, ordered, and chosen for the Methods-table
#' audience — they cover the major chemical structural ids
#' (\code{SMILES}, \code{InChIKey}), the chemical names (\code{Name}
#' combines \code{Name} + \code{Synonym}), the most-populated
#' resource-specific chemical ids, and the gene / protein ids the
#' manuscript needs (\code{UniProt}, \code{Entrez}, \code{Ensembl},
#' \code{Gene Name}, \code{HGNC}). Any id type not matching is rolled
#' into \code{Misc} by the consuming query.
#'
#' @return A list of named lists, each with \code{label} (column
#'     header) and \code{values} (raw \code{vit.name} strings).
#'
#' @keywords internal
#' @export
tbl_id_resolving_buckets <- function() {
    list(
        list(label = "SMILES",
             values = "Smiles:MI:0239"),
        list(label = "InChIKey",
             values = "Standard Inchi Key:MI:1101"),
        list(label = "Name",
             values = c(
                 "Name:OM:0202",
                 "Synonym:OM:0203",
                 "Iupac Name:OM:0210",
                 "Iupac Traditional Name:OM:0211"
             )),
        list(label = "ChEMBL",
             values = c("Chembl Compound:MI:0967",
                        "Chembl Target:MI:1348")),
        list(label = "PubChem",
             values = c("Pubchem Compound:OM:0002",
                        "Pubchem:MI:0730")),
        list(label = "ChEBI",
             values = "Chebi:MI:0474"),
        list(label = "HMDB",
             values = "Hmdb:OM:0004"),
        list(label = "SwissLipids",
             values = "Swisslipids:OM:0009"),
        list(label = "LIPID MAPS",
             values = "Lipidmaps:OM:0003"),
        list(label = "MetaNetX",
             values = "Metanetx:OM:0005"),
        list(label = "UniProt",
             values = c("Uniprot:MI:1097", "Uniprot Trembl:MI:1099",
                        "Uniprot Entry Name:OM:0221")),
        list(label = "Entrez",
             values = "Entrez:MI:0477"),
        list(label = "Gene Symbol",
             values = c("Gene Name Primary:OM:0200",
                        "Gene Name Synonym:OM:0201"))
    )
}


#' Identifier-type vocab values to exclude from the FR-014 query
#'
#' Internal scaffold keys the cycle-001 build emits to glue
#' identifier_evidence rows together (\code{omnipath:unresolved_entity_key},
#' \code{omnipath:reaction_member_hash}) inflate the \code{Misc}
#' bucket if not filtered. They are pipeline-internal and have no
#' manuscript meaning.
#'
#' @return Character vector of vocab_identifier_type.name values.
#'
#' @keywords internal
#' @export
tbl_id_resolving_excluded_types <- function() {
    c(
        "omnipath:unresolved_entity_key",
        "omnipath:reaction_member_hash"
    )
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
#' @param resource_labels Named character: optional resource_id →
#'     display label map (as returned by
#'     \code{\link{resources_label_map}}). When supplied, the
#'     \code{resource} column is replaced with the resource_short
#'     display label so the Methods table shows e.g. "ChEMBL" rather
#'     than the lowercase \code{chembl} slug.
#'
#' @importFrom dplyr arrange desc group_by summarise ungroup mutate select
#' @importFrom dplyr left_join across all_of
#' @importFrom tidyr pivot_wider replace_na
#' @importFrom rlang .data
#' @export
tbl_id_resolving_wide <- function(long_tibble,
                                  resource_labels = NULL) {

    if (nrow(long_tibble) == 0L) {
        return(tibble::tibble(resource = character(), Total = integer()))
    }

    # Bucket order follows tbl_id_resolving_buckets() declaration —
    # SMILES, InChIKey, Name, then chemical IDs by resource then
    # gene/protein IDs; Misc pinned last.
    bucket_labels <- vapply(
        tbl_id_resolving_buckets(),
        function(b) b$label,
        character(1L)
    )
    type_order <- intersect(c(bucket_labels, "Misc"), unique(long_tibble$id_type))

    if (!is.null(resource_labels)) {
        long_tibble <- dplyr::mutate(
            long_tibble,
            resource = unname(
                ifelse(
                    is.na(resource_labels[.data$resource]),
                    .data$resource,
                    resource_labels[.data$resource]
                )
            )
        )
    }

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
