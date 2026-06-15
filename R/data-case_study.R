#' Load the Shorthouse 2022 differential-metabolomics output
#'
#' Reads one limma sheet (`KRAS_filt_limma` or `EGFR_filt_limma`) from
#' the vendored
#' \code{omnipath_metabo_case1/Results/Differential_Analysis/}
#' ExtendedDataTable xlsx and returns a tibble suitable for the
#' Figure 6 Panel A / B volcanos and Panels E / F top-DEM ChEBI
#' lookup.
#'
#' The xlsx is treated as a vendored input pending the full
#' \code{Differential_Analysis.Rmd} refactor (T058–T061). The loader
#' attaches a \code{"source_path"} attribute and a SHA-256
#' \code{"fingerprint"} attribute for sidecar recording.
#'
#' @param contrast Character: one of \code{"KRAS"} or \code{"EGFR"}.
#' @param xlsx_path Character: path to the xlsx file. Defaults to the
#'     repo-relative vendored location.
#'
#' @return A tibble with at least \code{logFC}, \code{P.Value},
#'     \code{adj.P.Val}, \code{name}, \code{chebi}, \code{hmdb}, plus
#'     a derived \code{neg_log10_p = -log10(P.Value)} column and a
#'     \code{contrast} column carrying the contrast label.
#'
#' @examples
#' \dontrun{
#' kras <- case_study_differential("KRAS")
#' egfr <- case_study_differential("EGFR")
#' }
#'
#' @importFrom rlang abort
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate
#' @importFrom digest digest
#' @export
case_study_differential <- function(
    contrast,
    xlsx_path = file.path(
        "omnipath_metabo_case1",
        "Results",
        "Differential_Analysis",
        "ExtendedDataTable_DifferentialAnalysis_Shorthouse_LungMutations.xlsx"
    )
) {

    # NSE vs. R CMD check workaround
    P.Value <- NULL

    if (!contrast %in% c("KRAS", "EGFR")) {
        rlang::abort(sprintf(
            "Unknown contrast '%s'. Expected one of: KRAS, EGFR.",
            contrast
        ))
    }

    if (!file.exists(xlsx_path)) {
        rlang::abort(sprintf(
            paste0(
                "Differential-analysis xlsx not found at %s. ",
                "This file is a vendored output of the legacy ",
                "Differential_Analysis.Rmd pending the full refactor ",
                "(spec tasks T058 / T061)."
            ),
            xlsx_path
        ))
    }

    rlang::check_installed("readxl")

    sheet <- paste0(contrast, "_filt_limma")
    rows <- tibble::as_tibble(readxl::read_xlsx(xlsx_path, sheet = sheet))

    rows <- dplyr::mutate(
        rows,
        contrast = contrast,
        neg_log10_p = -log10(P.Value)
    )

    fingerprint <- digest::digest(file = xlsx_path, algo = "sha256")
    attr(rows, "source_path") <- xlsx_path
    attr(rows, "fingerprint") <- fingerprint
    attr(rows, "sheet") <- sheet

    rows
}


#' Load the vendored COSMOS PKN CSV fixtures
#'
#' Reads \code{omnipath_metabo_case1/data/pkn_allosteric.csv} and
#' \code{omnipath_metabo_case1/data/pkn_enzyme_metabolite.csv} produced
#' by the legacy \code{01_cosmos_pkn.py} script. These are treated as
#' a fixture pending the T066 refactor that will fetch the COSMOS PKN
#' from the OmniPath Metabo API via \code{omnipath-client}.
#'
#' Each entry is a tibble with the original 11 columns
#' (\code{source}, \code{target}, \code{source_type}, \code{target_type},
#' \code{id_type_a}, \code{id_type_b}, \code{interaction_type},
#' \code{resource}, \code{mor}, \code{locations}, \code{attrs}) and
#' attached \code{"source_path"} + SHA-256 \code{"fingerprint"} for
#' sidecar recording.
#'
#' @param pkn_dir Character: directory holding the two PKN CSVs.
#'     Defaults to the repo-relative vendored location.
#'
#' @return A named list with elements \code{allosteric} and
#'     \code{enzyme_metabolite}.
#'
#' @examples
#' \dontrun{
#' pkn <- case_study_cosmos_pkn()
#' head(pkn$allosteric)
#' head(pkn$enzyme_metabolite)
#' }
#'
#' @importFrom rlang abort
#' @importFrom readr read_csv cols col_character col_double
#' @importFrom digest digest
#' @export
case_study_cosmos_pkn <- function(
    pkn_dir = file.path("omnipath_metabo_case1", "data")
) {

    files <- list(
        allosteric        = file.path(pkn_dir, "pkn_allosteric.csv"),
        enzyme_metabolite = file.path(pkn_dir, "pkn_enzyme_metabolite.csv")
    )

    missing <- files[!vapply(files, file.exists, logical(1L))]
    if (length(missing) > 0L) {
        rlang::abort(sprintf(
            paste0(
                "COSMOS PKN fixture missing: %s. ",
                "Vendored from the legacy 01_cosmos_pkn.py output ",
                "pending the T066 refactor."
            ),
            paste(unlist(missing), collapse = ", ")
        ))
    }

    col_spec <- readr::cols(
        source           = readr::col_character(),
        target           = readr::col_character(),
        source_type      = readr::col_character(),
        target_type      = readr::col_character(),
        id_type_a        = readr::col_character(),
        id_type_b        = readr::col_character(),
        interaction_type = readr::col_character(),
        resource         = readr::col_character(),
        mor              = readr::col_double(),
        locations        = readr::col_character(),
        attrs            = readr::col_character()
    )

    pkn <- lapply(names(files), function(name) {
        path <- files[[name]]
        rows <- readr::read_csv(
            path,
            col_types = col_spec,
            progress = FALSE
        )
        attr(rows, "source_path") <- path
        attr(rows, "fingerprint") <- digest::digest(
            file = path,
            algo = "sha256"
        )
        rows
    })
    names(pkn) <- names(files)

    pkn
}


#' Top-N up-and-down DEMs by limma t-statistic for a contrast
#'
#' Reproduces the legacy
#' \code{omnipath_metabo_case1/scripts/02_connect_dem_pkn.R}
#' \code{top10_dem()} selection: top-N by largest positive t plus
#' top-N by smallest negative t, returned sorted by descending t.
#'
#' @param diff_tibble The output of
#'     \code{\link{case_study_differential}}.
#' @param n Integer: number per direction. Default \code{10L}.
#'
#' @return Tibble subset of \code{diff_tibble}.
#'
#' @importFrom dplyr slice_max slice_min bind_rows arrange desc
#' @importFrom rlang abort
#' @export
case_study_top_dems <- function(diff_tibble, n = 10L) {

    # NSE vs. R CMD check workaround
    t <- NULL

    if (!"t" %in% names(diff_tibble)) {
        rlang::abort(
            "diff_tibble must carry a 't' column (limma t-statistic)"
        )
    }
    if (!is.integer(n) || length(n) != 1L || n < 1L) {
        rlang::abort("`n` must be a positive integer scalar")
    }

    up <- dplyr::slice_max(diff_tibble, t, n = n)
    down <- dplyr::slice_min(diff_tibble, t, n = n)
    dplyr::arrange(dplyr::bind_rows(up, down), dplyr::desc(t))
}


#' Expand the semicolon-separated ChEBI column to a unique ID vector
#'
#' Reproduces the legacy
#' \code{omnipath_metabo_case1/scripts/02_connect_dem_pkn.R}
#' \code{extract_chebi()} helper.
#'
#' @param diff_tibble Tibble with a \code{chebi} column.
#'
#' @return Character vector of unique trimmed ChEBI IDs.
#'
#' @importFrom dplyr filter pull
#' @export
case_study_extract_chebi <- function(diff_tibble) {

    # NSE vs. R CMD check workaround
    chebi <- NULL

    diff_tibble |>
        dplyr::filter(!is.na(chebi)) |>
        dplyr::pull(chebi) |>
        strsplit(";") |>
        unlist() |>
        trimws() |>
        unique()
}
