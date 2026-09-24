#' Load the azelate Panel C source table
#'
#' Reads
#' \code{data/derived/azelate/csvs/interaction_types_by_source.csv}
#' — the resolved-only interaction-type × source × relation-count
#' table produced by
#' \code{analyses/azelate/azelaic_acid_query_tables.py}.
#' Drives Figure 6 Panel C (interaction-type composition for the
#' resolved \emph{Azelate} entity, stacked by upstream source).
#'
#' @param csv_path Character: path to \code{interaction_types_by_source.csv}.
#'     Defaults to the repo-relative vendored location.
#'
#' @return A tibble with columns \code{interaction_type},
#'     \code{source}, \code{query_entity_label},
#'     \code{entity_resolution_status}, \code{relation_count}.
#'     Carries \code{source_path} + SHA-256 \code{fingerprint}
#'     attributes for sidecar recording.
#'
#' @examples
#' \dontrun{
#' azelate_panel_c_data()
#' }
#'
#' @importFrom readr read_csv cols col_character col_integer
#' @importFrom rlang abort
#' @importFrom digest digest
#' @export
azelate_panel_c_data <- function(
    csv_path = file.path(
        "data",
        "derived",
        "azelate",
        "csvs",
        "interaction_types_by_source.csv"
    )
) {

    if (!file.exists(csv_path)) {
        rlang::abort(sprintf(
            paste0(
                "Azelate Panel C source CSV not found at %s. ",
                "Regenerate via azelaic_acid_query_tables.py per ",
                "analyses/azelate/README.md."
            ),
            csv_path
        ))
    }

    rows <- readr::read_csv(
        csv_path,
        col_types = readr::cols(
            interaction_type         = readr::col_character(),
            source                   = readr::col_character(),
            query_entity_label       = readr::col_character(),
            entity_resolution_status = readr::col_character(),
            relation_count           = readr::col_integer()
        ),
        progress = FALSE
    )

    attr(rows, "source_path") <- csv_path
    attr(rows, "fingerprint") <- digest::digest(
        file = csv_path,
        algo = "sha256"
    )

    rows
}


#' Load the azelate Panel D source table
#'
#' Reads
#' \code{data/derived/azelate/csvs/cancer_assoc_by_sample_type.csv}
#' — the resolved-only cancer-association table grouped by disease
#' type and sample (tissue) type, produced by the same azelate
#' extraction notebook as Panel C. Drives Figure 6 Panel D.
#'
#' @param csv_path Character: path to \code{cancer_assoc_by_sample_type.csv}.
#'     Defaults to the repo-relative vendored location.
#'
#' @return A tibble with columns \code{disease_type},
#'     \code{sample_type}, \code{evidence_count}. Carries
#'     \code{source_path} + SHA-256 \code{fingerprint} attributes.
#'
#' @examples
#' \dontrun{
#' azelate_panel_d_data()
#' }
#'
#' @importFrom readr read_csv cols col_character col_integer
#' @importFrom rlang abort
#' @importFrom digest digest
#' @export
azelate_panel_d_data <- function(
    csv_path = file.path(
        "data",
        "derived",
        "azelate",
        "csvs",
        "cancer_assoc_by_sample_type.csv"
    )
) {

    if (!file.exists(csv_path)) {
        rlang::abort(sprintf(
            paste0(
                "Azelate Panel D source CSV not found at %s. ",
                "Regenerate via azelaic_acid_query_tables.py per ",
                "analyses/azelate/README.md."
            ),
            csv_path
        ))
    }

    rows <- readr::read_csv(
        csv_path,
        col_types = readr::cols(
            disease_type   = readr::col_character(),
            sample_type    = readr::col_character(),
            evidence_count = readr::col_integer()
        ),
        progress = FALSE
    )

    attr(rows, "source_path") <- csv_path
    attr(rows, "fingerprint") <- digest::digest(
        file = csv_path,
        algo = "sha256"
    )

    rows
}
