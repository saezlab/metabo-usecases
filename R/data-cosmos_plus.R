#' Load COSMOS+ prior-knowledge network data from vendored CSVs
#'
#' Reads the vendored COSMOS+ CSVs produced by T052a from
#' \code{inst/extdata/cosmos/cosmos_plus_human.csv} and
#' \code{inst/extdata/cosmos/cosmos_plus_mouse.csv} and derives three
#' shapes for Figure 4:
#'
#' \enumerate{
#'   \item \strong{by_type_species}: interaction counts by
#'     \code{interaction_type} and \code{species} (Panel A).
#'   \item \strong{by_compartment}: interaction counts per compartment
#'     (Panel B); if all edges have empty \code{locations}, returns a
#'     one-row tibble with \code{compartment = "unannotated"}.
#'   \item \strong{by_resource}: unique metabolite (ChEBI), protein
#'     (UniProt), and interaction counts per \code{resource} (Panel C).
#' }
#'
#' SHA-256 fingerprints of both input CSVs are recorded as
#' \code{external_inputs} attributes for the provenance sidecar.
#'
#' @param human_csv Path to the human COSMOS+ CSV.
#' @param mouse_csv Path to the mouse COSMOS+ CSV. Pass \code{NULL} to
#'     skip mouse and report human-only shapes.
#'
#' @return A named list with elements \code{by_type_species},
#'     \code{by_compartment}, \code{by_resource}, and
#'     \code{external_inputs} (a list of two sidecar records).
#'
#' @importFrom readr read_csv
#' @importFrom dplyr filter mutate bind_rows count n_distinct group_by
#'     summarise arrange desc
#' @importFrom tidyr separate_rows
#' @importFrom tibble tibble
#' @importFrom digest digest
#' @importFrom rlang abort
#' @export
cosmos_plus_data <- function(
    human_csv = cosmos_plus_csv_path("human"),
    mouse_csv = cosmos_plus_csv_path("mouse")
) {
    interaction_type <- species <- n_interactions <- NULL
    locations <- compartment <- resource <- NULL
    source_type <- target_type <- source <- target <- NULL

    if (!file.exists(human_csv)) {
        rlang::abort(paste0(
            "COSMOS+ human CSV not found at: ", human_csv, ". ",
            "Re-vendor by running T052a: cosmos-pkn export --all-columns ",
            "--organism 9606 --output ", human_csv
        ))
    }

    human <- readr::read_csv(human_csv, show_col_types = FALSE) |>
        dplyr::filter(interaction_type != "connector") |>
        dplyr::mutate(species = "human")

    human_fp <- digest::digest(file = human_csv, algo = "sha256")

    if (!is.null(mouse_csv) && file.exists(mouse_csv)) {
        mouse <- readr::read_csv(mouse_csv, show_col_types = FALSE) |>
            dplyr::filter(interaction_type != "connector") |>
            dplyr::mutate(species = "mouse")
        mouse_fp <- digest::digest(file = mouse_csv, algo = "sha256")
        combined <- dplyr::bind_rows(human, mouse)
    } else {
        mouse_fp <- NA_character_
        combined <- human
    }

    # Shape (a): counts by interaction_type × species
    by_type_species <- combined |>
        dplyr::count(interaction_type, species, name = "n_interactions") |>
        dplyr::arrange(dplyr::desc(n_interactions))

    # Shape (b): counts per compartment
    locs <- combined$locations
    has_locations <- any(!is.na(locs) & nzchar(trimws(locs)))
    if (has_locations) {
        by_compartment <- combined |>
            dplyr::filter(!is.na(locations), nzchar(trimws(locations))) |>
            tidyr::separate_rows(locations, sep = ",\\s*") |>
            dplyr::filter(nzchar(trimws(locations))) |>
            dplyr::mutate(compartment = trimws(locations)) |>
            dplyr::count(compartment, name = "n_interactions") |>
            dplyr::arrange(dplyr::desc(n_interactions))
    } else {
        by_compartment <- tibble::tibble(
            compartment    = "unannotated",
            n_interactions = nrow(combined)
        )
    }

    # Shape (c): metabolites, proteins, interactions per resource
    by_resource <- combined |>
        dplyr::group_by(resource) |>
        dplyr::summarise(
            n_metabolites  = dplyr::n_distinct(
                source[source_type == "metabolite"],
                target[target_type == "metabolite"]
            ),
            n_proteins     = dplyr::n_distinct(
                source[source_type == "protein"],
                target[target_type == "protein"]
            ),
            n_interactions = dplyr::n_distinct(source, target),
            .groups        = "drop"
        ) |>
        dplyr::arrange(dplyr::desc(n_interactions))

    external_inputs <- list(
        list(
            kind        = "cosmos-plus-csv",
            species     = "human",
            path        = human_csv,
            fingerprint = human_fp
        ),
        list(
            kind        = "cosmos-plus-csv",
            species     = "mouse",
            path        = mouse_csv %||% "not-vendored",
            fingerprint = mouse_fp
        )
    )

    list(
        by_type_species  = by_type_species,
        by_compartment   = by_compartment,
        by_resource      = by_resource,
        external_inputs  = external_inputs
    )
}


#' Resolve the vendored COSMOS+ CSV path
#'
#' @param species \code{"human"} or \code{"mouse"}.
#' @return Character scalar path.
#' @keywords internal
#' @noRd
cosmos_plus_csv_path <- function(species = c("human", "mouse")) {
    species <- match.arg(species)
    filename <- paste0("cosmos_plus_", species, ".csv")
    installed <- system.file(
        "extdata", "cosmos", filename,
        package = "metabo.figures"
    )
    if (nzchar(installed) && file.exists(installed)) {
        return(installed)
    }
    file.path("inst", "extdata", "cosmos", filename)
}
