#' Normalize a metabolite-protein interaction resource snapshot
#'
#' Harmonizes a resource-specific interaction table onto the Figure 3
#' comparison basis: HMDB ids for metabolites, UniProt ids for proteins,
#' one row per interaction observation, explicit source / relation-type
#' columns, and optional metabolite / protein class annotations.
#'
#' @param resource Character scalar naming the resource.
#' @param interactions Tibble or data.frame containing the interaction rows.
#' @param metabolite_key Character scalar naming the HMDB-id column.
#' @param protein_key Character scalar naming the UniProt-id column.
#' @param source_col Optional character scalar naming the source column.
#' @param relation_type_col Optional character scalar naming the relation-type column.
#' @param metabolite_class Optional lookup table with metabolite classes.
#' @param metabolite_class_key Character scalar naming the HMDB-id column in
#'     `metabolite_class`.
#' @param metabolite_class_col Character scalar naming the class column in
#'     `metabolite_class`.
#' @param protein_class Optional lookup table with protein classes.
#' @param protein_class_key Character scalar naming the UniProt-id column in
#'     `protein_class`.
#' @param protein_class_col Character scalar naming the class column in
#'     `protein_class`.
#' @param evidence_cols Named list describing optional evidence columns:
#'     `source_count`, `citation_count`, `affinity_value`, `curation_mode`.
#' @param interaction_definition Character scalar describing what counts as one
#'     interaction row for this resource.
#'
#' @return Tibble on the harmonized Figure 3 basis.
#'
#' @importFrom dplyr left_join mutate filter transmute group_by summarise n_distinct
#' @importFrom tibble as_tibble tibble
#' @importFrom rlang abort .data
#' @export
normalize_mpi_resource <- function(resource,
                                   interactions,
                                   metabolite_key,
                                   protein_key,
                                   source_col = NULL,
                                   relation_type_col = NULL,
                                   metabolite_class = NULL,
                                   metabolite_class_key = NULL,
                                   metabolite_class_col = NULL,
                                   protein_class = NULL,
                                   protein_class_key = NULL,
                                   protein_class_col = NULL,
                                   evidence_cols = list(),
                                   interaction_definition = 'hmdb-uniprot pair') {

    hmdb_id <- uniprot_id <- source <- relation_type <- metabolite_class_label <- NULL
    protein_class_label <- source_count <- citation_count <- affinity_value <- NULL
    curation_mode <- interaction_id <- NULL

    interactions <- tibble::as_tibble(interactions)

    required <- c(metabolite_key, protein_key)
    missing_required <- setdiff(required, names(interactions))

    if (length(missing_required) > 0L) {
        rlang::abort(sprintf(
            "Resource '%s' is missing required columns: %s",
            resource,
            paste(missing_required, collapse = ', ')
        ))
    }

    data <- dplyr::transmute(
        interactions,
        hmdb_id = .data[[metabolite_key]],
        uniprot_id = .data[[protein_key]],
        source = if (!is.null(source_col) && source_col %in% names(interactions)) {
            .data[[source_col]]
        } else {
            resource
        },
        relation_type = if (!is.null(relation_type_col) &&
            relation_type_col %in% names(interactions)) {
            .data[[relation_type_col]]
        } else {
            'interaction'
        },
        source_count = normalize_numeric_column(interactions, evidence_cols$source_count, 1),
        citation_count = normalize_numeric_column(interactions, evidence_cols$citation_count, NA_real_),
        affinity_value = normalize_numeric_column(interactions, evidence_cols$affinity_value, NA_real_),
        curation_mode = normalize_character_column(interactions, evidence_cols$curation_mode, NA_character_)
    )

    data <- dplyr::mutate(
        data,
        resource = resource,
        hmdb_id = as.character(.data$hmdb_id),
        uniprot_id = as.character(.data$uniprot_id),
        source = as.character(.data$source),
        relation_type = as.character(.data$relation_type),
        curation_mode = as.character(.data$curation_mode),
        interaction_definition = interaction_definition
    )
    data <- dplyr::filter(data, !is.na(.data$hmdb_id), !is.na(.data$uniprot_id))

    if (!is.null(metabolite_class)) {
        metabolite_class <- tibble::as_tibble(metabolite_class)
        met_lookup <- dplyr::distinct(
            dplyr::transmute(
                metabolite_class,
                hmdb_id = .data[[metabolite_class_key]],
                metabolite_class_label = .data[[metabolite_class_col]]
            ),
            .data$hmdb_id,
            .keep_all = TRUE
        )
        data <- dplyr::left_join(data, met_lookup, by = 'hmdb_id')
    } else {
        data$metabolite_class_label <- NA_character_
    }

    if (!is.null(protein_class)) {
        protein_class <- tibble::as_tibble(protein_class)
        prot_lookup <- dplyr::distinct(
            dplyr::transmute(
                protein_class,
                uniprot_id = .data[[protein_class_key]],
                protein_class_label = .data[[protein_class_col]]
            ),
            .data$uniprot_id,
            .keep_all = TRUE
        )
        data <- dplyr::left_join(data, prot_lookup, by = 'uniprot_id')
    } else {
        data$protein_class_label <- NA_character_
    }

    dplyr::mutate(
        data,
        metabolite_class_label = fill_missing_label(.data$metabolite_class_label),
        protein_class_label = fill_missing_label(.data$protein_class_label),
        relation_type = fill_missing_label(.data$relation_type),
        source = fill_missing_label(.data$source),
        interaction_id = paste(.data$hmdb_id, .data$uniprot_id, .data$relation_type, .data$source, sep = '::')
    )
}


#' Load a vendored MPI baseline snapshot from data/vendored
#'
#' Reads one of the manually-normalized Figure 3 baseline files from
#' `data/vendored/mpi-baselines/`. Supported file formats are CSV,
#' TSV, and JSON.
#'
#' @param resource Character scalar naming the resource.
#' @param root Character path to the snapshot directory.
#' @param required Logical: if `TRUE`, abort when the resource is missing.
#'
#' @return Tibble of normalized rows, or an empty tibble when missing and
#'     `required = FALSE`.
#'
#' @importFrom fs path file_exists
#' @importFrom readr read_csv read_tsv
#' @importFrom jsonlite fromJSON
#' @importFrom tibble as_tibble tibble
#' @importFrom rlang abort
#' @export
load_vendored_mpi_snapshot <- function(resource,
                                       root = fig03_baseline_dir(),
                                       required = FALSE) {

    slug <- snapshot_slug(resource)
    candidates <- c(
        fs::path(root, paste0(slug, '.csv')),
        fs::path(root, paste0(slug, '.tsv')),
        fs::path(root, paste0(slug, '.json'))
    )
    found <- candidates[fs::file_exists(candidates)]

    if (length(found) == 0L) {
        if (isTRUE(required)) {
            rlang::abort(sprintf(
                "Vendored Figure 3 baseline snapshot for '%s' is missing in %s",
                resource,
                root
            ))
        }

        return(tibble::tibble())
    }

    path <- found[[1L]]
    if (grepl('\\.csv$', path)) {
        readr::read_csv(path, show_col_types = FALSE)
    } else if (grepl('\\.tsv$', path)) {
        readr::read_tsv(path, show_col_types = FALSE)
    } else {
        tibble::as_tibble(jsonlite::fromJSON(path))
    }
}


#' Directory containing vendored Figure 3 baseline snapshots
#'
#' @return Character scalar path.
#' @keywords internal
#' @noRd
fig03_baseline_dir <- function() {
    file.path('data', 'vendored', 'mpi-baselines')
}


#' Summarize harmonized MPI rows to the Figure 3 coverage basis
#'
#' @param data Tibble from `normalize_mpi_resource()`.
#'
#' @return Tibble with one row per resource.
#' @importFrom dplyr group_by summarise n_distinct
#' @importFrom rlang .data
#' @export
summarize_mpi_coverage <- function(data) {

    resource <- interaction_id <- hmdb_id <- uniprot_id <- NULL

    dplyr::summarise(
        dplyr::group_by(data, .data$resource),
        n_interactions = dplyr::n_distinct(.data$interaction_id),
        n_metabolites = dplyr::n_distinct(.data$hmdb_id),
        n_proteins = dplyr::n_distinct(.data$uniprot_id),
        .groups = 'drop'
    )
}


#' Normalize a numeric evidence column or fall back to a default
#'
#' @param data Input interaction table.
#' @param column Column name or NULL.
#' @param default Default numeric value.
#'
#' @return Numeric vector.
#' @keywords internal
#' @noRd
normalize_numeric_column <- function(data, column, default) {
    if (!is.null(column) && column %in% names(data)) {
        as.numeric(data[[column]])
    } else {
        rep(default, nrow(data))
    }
}


#' Normalize a character evidence column or fall back to a default
#'
#' @param data Input interaction table.
#' @param column Column name or NULL.
#' @param default Default character value.
#'
#' @return Character vector.
#' @keywords internal
#' @noRd
normalize_character_column <- function(data, column, default) {
    if (!is.null(column) && column %in% names(data)) {
        as.character(data[[column]])
    } else {
        rep(default, nrow(data))
    }
}


#' Replace missing/blank class labels with an explicit fallback
#'
#' @param x Character vector.
#' @return Character vector.
#' @keywords internal
#' @noRd
fill_missing_label <- function(x) {
    x <- as.character(x)
    x[is.na(x) | trimws(x) == ''] <- 'Unclassified'
    x
}


#' Slugify a resource label to a baseline snapshot filename stem
#'
#' @param resource Character scalar.
#' @return Character scalar.
#' @keywords internal
#' @noRd
snapshot_slug <- function(resource) {
    slug <- tolower(gsub('[^a-zA-Z0-9]+', '-', resource))
    gsub('(^-+|-+$)', '', slug)
}
