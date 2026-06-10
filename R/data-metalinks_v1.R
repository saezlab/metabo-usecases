#' Load the MetaLinksDB v1 baseline via OmnipathR
#'
#' Reads the SQLite snapshot shipped through `OmnipathR`, joins the
#' `edges`, `metabolites`, and `proteins` tables, and returns a
#' harmonized Figure 3 snapshot on the HMDB / UniProt counting basis.
#'
#' @return List with `data` (normalized tibble) and `external_input`
#'     (metadata for the provenance sidecar).
#'
#' @importFrom DBI dbDisconnect dbGetInfo
#' @importFrom digest digest
#' @importFrom dplyr select
#' @importFrom rlang check_installed
#' @export
metalinks_v1_snapshot <- function() {

    rlang::check_installed("OmnipathR")
    rlang::check_installed("RSQLite")

    sqlite_con <- OmnipathR::metalinksdb_sqlite()
    on.exit(DBI::dbDisconnect(sqlite_con), add = TRUE)

    sqlite_path <- normalizePath(
        DBI::dbGetInfo(sqlite_con)$dbname,
        winslash = "/",
        mustWork = FALSE
    )

    edges <- OmnipathR::metalinksdb_table("edges")
    metabolites <- OmnipathR::metalinksdb_table("metabolites")
    proteins <- OmnipathR::metalinksdb_table("proteins")

    normalized <- normalize_mpi_resource(
        resource = "MetaLinksDB v1.0",
        interactions = edges,
        metabolite_key = "hmdb",
        protein_key = "uniprot",
        source_col = "source",
        relation_type_col = "type",
        metabolite_class = dplyr::select(
            metabolites,
            hmdb,
            metabolite_subclass
        ),
        metabolite_class_key = "hmdb",
        metabolite_class_col = "metabolite_subclass",
        protein_class = dplyr::select(
            proteins,
            uniprot,
            protein_type
        ),
        protein_class_key = "uniprot",
        protein_class_col = "protein_type",
        evidence_cols = list(
            source_count = NULL,
            citation_count = NULL,
            affinity_value = NULL,
            curation_mode = NULL
        ),
        interaction_definition = "one HMDB-UniProt-source-type row from the OmnipathR MetaLinksDB v1 SQLite"
    )

    list(
        data = normalized,
        external_input = list(
            kind = "metalinksdb-v1-sqlite",
            path = sqlite_path,
            source = sprintf(
                "OmnipathR %s",
                as.character(utils::packageVersion("OmnipathR"))
            ),
            fingerprint = digest::digest(
                file = sqlite_path,
                algo = "sha256"
            ),
            version = as.character(utils::packageVersion("OmnipathR"))
        )
    )
}
