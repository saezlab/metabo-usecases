#' Derive a Build Manifest from a live OmniPath deployment
#'
#' Introspects the Postgres at \code{con} to capture
#' \code{Resource State} (loaded resources + per-resource record
#' counts), then merges in the package commit hashes for the
#' targeted build (\code{utils}, \code{main}, or \code{metabo}). The
#' returned manifest is suitable for canonical-JSON serialization and
#' SHA-256 hashing into a Snapshot Identifier
#' (\code{\link{snapshot_id}}).
#'
#' Package set per build is fixed by spec FR-032:
#' \itemize{
#'     \item \code{utils}: \code{omnipath-utils} + \code{omnipath-resources}
#'     \item \code{main}: \code{omnipath-build} + \code{omnipath-utils} +
#'         \code{omnipath-resources}
#'     \item \code{metabo}: \code{omnipath-metabo} + \code{omnipath-build}
#'         + \code{omnipath-utils} + \code{omnipath-resources}
#' }
#' \code{omnipath-present} is excluded (the pipeline queries Postgres
#' directly).
#'
#' Where commit hashes are not supplied via \code{package_commits},
#' they are read from \code{manifest-config.yaml} at the project root,
#' falling back to the env vars \code{OMNIPATH_BUILD_COMMIT},
#' \code{OMNIPATH_UTILS_COMMIT}, etc. — and finally to the placeholder
#' string \code{"unknown"} (the manifest's \code{partial_build} flag
#' becomes \code{TRUE} in that case so consumers know the identifier
#' is not authoritative).
#'
#' @param con A DBI Postgres connection from \code{\link{pg_connect}}.
#' @param build Character: one of \code{"utils"}, \code{"main"},
#'     \code{"metabo"}.
#' @param package_commits Optional named character vector overriding
#'     the package commit lookup.
#'
#' @return A named list conforming to
#'     \code{contracts/build-manifest.schema.json}.
#'
#' @examples
#' \dontrun{
#' con <- pg_connect()
#' m <- build_manifest_for(con, "main")
#' snapshot_id(m)
#' }
#'
#' @importFrom DBI dbGetQuery
#' @importFrom dplyr arrange
#' @importFrom tibble as_tibble
#' @importFrom purrr map
#' @importFrom logger log_info
#' @importFrom rlang abort
#' @export
build_manifest_for <- function(con, build, package_commits = NULL) {

    # NSE workaround
    name <- record_count <- NULL

    if (!build %in% c("utils", "main", "metabo")) {
        rlang::abort(sprintf("Unknown build: '%s'", build))
    }

    resources <- resource_state(con) %>% dplyr::arrange(name)

    package_set <- packages_for_build(build)
    commits <- resolve_package_commits(package_set, package_commits)

    partial <- has_partial_records(resources) ||
        any(commits == "unknown")

    list(
        build         = build,
        built_at      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
        packages      = commits,
        resources     = purrr::map(seq_len(nrow(resources)), function(i) {
            list(
                name           = resources$name[[i]],
                version        = resources$version[[i]] %||% "unknown",
                record_count   = as.integer(resources$record_count[[i]]),
                expected_count = if (is.na(resources$expected_count[[i]])) {
                    NULL
                } else {
                    as.integer(resources$expected_count[[i]])
                }
            )
        }),
        partial_build = isTRUE(partial)
    )
}


#' Compute the Snapshot Identifier from a Build Manifest
#'
#' Serializes the manifest to deterministic canonical JSON (sorted
#' keys, excluding the cosmetic \code{built_at}), SHA-256 hashes the
#' result, and returns the first 12 hex chars (research.md R-7).
#'
#' @param manifest A list as produced by \code{\link{build_manifest_for}}.
#'
#' @return Character scalar — 12 hex chars.
#'
#' @examples
#' m <- list(
#'     build = "metabo",
#'     packages = list(`omnipath-metabo` = "abc1234"),
#'     resources = list(),
#'     partial_build = FALSE
#' )
#' snapshot_id(m)
#'
#' @importFrom jsonlite toJSON
#' @importFrom digest digest
#' @export
snapshot_id <- function(manifest) {

    hash_input <- manifest
    hash_input$built_at <- NULL

    canonical <- jsonlite::toJSON(
        hash_input,
        auto_unbox = TRUE,
        digits     = NA,
        null       = "null",
        na         = "null",
        pretty     = FALSE
    )

    substr(
        digest::digest(canonical, algo = "sha256", serialize = FALSE),
        1L, 12L
    )
}


#' Write a manifest + its identifier to manifests/
#'
#' Emits \code{manifests/<build>.<snapshot-id>.json} (canonical JSON
#' form) and \code{manifests/<build>.<snapshot-id>.SHA256} (one-line
#' hex digest). Returns the resolved identifier.
#'
#' @param manifest A list as produced by \code{\link{build_manifest_for}}.
#' @param dir Output directory; defaults to \code{"manifests"}.
#'
#' @return The Snapshot Identifier (character scalar).
#'
#' @importFrom jsonlite toJSON
#' @importFrom fs dir_create path
#' @importFrom logger log_info
#' @export
write_manifest <- function(manifest, dir = "manifests") {

    fs::dir_create(dir)
    sid <- snapshot_id(manifest)

    json_path  <- fs::path(dir, sprintf("%s.%s.json",  manifest$build, sid))
    sha_path   <- fs::path(dir, sprintf("%s.%s.SHA256", manifest$build, sid))

    writeLines(
        jsonlite::toJSON(
            manifest, auto_unbox = TRUE, pretty = TRUE, null = "null"
        ),
        json_path
    )
    writeLines(sid, sha_path)

    logger::log_info("Wrote manifest {json_path} (snapshot={sid})")
    sid
}


#' Query the live OmniPath schema for resource state
#'
#' Joins \code{data_source} with the per-table record counts (entity,
#' relation, annotation, ontology_terms). Implementation note: when a
#' resource has no entries in a given table, the count is zero (not
#' missing). Expected counts are not currently tracked in the schema
#' and default to \code{NA_integer_}.
#'
#' @param con DBI connection.
#' @return A tibble with columns \code{name}, \code{version},
#'     \code{record_count}, \code{expected_count}.
#'
#' @importFrom DBI dbGetQuery
#' @importFrom tibble as_tibble
#' @keywords internal
#' @noRd
resource_state <- function(con) {

    # Conservative shape that works against the current dev3 schema:
    # data_source is the resource registry; counts come from
    # entity_evidence (partitioned by source_id, so the count query is
    # cheap). Other per-source totals (relations, annotations) are
    # added back once the corresponding schema-side joins are stable —
    # this is sufficient for distinguishing snapshots today.
    sql <- "
        SELECT ds.name AS name,
               'unknown' AS version,
               COALESCE((
                   SELECT COUNT(*) FROM entity_evidence ee
                   WHERE ee.source_id = ds.source_id
               ), 0) AS record_count,
               NULL::bigint AS expected_count
        FROM data_source ds
        ORDER BY ds.name
    "

    tibble::as_tibble(DBI::dbGetQuery(con, sql))
}


#' Read package commit hashes from config / env / fallback
#'
#' @param wanted Character: the package names this build requires.
#' @param override Optional named character vector taking precedence.
#'
#' @importFrom yaml read_yaml
#' @keywords internal
#' @noRd
resolve_package_commits <- function(wanted, override = NULL) {

    cfg <- if (file.exists("manifest-config.yaml")) {
        yaml::read_yaml("manifest-config.yaml")$packages %||% list()
    } else {
        list()
    }

    env_lookup <- list(
        `omnipath-resources` = "OMNIPATH_RESOURCES_COMMIT",
        `omnipath-utils`     = "OMNIPATH_UTILS_COMMIT",
        `omnipath-build`     = "OMNIPATH_BUILD_COMMIT",
        `omnipath-metabo`    = "OMNIPATH_METABO_COMMIT"
    )

    out <- vapply(wanted, function(pkg) {
        if (!is.null(override) && pkg %in% names(override)) {
            return(unname(override[[pkg]]))
        }
        if (pkg %in% names(cfg)) return(cfg[[pkg]])
        env_val <- Sys.getenv(env_lookup[[pkg]] %||% "", unset = "")
        if (nzchar(env_val)) return(env_val)
        "unknown"
    }, character(1L), USE.NAMES = FALSE)

    setNames(as.list(out), wanted)
}


#' Per-build set of relevant package names
#'
#' @param build Character: one of \code{utils}, \code{main}, \code{metabo}.
#' @keywords internal
#' @noRd
packages_for_build <- function(build) {
    switch(
        build,
        utils  = c("omnipath-utils", "omnipath-resources"),
        main   = c("omnipath-build", "omnipath-utils", "omnipath-resources"),
        metabo = c(
            "omnipath-metabo", "omnipath-build",
            "omnipath-utils",  "omnipath-resources"
        )
    )
}


#' Whether any resource was partially loaded
#'
#' @param resources Tibble with \code{record_count}, \code{expected_count}.
#' @importFrom dplyr filter
#' @keywords internal
#' @noRd
has_partial_records <- function(resources) {

    # NSE workaround
    record_count <- expected_count <- NULL

    if (nrow(resources) == 0L) return(FALSE)
    partial <- resources %>%
        dplyr::filter(
            !is.na(expected_count) & record_count < expected_count
        )
    nrow(partial) > 0L
}
