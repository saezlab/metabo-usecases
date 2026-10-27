#' Read a Build Manifest natively from a live OmniPath deployment
#'
#' Reads the single-row \code{build_manifest} table on the target
#' deployment (cycle-001 deliverable from
#' \code{omnipath-build}). Returns a list shaped to
#' \code{contracts/build-manifest.schema.json}. The
#' \code{build_id} (12-hex SHA-256 truncation of canonicalised
#' \code{package_commits + resources}) is read from the column, not
#' re-derived — by construction it equals what the FR-032a
#' standalone extractor would compute.
#'
#' Per the 2026-06-09 cycle-001 + 002 handover, every deployment the
#' figure pipeline targets (\code{dev3} + \code{dev4}) carries this
#' table; the legacy post-hoc derivation track (the original FR-032b
#' fallback) is not implemented here per developer's direction.
#'
#' @param con A DBI Postgres connection from
#'     \code{\link{pg_connect_panel}} (or \code{\link{pg_connect}}).
#'
#' @return A named list with elements \code{build_id},
#'     \code{built_at}, \code{build}, \code{packages},
#'     \code{resources}, \code{partial_build}. The \code{build} field
#'     is inferred from \code{packages}: \code{"metabo"} when
#'     \code{omnipath-metabo} is present; \code{"main"} when
#'     \code{omnipath-build} is present; \code{"utils"} otherwise.
#'
#' @examples
#' \dontrun{
#' con <- pg_connect_panel("dev3")
#' m <- build_manifest_for(con)
#' m$build_id        # e.g. "a3f9c2e74b81"
#' }
#'
#' @importFrom DBI dbGetQuery
#' @importFrom jsonlite fromJSON
#' @importFrom logger log_info
#' @importFrom rlang abort
#' @export
build_manifest_for <- function(con) {

    rows <- DBI::dbGetQuery(con, "SELECT * FROM build_manifest LIMIT 1")

    if (nrow(rows) == 0L) {
        rlang::abort(paste0(
            "build_manifest table is empty on the target deployment — ",
            "expected the cycle-001 single-row schema. Re-run the ",
            "database build's manifest-emit step."
        ))
    }

    packages_raw <- parse_jsonb(rows$package_commits[[1L]])
    packages <- flatten_package_commits(packages_raw)
    resources <- parse_jsonb(rows$resources[[1L]])
    build_kind <- infer_build_kind(packages)

    manifest <- list(
        build_id      = as.character(rows$build_id[[1L]]),
        built_at      = format(rows$built_at[[1L]], "%Y-%m-%dT%H:%M:%S%z"),
        build         = build_kind,
        packages      = packages,
        packages_raw  = packages_raw,
        resources     = resources,
        partial_build = isTRUE(as.logical(rows$partial_build[[1L]]))
    )

    logger::log_info(
        "build_manifest read: build_id={manifest$build_id} ",
        "(build={build_kind}, resources={length(resources)}, ",
        "partial_build={manifest$partial_build})"
    )

    manifest
}


#' Snapshot Identifier accessor
#'
#' Returns the manifest's \code{build_id}. The cycle-001
#' \code{build_manifest} column is defined as the SHA-256 of
#' canonicalised \code{package_commits + resources} truncated to 12
#' hex chars (FR-032a), so this is the Snapshot Identifier by
#' construction — no in-pipeline re-derivation.
#'
#' The function is kept for backwards compatibility; new code should
#' read \code{manifest$build_id} directly.
#'
#' @param manifest A list as produced by \code{\link{build_manifest_for}}.
#'
#' @return Character scalar — 12 hex chars.
#'
#' @examples
#' \dontrun{
#' snapshot_id(build_manifest_for(con))
#' }
#'
#' @export
snapshot_id <- function(manifest) {

    if (is.null(manifest$build_id)) {
        rlang::abort(
            "Manifest is missing build_id (post-cycle-001 schema)"
        )
    }

    as.character(manifest$build_id)
}


#' Archive a manifest + its identifier to build/manifests/
#'
#' Emits \code{build/manifests/<deployment>.<build_id>.json} (canonical JSON
#' form) and \code{build/manifests/<deployment>.<build_id>.SHA256} (one-line
#' hex digest). Returns the build_id.
#'
#' @param manifest A list as produced by \code{\link{build_manifest_for}}.
#' @param deployment Character: deployment label that produced the
#'     manifest (e.g. \code{"dev3"}, \code{"dev4"}). Used as the
#'     filename prefix so a multi-deployment rebuild keeps each
#'     archive distinct.
#' @param dir Output directory; defaults to \code{"build/manifests"}.
#'
#' @return The build_id (character scalar).
#'
#' @importFrom jsonlite toJSON
#' @importFrom fs dir_create path
#' @importFrom logger log_info
#' @export
write_manifest <- function(manifest, deployment,
                           dir = file.path("build", "manifests")) {

    fs::dir_create(dir)
    sid <- snapshot_id(manifest)

    json_path <- fs::path(
        dir, sprintf("%s.%s.json", deployment, sid)
    )
    sha_path  <- fs::path(
        dir, sprintf("%s.%s.SHA256", deployment, sid)
    )

    # Archive the FAITHFUL DB shape — rich packages (with dirty flag
    # when present). The sidecar's package_commits stays flat (per the
    # FR-032 schema), the archive preserves the dirty-flag detail.
    archive <- manifest
    if (!is.null(manifest$packages_raw)) {
        archive$packages <- manifest$packages_raw
    }
    archive$packages_raw <- NULL

    writeLines(
        jsonlite::toJSON(
            archive, auto_unbox = TRUE, pretty = TRUE, null = "null"
        ),
        json_path
    )
    writeLines(sid, sha_path)

    logger::log_info(
        "Wrote manifest {json_path} (build_id={sid})"
    )

    sid
}


#' Per-build set of relevant package names (FR-032)
#'
#' @param build Character: one of \code{utils}, \code{main}, \code{metabo}.
#'
#' @return Character vector of package names.
#'
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


#' Infer build kind from a package_commits map
#'
#' Mirrors the per-build package sets in FR-032: \code{metabo} when
#' \code{omnipath-metabo} appears; otherwise \code{main} when
#' \code{omnipath-build} appears; otherwise \code{utils}. The
#' canonical FR-032 form uses hyphens (\code{omnipath-build}); the
#' live cycle-001 \code{build_manifest} stores keys with underscores
#' (\code{omnipath_build}). Both forms are accepted — keys are
#' normalised to hyphen-form before classification.
#'
#' @param packages Named list of package → commit hash.
#'
#' @return Character scalar: \code{"utils"} | \code{"main"} |
#'     \code{"metabo"}.
#'
#' @keywords internal
#' @noRd
infer_build_kind <- function(packages) {

    names_set <- gsub("_", "-", names(packages))
    if ("omnipath-metabo" %in% names_set) {
        "metabo"
    } else if ("omnipath-build" %in% names_set) {
        "main"
    } else {
        "utils"
    }
}


#' Flatten the build_manifest.package_commits column to name → hash
#'
#' The cycle-001 \code{build_manifest} stores each commit as
#' \code{{"commit": "<sha>", "dirty": <bool>}}; the FR-032 contract
#' (build-manifest.schema.json) and the provenance sidecar schema
#' both require a flat name → hash string mapping. This helper
#' flattens the rich shape to the contract's flat shape so the
#' sidecar's \code{package_commits} field is schema-compliant. The
#' rich shape is preserved as \code{packages_raw} on the manifest so
#' the archived manifest JSON in \code{build/manifests/} keeps the
#' \code{dirty} flag.
#'
#' Tolerates the legacy flat shape: when a value is already a
#' character scalar, it is returned unchanged.
#'
#' @param packages_raw A named list as produced by
#'     \code{\link{parse_jsonb}} over the \code{package_commits}
#'     column.
#'
#' @return Named list of name → commit hash string.
#'
#' @keywords internal
#' @noRd
flatten_package_commits <- function(packages_raw) {

    lapply(packages_raw, function(value) {
        if (is.character(value) && length(value) == 1L) {
            value
        } else if (is.list(value) && !is.null(value$commit)) {
            as.character(value$commit)
        } else {
            # Defensive fall-through: jsonlite serializer represents
            # this as the string "unknown" so downstream schema
            # validation can still pass.
            "unknown"
        }
    })
}


#' Parse a jsonb column value into an R list
#'
#' RPostgres returns jsonb as character; this helper trims whitespace,
#' parses with \code{simplifyVector = FALSE}, and returns an empty
#' list on null / empty input.
#'
#' @param value Character or list (raw column value).
#'
#' @return List — empty when the column was null / empty.
#'
#' @importFrom jsonlite fromJSON
#' @keywords internal
#' @noRd
parse_jsonb <- function(value) {

    if (is.null(value) || (is.character(value) && !nzchar(value))) {
        return(list())
    }

    if (is.list(value)) {
        return(value)
    }

    jsonlite::fromJSON(as.character(value), simplifyVector = FALSE)
}
