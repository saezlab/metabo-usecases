#' Write a provenance sidecar for a pipeline artifact
#'
#' Serializes the run-level + per-query inputs into the JSON shape
#' documented in
#' \code{specs/001-figures-pipeline/contracts/provenance-sidecar.schema.json}.
#' Every pipeline figure / table / extra plot emits one sidecar
#' alongside its output file at \code{<artifact>.provenance.json}.
#'
#' Per the 2026-06-09 handover the sidecar carries
#' \code{deployments: array} — a panel that combines records from
#' \code{dev3} and \code{dev4} (e.g. FR-007a Entities facet from dev3
#' joined with the Structures facet from dev4) records BOTH
#' deployments here. Each entry carries the \code{build_id} read
#' natively from that deployment's \code{build_manifest} table
#' (FR-032, FR-032a).
#'
#' @param artifact_id Character: stable id, e.g.
#'     \code{"fig01-overview/panelB"}.
#' @param artifact_path Character: path of the artifact whose sidecar
#'     this is — the sidecar file is written to
#'     \code{paste0(artifact_path, ".provenance.json")}.
#' @param deployments List of deployment-provenance records (see
#'     \code{\link{deployment_provenance}}). Each entry has
#'     \code{name}, \code{url}, \code{db_host}, \code{db_port},
#'     \code{db_name}, \code{build_id}, and optionally
#'     \code{build_manifest_source}. At least one entry is required.
#' @param manifests List of build manifests (one per deployment in
#'     \code{deployments}, parallel index). The sidecar's
#'     \code{package_commits} is the union of every manifest's
#'     \code{packages}; \code{partial_build = any}.
#' @param script_path Character: path to the generating script
#'     relative to the repo root.
#' @param queries List of query records — each list element has
#'     \code{sql}, \code{row_count}, \code{result_hash}. Typically
#'     accumulated via the \code{attr(rows, "sql")} +
#'     \code{attr(rows, "result_hash")} produced by \code{\link{pg_query}}
#'     or \code{\link{pg_query_panel}}.
#' @param external_inputs Optional list of non-Postgres inputs (manual
#'     assets, palettes, vendored RData, the Figure 1 architecture
#'     asset). Each entry has \code{kind}, \code{path}, optional
#'     \code{source}, and \code{fingerprint}.
#' @param parameters Optional named list of filter / contrast /
#'     threshold values used by the script.
#' @param seed Optional integer if the script used randomness.
#' @param caption Optional named list with the FR-040..FR-041b
#'     caption-pipeline metadata: \code{source_path},
#'     \code{source_kind}, \code{panel_letter_count},
#'     \code{composite_panel_count}. Typically the
#'     \code{\link{compose_caption}} return value.
#'
#' @return Invisibly the sidecar list (same shape as written to disk).
#'
#' @examples
#' \dontrun{
#' con3 <- pg_connect_panel("dev3")
#' dep3 <- deployment_provenance(con3, "dev3")
#' write_sidecar(
#'     artifact_id   = "fig01-overview/panelB",
#'     artifact_path = "figures/fig01-overview/out/panelB.pdf",
#'     deployments   = list(dep3$deployment),
#'     manifests     = list(dep3$manifest),
#'     script_path   = "figures/fig01-overview/build.R",
#'     queries       = list()
#' )
#' }
#'
#' @importFrom jsonlite toJSON
#' @importFrom purrr map_chr
#' @importFrom logger log_info
#' @importFrom rlang abort
#' @export
write_sidecar <- function(artifact_id,
                          artifact_path,
                          deployments,
                          manifests,
                          script_path,
                          queries        = list(),
                          external_inputs = list(),
                          parameters     = list(),
                          seed           = NULL,
                          caption        = NULL) {

    if (length(deployments) == 0L) {
        rlang::abort(
            "Provenance sidecar requires at least one deployment"
        )
    }
    if (length(manifests) == 0L) {
        rlang::abort(
            "Provenance sidecar requires at least one manifest"
        )
    }

    package_commits <- list()
    for (m in manifests) {
        package_commits <- c(package_commits, m$packages)
    }
    package_commits <- package_commits[!duplicated(names(package_commits))]

    sidecar <- list(
        artifact_id     = artifact_id,
        produced_at     = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
        deployments     = lapply(deployments, sidecar_deployment_entry),
        package_commits = package_commits,
        script_path     = script_path,
        script_commit   = git_commit_of(script_path),
        queries         = queries,
        external_inputs = external_inputs
    )

    # `parameters` and `seed` are schema-optional; emit only when set so
    # jsonlite can't accidentally serialize an empty named list as the
    # JSON array `[]` (which would fail the schema's `type: object`).
    if (length(parameters) > 0L) sidecar$parameters <- parameters
    if (!is.null(seed))          sidecar$seed       <- as.integer(seed)
    if (!is.null(caption))       sidecar$caption    <- list(
        source_path           = caption$source_path,
        source_kind           = caption$source_kind,
        panel_letter_count    = as.integer(caption$panel_letter_count),
        composite_panel_count = as.integer(caption$composite_panel_count)
    )

    out_path <- paste0(artifact_path, ".provenance.json")
    writeLines(
        jsonlite::toJSON(
            sidecar, auto_unbox = TRUE, pretty = TRUE, null = "null"
        ),
        out_path
    )

    logger::log_info(
        "Wrote sidecar {out_path} ",
        "(deployments={paste(purrr::map_chr(deployments, 'name'), collapse=',')})"
    )
    invisible(sidecar)
}


#' Build a sidecar's deployments-array entry from a deployment record
#'
#' Picks the schema-required fields off a deployment record and adds
#' the default \code{build_manifest_source = "table"} when missing.
#'
#' @param entry A list with \code{name}, \code{url}, \code{db_host},
#'     \code{db_port}, \code{db_name}, \code{build_id}, optionally
#'     \code{build_manifest_source}.
#'
#' @return List with the schema-required fields.
#'
#' @keywords internal
#' @noRd
sidecar_deployment_entry <- function(entry) {

    list(
        name                  = entry$name,
        url                   = entry$url,
        db_host               = entry$db_host,
        db_port               = as.integer(entry$db_port),
        db_name               = entry$db_name,
        build_id              = entry$build_id,
        build_manifest_source = entry$build_manifest_source %||% "table"
    )
}


#' Resolve a deployment + its native build_manifest in one shot
#'
#' Convenience helper for \code{build.R} orchestrators: opens / reuses
#' a panel-pool connection to \code{deployment}, reads
#' \code{build_manifest} natively (cycle-001+), archives the manifest
#' under \code{manifests/}, and returns both the augmented deployment
#' record (ready for the sidecar's \code{deployments} array) and the
#' manifest itself (for \code{write_sidecar}'s parallel
#' \code{manifests} arg).
#'
#' @param deployment Character: deployment label
#'     (e.g. \code{"dev3"}, \code{"dev4"}).
#' @param archive_dir Character: directory for the archived
#'     manifest JSON / SHA256 files; defaults to \code{"manifests"}.
#'
#' @return A list with two elements:
#'   \itemize{
#'     \item \code{deployment} — deployment record augmented with
#'       \code{build_id} + \code{build_manifest_source} ready for the
#'       sidecar.
#'     \item \code{manifest} — the full Build Manifest as returned by
#'       \code{\link{build_manifest_for}}.
#'   }
#'
#' @examples
#' \dontrun{
#' dep3 <- deployment_provenance("dev3")
#' dep4 <- deployment_provenance("dev4")
#' write_sidecar(
#'     artifact_id = "fig01-overview",
#'     artifact_path = "figures/fig01-overview/out/fig01-overview.pdf",
#'     deployments = list(dep3$deployment, dep4$deployment),
#'     manifests   = list(dep3$manifest,   dep4$manifest),
#'     script_path = "figures/fig01-overview/build.R"
#' )
#' }
#'
#' @importFrom logger log_info
#' @export
deployment_provenance <- function(deployment, archive_dir = "manifests") {

    record <- load_connection(deployment = deployment)
    con <- pg_connect_panel(deployment)
    manifest <- build_manifest_for(con)
    write_manifest(manifest, deployment = deployment, dir = archive_dir)

    list(
        deployment = list(
            name                  = record$name,
            url                   = record$url,
            db_host               = record$db_host,
            db_port               = record$db_port,
            db_name               = record$db_name,
            build_id              = manifest$build_id,
            build_manifest_source = "table"
        ),
        manifest = manifest
    )
}


#' Build a query record from a pg_query tibble's attributes
#'
#' Helper to convert one of the tibbles returned by
#' \code{\link{pg_query}} (or \code{\link{pg_query_panel}}) into the
#' dict shape expected by \code{\link{write_sidecar}}'s
#' \code{queries} argument. When the input tibble carries a
#' \code{"deployment"} attribute (always set by
#' \code{\link{pg_query_panel}}), it is preserved on the record so
#' the orchestrator can audit which deployment each query touched.
#'
#' @param rows A tibble returned by \code{\link{pg_query}} or
#'     \code{\link{pg_query_panel}}.
#'
#' @return A named list with \code{sql}, \code{row_count},
#'     \code{result_hash}, and optionally \code{deployment}.
#'
#' @examples
#' \dontrun{
#' rows <- pg_query(con, "SELECT 1")
#' query_record(rows)
#' }
#'
#' @export
query_record <- function(rows) {

    rec <- list(
        sql         = attr(rows, "sql") %||% "",
        row_count   = nrow(rows),
        result_hash = attr(rows, "result_hash") %||% ""
    )
    deployment <- attr(rows, "deployment")
    if (!is.null(deployment)) rec$deployment <- deployment
    rec
}


#' Resolve the git commit hash for a path
#'
#' Calls \code{git log -1 --format=%h -- <path>}. Returns the literal
#' \code{"untracked"} when the file is not tracked, or
#' \code{"no-git"} when git is unavailable.
#'
#' @param path Character: file path.
#'
#' @return Character scalar.
#'
#' @keywords internal
#' @noRd
git_commit_of <- function(path) {

    if (Sys.which("git") == "") return("no-git")
    res <- suppressWarnings(system2(
        "git",
        c("log", "-1", "--format=%h", "--", path),
        stdout = TRUE, stderr = FALSE
    ))
    if (length(res) == 0L) "untracked" else res[[1L]]
}
