#' Write a provenance sidecar for a pipeline artifact
#'
#' Serializes the run-level + per-query inputs into the JSON shape
#' documented in
#' \code{specs/001-figures-pipeline/contracts/provenance-sidecar.schema.json}.
#' Every pipeline figure / table / extra plot emits one sidecar
#' alongside its output file at \code{<artifact>.provenance.json}.
#'
#' @param artifact_id Character: stable id, e.g.
#'     \code{"fig01-overview/panelB"}.
#' @param artifact_path Character: path of the artifact whose sidecar
#'     this is -- the sidecar file is written to
#'     \code{paste0(artifact_path, ".provenance.json")}.
#' @param deployment Deployment record from
#'     \code{\link{load_connection}}.
#' @param manifests List of Build Manifests (one per build the artifact
#'     touched).
#' @param script_path Character: path to the generating script
#'     relative to the repo root.
#' @param queries List of query records -- each list element has
#'     \code{sql}, \code{row_count}, \code{result_hash}. Typically
#'     accumulated via the \code{attr(rows, "sql")} +
#'     \code{attr(rows, "result_hash")} produced by \code{\link{pg_query}}.
#' @param external_inputs Optional list of non-Postgres inputs (manual
#'     assets, palettes, vendored RData). Each entry has \code{kind},
#'     \code{path}, \code{fingerprint}, optional \code{source}.
#' @param parameters Optional named list of filter / contrast / threshold
#'     values used by the script.
#' @param seed Optional integer if the script used randomness.
#'
#' @return Invisibly the sidecar list (same shape as written to disk).
#'
#' @examples
#' \dontrun{
#' sc <- write_sidecar(
#'     artifact_id  = "fig01-overview/panelB",
#'     artifact_path = "figures/fig01-overview/out/panelB.pdf",
#'     deployment   = load_connection(),
#'     manifests    = list(main_manifest),
#'     script_path  = "figures/fig01-overview/build.R",
#'     queries      = list()
#' )
#' }
#'
#' @importFrom jsonlite toJSON
#' @importFrom purrr map_chr map
#' @importFrom logger log_info
#' @importFrom rlang abort
#' @export
write_sidecar <- function(artifact_id,
                          artifact_path,
                          deployment,
                          manifests,
                          script_path,
                          queries        = list(),
                          external_inputs = list(),
                          parameters     = list(),
                          seed           = NULL) {

    if (length(manifests) == 0L) {
        rlang::abort(
            "Provenance sidecar requires at least one manifest"
        )
    }

    snapshot_ids <- list()
    package_commits <- list()
    for (m in manifests) {
        snapshot_ids[[m$build]]   <- snapshot_id(m)
        package_commits           <- c(package_commits, m$packages)
    }
    package_commits <- package_commits[!duplicated(names(package_commits))]

    sidecar <- list(
        artifact_id      = artifact_id,
        produced_at      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
        deployment       = list(
            name    = deployment$name,
            url     = deployment$url,
            db_host = deployment$db_host,
            db_port = deployment$db_port,
            db_name = deployment$db_name
        ),
        snapshot_ids     = snapshot_ids,
        package_commits  = package_commits,
        script_path      = script_path,
        script_commit    = git_commit_of(script_path),
        queries          = queries,
        external_inputs  = external_inputs,
        parameters       = parameters
    )

    if (!is.null(seed)) sidecar$seed <- as.integer(seed)

    out_path <- paste0(artifact_path, ".provenance.json")
    writeLines(
        jsonlite::toJSON(
            sidecar, auto_unbox = TRUE, pretty = TRUE, null = "null"
        ),
        out_path
    )

    logger::log_info("Wrote sidecar {out_path}")
    invisible(sidecar)
}


#' Build a query record from a pg_query tibble's attributes
#'
#' Helper to convert one of the tibbles returned by
#' \code{\link{pg_query}} into the dict shape expected by
#' \code{\link{write_sidecar}}'s \code{queries} argument.
#'
#' @param rows A tibble returned by \code{\link{pg_query}}.
#'
#' @return A named list with \code{sql}, \code{row_count},
#'     \code{result_hash}.
#'
#' @examples
#' \dontrun{
#' rows <- pg_query(con, "SELECT 1")
#' query_record(rows)
#' }
#'
#' @export
query_record <- function(rows) {
    list(
        sql         = attr(rows, "sql") %||% "",
        row_count   = nrow(rows),
        result_hash = attr(rows, "result_hash") %||% ""
    )
}


#' Resolve the git commit hash for a path
#'
#' Calls \code{git log -1 --format=%h -- <path>}. Returns the literal
#' \code{"untracked"} when the file is not tracked, or
#' \code{"no-git"} when git is unavailable.
#'
#' @param path Character: file path.
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
