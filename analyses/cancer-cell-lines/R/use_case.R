# Helpers shared by the Rmds of the cancer cell lines use case.

#' Read the use-case configuration
#'
#' @param file Path of the YAML file, relative to the repository root.
#' @return The configuration as a nested list.
load_config <- function(file = "analyses/cancer-cell-lines/config/config.yml") {
    yaml::read_yaml(here::here(file))
}

#' Directory for one kind of file, created if missing
#'
#' @param config Configuration from `load_config()`.
#' @param kind One of the entries of `config$paths` ("raw", "derived", "output").
#' @param ... Further path components below that directory.
#' @return Absolute path.
uc_dir <- function(config, kind, ...) {
    path <- here::here(config$paths[[kind]], ...)
    if (kind != "raw") dir.create(path, recursive = TRUE, showWarnings = FALSE)
    path
}

#' Record a filtering step
#'
#' Every step that removes samples or features goes through this function,
#' so counts before and after end up in one table per Rmd.
#'
#' @param log Data frame of earlier steps (or NULL).
#' @param step Short name of the step.
#' @param what "injections", "lines", "ions", ...
#' @param before,after Counts before and after the step.
#' @param reason Why the step is done.
#' @return `log` with one row appended.
log_filter <- function(log, step, what, before, after, reason) {
    rbind(
        log,
        data.frame(step = step, what = what, before = before, after = after,
                   removed = before - after, reason = reason)
    )
}

#' Write a manifest describing the files of a processing step
#'
#' @param files Named character vector of written files.
#' @param inputs Named character vector of input files.
#' @param path Where to write the JSON manifest.
#' @param extra List of further fields (e.g. parameters).
write_manifest <- function(files, inputs, path, extra = list()) {
    md5 <- function(f) unname(tools::md5sum(f))
    git_sha <- tryCatch(
        system2("git", c("-C", here::here(), "rev-parse", "HEAD"), stdout = TRUE),
        error = function(e) NA_character_
    )
    pkgs <- c("MetaProViz", "OmnipathR", "dplyr", "readxl", "rmarkdown")
    manifest <- c(
        list(
            created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
            git_commit = git_sha,
            r_version = R.version.string,
            packages = lapply(setNames(pkgs, pkgs), function(p) {
                if (requireNamespace(p, quietly = TRUE)) {
                    d <- utils::packageDescription(p)
                    list(version = d$Version, remote_sha = d$RemoteSha)
                } else NULL
            }),
            inputs = lapply(inputs, function(f) list(file = basename(f), md5 = md5(f))),
            outputs = lapply(files, function(f) list(file = basename(f), md5 = md5(f)))
        ),
        extra
    )
    jsonlite::write_json(manifest, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
    invisible(path)
}
