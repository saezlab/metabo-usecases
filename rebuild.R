#!/usr/bin/env Rscript

# rebuild.R
#
# Top-level rebuild entrypoint for the metabo.figures pipeline. See
# specs/001-figures-pipeline/contracts/rebuild-cli.md for the
# authoritative CLI surface.

suppressPackageStartupMessages({
    library(metabo.figures)
})

# ---- CLI parsing -----------------------------------------------------------

parse_args <- function(argv) {

    out <- list(
        targets       = character(0),
        deployment    = NA_character_,
        snapshot_id   = NA_character_,
        png           = FALSE,
        bundle        = FALSE,
        dry_run       = FALSE,
        check         = FALSE,
        jobs          = 1L,
        verbose       = FALSE,
        quiet         = FALSE
    )

    i <- 1L
    while (i <= length(argv)) {
        a <- argv[[i]]
        if (a %in% c("-h", "--help")) {
            cat(usage(), sep = "\n")
            quit(save = "no", status = 0L)
        } else if (a == "--deployment") {
            out$deployment <- argv[[i + 1L]]; i <- i + 2L
        } else if (a == "--snapshot-id") {
            out$snapshot_id <- argv[[i + 1L]]; i <- i + 2L
        } else if (a == "--png") {
            out$png <- TRUE; i <- i + 1L
        } else if (a == "--no-png") {
            out$png <- FALSE; i <- i + 1L
        } else if (a == "--bundle") {
            out$bundle <- TRUE; i <- i + 1L
        } else if (a == "--dry-run") {
            out$dry_run <- TRUE; i <- i + 1L
        } else if (a == "--check") {
            out$check <- TRUE; i <- i + 1L
        } else if (a == "--jobs") {
            out$jobs <- as.integer(argv[[i + 1L]]); i <- i + 2L
        } else if (a %in% c("-v", "--verbose")) {
            out$verbose <- TRUE; i <- i + 1L
        } else if (a %in% c("-q", "--quiet")) {
            out$quiet <- TRUE; i <- i + 1L
        } else if (startsWith(a, "--")) {
            stop(sprintf("Unknown flag: %s (try --help)", a), call. = FALSE)
        } else {
            out$targets <- c(out$targets, a)
            i <- i + 1L
        }
    }

    out
}


usage <- function() c(
    "Usage: ./rebuild.sh [OPTIONS] [TARGET ...]",
    "",
    "Targets:",
    "  fig01-overview                 a single figure folder by id",
    "  figures                        every figure under figures/",
    "  tables                         every table under tables/",
    "  main                           every main-text artifact",
    "  supplementary                  every supplementary artifact",
    "  fig01-overview/panelB          a single panel by panel-id",
    "  fig01-overview:depends         the artifact plus its dependents",
    "",
    "Options:",
    "  --deployment NAME              override the beauty deployment",
    "  --snapshot-id HEX12            pin an existing Snapshot Identifier",
    "  --no-png / --png               PNG companion (FR-027)",
    "  --bundle                       assemble the concatenated PDF bundle",
    "  --dry-run                      print work plan; do not execute",
    "  --check                        run gates after the rebuild",
    "  --jobs N                       parallel script-execution count",
    "  -v / --verbose                 stream per-script logs to stdout",
    "  -q / --quiet                   suppress per-script logs",
    "",
    "See specs/001-figures-pipeline/contracts/rebuild-cli.md."
)


# ---- log lifecycle ---------------------------------------------------------

setup_run_log <- function(args, snapshot_label = "pending") {

    log_dir <- file.path(getwd(), "logs")
    dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
    ts <- format(Sys.time(), "%Y%m%dT%H%M%S")
    log_path <- file.path(
        log_dir,
        sprintf("rebuild-%s-%s.log", snapshot_label, ts)
    )

    Sys.setenv(METABO_FIGURES_LOG = log_path)

    latest <- file.path(log_dir, "latest.log")
    suppressWarnings(file.remove(latest))
    suppressWarnings(file.symlink(basename(log_path), latest))

    setup_pipeline_log("rebuild")

    dep_label <- if (is.na(args$deployment)) "(default)" else args$deployment
    logger::log_info(
        "Starting rebuild (deployment={dep_label} jobs={args$jobs})"
    )
    if (length(args$targets) > 0L) {
        logger::log_info(
            "Targets: {paste(args$targets, collapse = ', ')}"
        )
    }

    log_path
}


# ---- work-plan discovery ---------------------------------------------------

discover_targets <- function(targets) {

    builds <- c(
        list.files("figures", pattern = "build\\.R$", full.names = TRUE,
                   recursive = TRUE),
        list.files("tables",  pattern = "build\\.R$", full.names = TRUE,
                   recursive = TRUE)
    )
    builds <- sort(builds)

    if (length(targets) == 0L) return(builds)

    selected <- character(0)
    for (t in targets) {
        sub <- if (t == "figures") {
            builds[grepl("^figures/",  builds)]
        } else if (t == "tables") {
            builds[grepl("^tables/",   builds)]
        } else if (t == "main" || t == "supplementary") {
            classify_builds(builds, t)
        } else {
            builds[grepl(t, builds, fixed = TRUE)]
        }
        selected <- c(selected, sub)
    }

    unique(selected)
}


classify_builds <- function(builds, kind) {

    keep <- vapply(builds, function(b) {
        readme <- file.path(dirname(b), "README.md")
        if (!file.exists(readme)) return(TRUE)  # default: main
        content <- paste(readLines(readme), collapse = "\n")
        if (kind == "supplementary") {
            grepl("supplementary", content, ignore.case = TRUE)
        } else {
            !grepl("supplementary", content, ignore.case = TRUE)
        }
    }, logical(1L))

    builds[keep]
}


# ---- main ------------------------------------------------------------------

main <- function() {

    argv <- commandArgs(trailingOnly = TRUE)
    args <- parse_args(argv)

    log_path <- setup_run_log(args)

    builds <- discover_targets(args$targets)
    if (length(builds) == 0L) {
        logger::log_warn("No build.R scripts matched the targets")
        return(invisible(0L))
    }

    if (isTRUE(args$dry_run)) {
        for (b in builds) {
            cat(sprintf("[plan] %s\n", sub("/build\\.R$", "", b)))
        }
        return(invisible(0L))
    }

    started <- Sys.time()
    errors <- 0L

    for (b in builds) {
        artifact_id <- sub("/build\\.R$", "", b)
        logger::log_info("→ {artifact_id}")
        rc <- tryCatch({
            source(
                normalizePath(b, mustWork = TRUE),
                local = new.env(parent = globalenv())
            )
            0L
        }, error = function(e) {
            # Use paste0 not glue: error messages may contain {…}
            # tokens that confuse the default glue layout (e.g. SQL,
            # nested error formatters).
            logger::log_error(paste0(
                artifact_id, " failed: ", conditionMessage(e)
            ))
            1L
        })
        if (rc != 0L) errors <- errors + 1L
    }

    if (isTRUE(args$bundle) && errors == 0L) {
        pdfs <- list.files(
            c("figures", "tables"),
            pattern = "\\.pdf$",
            full.names = TRUE, recursive = TRUE
        )
        if (length(pdfs) > 0L) {
            assemble_bundle(sort(pdfs), output = "out/manuscript-bundle.pdf")
        }
    }

    elapsed <- difftime(Sys.time(), started, units = "secs")
    logger::log_info(sprintf(
        "Done. %d builds, %d errors. (%.1fs)",
        length(builds), errors, as.numeric(elapsed)
    ))

    invisible(if (errors == 0L) 0L else 2L)
}


# Allow the file to be sourced for testing without auto-executing.
if (sys.nframe() == 0L) {
    rc <- main()
    quit(save = "no", status = as.integer(rc))
}
