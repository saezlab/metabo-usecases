#' Locate the Panel A statistics schema bundled with the package
#'
#' Returns the path to the vendored
#' \code{panel-a-stats.schema.json}. The canonical source is the
#' spec-kit repository
#' (\code{saezverse/ai/specifications/omnipath-metabo-figures/}); the
#' file under \code{inst/extdata/contracts/} is a copy kept in sync
#' with that spec-kit revision. The two locations are checked in
#' order so a sibling spec-kit checkout is preferred during
#' development.
#'
#' @return Character.
#'
#' @keywords internal
#' @noRd
panel_a_stats_schema_path <- function() {

    spec_kit_candidate <- file.path(
        "..", "..",
        "saezverse", "ai", "specifications",
        "omnipath-metabo-figures",
        "specs", "001-figures-pipeline",
        "contracts", "panel-a-stats.schema.json"
    )
    if (file.exists(spec_kit_candidate)) {
        return(normalizePath(spec_kit_candidate))
    }

    vendored <- system.file(
        "extdata", "contracts", "panel-a-stats.schema.json",
        package = "metabo.figures",
        mustWork = FALSE
    )
    if (nzchar(vendored) && file.exists(vendored)) {
        return(vendored)
    }

    # Fallback for in-repo runs before the package is installed.
    in_repo <- file.path(
        "inst", "extdata", "contracts",
        "panel-a-stats.schema.json"
    )
    if (file.exists(in_repo)) {
        return(normalizePath(in_repo))
    }

    rlang::abort(paste0(
        "panel-a-stats.schema.json not found. Looked in: ",
        spec_kit_candidate, ", ", vendored, ", ", in_repo
    ))
}


#' Section metadata — fixed five-section ordering
#'
#' Maps each section id to the human-readable name the schema's
#' \code{sections[].name} enum requires.
#'
#' @return Named character vector keyed by section id.
#'
#' @keywords internal
#' @noRd
panel_a_section_names <- function() {
    c(
        `1` = "Entities",
        `2` = "Metabolite-Protein Interactions",
        `3` = "Interactions",
        `4` = "Structures",
        `5` = "Annotation"
    )
}


#' Build the Figure 1 Panel A statistics digest (FR-043 orchestrator)
#'
#' Top-level entry point. Loads the digest config, runs every
#' section, validates the FR-043b subset constraint (inside the
#' Section-2 builder), validates the FR-043h snapshot binding,
#' validates the final \code{stats.json} against
#' \code{contracts/panel-a-stats.schema.json}, and writes the
#' four output files plus the provenance sidecar.
#'
#' @param snapshot_id Character or \code{NULL}: the snapshot
#'     identifier Figure 1's composite is being built against. If
#'     non-\code{NULL}, FR-043h enforces that the value matches the
#'     digest's \code{metadata.snapshot_id}; mismatch aborts. If
#'     \code{NULL}, the digest reads the snapshot id from
#'     \code{dev5.build_manifest} and records it.
#' @param out_dir Character: output directory. Defaults to
#'     \code{figures/fig01-overview/panel-a-stats/} relative to the
#'     repo root.
#' @param config_path Character or \code{NULL}: digest-config.yaml
#'     path. \code{NULL} → default location next to \code{out_dir}.
#' @param caption_sty Character: shared caption layout used for the
#'     typeset PDF. Defaults to the in-repo
#'     \code{tex/caption.sty}.
#' @param validate_schema Logical: whether to validate the assembled
#'     \code{stats.json} against the vendored schema before writing.
#'     Default \code{TRUE}; pass \code{FALSE} only in tests against a
#'     deliberately malformed fixture.
#'
#' @return Named list:
#'     \itemize{
#'       \item \code{out_dir}: output directory.
#'       \item \code{files}: named list (\code{json, csv, md, pdf}).
#'       \item \code{sidecar}: sidecar path.
#'       \item \code{config_path}, \code{config_sha256}: config
#'             provenance.
#'       \item \code{snapshot_id}: the snapshot id recorded in the
#'             digest.
#'     }
#'
#' @importFrom dplyr bind_rows
#' @importFrom fs dir_create file_exists path
#' @importFrom jsonlite toJSON write_json
#' @importFrom logger log_info log_warn
#' @importFrom purrr map map_chr pluck
#' @importFrom readr write_csv
#' @importFrom rlang abort
#' @importFrom tibble tibble
#' @export
build_panel_a_digest <- function(
    snapshot_id = NULL,
    out_dir = "figures/fig01-overview/panel-a-stats",
    config_path = NULL,
    caption_sty = "tex/caption.sty",
    validate_schema = TRUE
) {

    setup_pipeline_log("stats:panel_a_digest")

    cfg <- digest_config(config_path)
    defs <- cfg$definitions
    rt   <- cfg$runtime

    fs::dir_create(out_dir, recurse = TRUE)

    logger::log_info(
        "build_panel_a_digest: out_dir={out_dir}, ",
        "config={cfg$path}"
    )

    # ─── Section metrics ─────────────────────────────────────────────
    entities    <- section_entities()
    mpi         <- section_mpi(definitions = defs)
    interactions <- section_interactions(
        definitions = defs,
        runtime     = rt
    )
    structures  <- section_structures(runtime = rt)
    annotation  <- section_annotation(runtime = rt)

    section_metrics <- list(
        `1` = entities,
        `2` = mpi,
        `3` = interactions,
        `4` = structures,
        `5` = annotation
    )

    # ─── Section-level contributing-resource enumeration ─────────────
    section_resources <- purrr::map(1:5, function(i) {
        section_contributing_resources(i)
    })
    names(section_resources) <- as.character(1:5)

    # ─── Chosen branches → metadata.definitions block ────────────────
    branches <- list(
        pathway          = attr(interactions, "branches")$pathway,
        reaction         = attr(interactions, "branches")$reaction,
        annotation_slot3 = attr(annotation,  "branches")$annotation_slot3
    )
    metadata_definitions <- mirror_to_metadata(defs, branches)

    # ─── Snapshot binding (FR-043h) ──────────────────────────────────
    deployment_name <- single_deployment_or_abort(section_metrics)
    snapshot_payload <- resolve_snapshot_id(
        snapshot_id,
        section_metrics,
        deployment_name
    )

    # ─── Assemble PanelAStatisticsDigest object ──────────────────────
    digest_obj <- list(
        metadata = list(
            schema_version  = "1.0.0",
            produced_at     = format(
                Sys.time(),
                "%Y-%m-%dT%H:%M:%S%z"
            ),
            deployments     = list(list(
                name        = deployment_name,
                section_ids = 1:5,
                build_id    = snapshot_payload$build_id
            )),
            snapshot_ids    = snapshot_payload$snapshot_ids,
            script_path     = "R/stats-panel_a_digest.R",
            script_commit   = git_commit_of("R/stats-panel_a_digest.R"),
            definitions     = metadata_definitions
        ),
        sections = purrr::map(1:5, function(i) {
            build_section_object(
                section_id = i,
                metrics    = section_metrics[[as.character(i)]],
                resources  = section_resources[[as.character(i)]]
            )
        })
    )

    # ─── Schema validation ───────────────────────────────────────────
    if (isTRUE(validate_schema)) {
        validate_digest_against_schema(digest_obj)
    }

    # ─── Output files ────────────────────────────────────────────────
    json_path <- fs::path(out_dir, "stats.json")
    csv_path  <- fs::path(out_dir, "stats.csv")
    md_path   <- fs::path(out_dir, "stats.md")
    pdf_path  <- fs::path(out_dir, "stats.pdf")
    sidecar_path <- fs::path(out_dir, "stats.provenance.json")

    write_digest_json(digest_obj, json_path)
    write_digest_csv(digest_obj, csv_path)
    write_digest_markdown(digest_obj, md_path)
    write_digest_pdf(md_path, pdf_path, caption_sty)
    write_digest_sidecar(
        digest_obj,
        section_metrics,
        section_resources,
        cfg,
        sidecar_path
    )

    logger::log_info(
        "Panel A digest emitted: {json_path}, {csv_path}, ",
        "{md_path}, {pdf_path}, {sidecar_path}"
    )

    list(
        out_dir       = out_dir,
        files         = list(
            json = json_path,
            csv  = csv_path,
            md   = md_path,
            pdf  = pdf_path
        ),
        sidecar       = sidecar_path,
        config_path   = cfg$path,
        config_sha256 = cfg$sha256,
        snapshot_id   = snapshot_payload$build_id
    )
}


#' Resolve and assert the snapshot id for FR-043h
#'
#' Reads the deployment's \code{build_manifest} row, picks the
#' \code{build_id}, and (when the caller passed an expected snapshot
#' id) asserts the values match.
#'
#' @param expected_snapshot_id Character or \code{NULL}.
#' @param section_metrics Named list of section tibbles (used to
#'     surface the recorded deployment name in error messages).
#' @param deployment_name Character.
#'
#' @return Named list \code{{build_id, snapshot_ids}}.
#'
#' @importFrom rlang abort
#' @keywords internal
#' @noRd
resolve_snapshot_id <- function(
    expected_snapshot_id,
    section_metrics,
    deployment_name
) {

    con <- pg_connect_panel(deployment_name)
    manifest <- build_manifest_for(con)
    actual_build_id <- manifest$build_id

    if (!is.null(expected_snapshot_id) &&
        !identical(actual_build_id, expected_snapshot_id)) {
        rlang::abort(sprintf(
            paste0(
                "FR-043h snapshot-id mismatch: digest was passed ",
                "snapshot_id='%s' but %s.build_manifest.build_id='%s'. ",
                "The composite and the digest MUST share a snapshot id."
            ),
            expected_snapshot_id,
            deployment_name,
            actual_build_id
        ))
    }

    snapshot_ids <- list()
    if (!is.null(manifest$build_kind)) {
        snapshot_ids[[manifest$build_kind]] <- actual_build_id
    } else {
        snapshot_ids$main <- actual_build_id
    }

    list(
        build_id     = actual_build_id,
        snapshot_ids = snapshot_ids
    )
}


#' Assert every section ran against a single deployment (FR-043h)
#'
#' Post-2026-06-14 dev5 integrated-build promotion the digest opens
#' one connection per rebuild; every metric should carry the same
#' \code{deployment} value. A mismatch indicates a stale per-panel
#' override and aborts the rebuild before any output is written.
#'
#' @param section_metrics Named list of section tibbles.
#'
#' @return Character: the single deployment name.
#'
#' @importFrom rlang abort
#' @keywords internal
#' @noRd
single_deployment_or_abort <- function(section_metrics) {

    deployments <- unique(unlist(purrr::map(
        section_metrics,
        function(m) unique(m$deployment)
    )))
    deployments <- deployments[!is.na(deployments)]

    if (length(deployments) != 1L) {
        rlang::abort(sprintf(
            paste0(
                "FR-043h: digest expected a single deployment ",
                "post-2026-06-14, found: %s. Check ",
                "R/utils-deployment.R::deployment_registry() for stale ",
                "per-panel overrides."
            ),
            paste(deployments, collapse = ", ")
        ))
    }

    deployments
}


#' Build the JSON section object for one section
#'
#' @param section_id Integer 1..5.
#' @param metrics Tibble of DigestMetric rows.
#' @param resources Tibble from
#'     \code{\link{section_contributing_resources}}.
#'
#' @return List shaped for \code{stats.json -> sections[i]}.
#'
#' @keywords internal
#' @noRd
build_section_object <- function(section_id, metrics, resources) {

    metric_objs <- purrr::map(seq_len(nrow(metrics)), function(i) {
        list(
            metric_name      = metrics$metric_name[[i]],
            value            = as.integer(metrics$value[[i]]),
            definition_label = metrics$definition_label[[i]],
            state            = metrics$state[[i]],
            deployment       = metrics$deployment[[i]],
            sql_hash         = metrics$sql_hash[[i]]
        )
    })

    list(
        section_id     = as.integer(section_id),
        name           = unname(panel_a_section_names()[
            as.character(section_id)
        ]),
        metrics        = metric_objs,
        resource_count = as.integer(nrow(resources)),
        resources      = as.character(resources$resource_name)
    )
}


#' Validate the assembled digest object against the JSON Schema
#'
#' Serialises the in-memory list to a transient JSON string, then
#' uses \code{jsonvalidate::json_validate} (Suggests). The schema is
#' located via \code{\link{panel_a_stats_schema_path}}.
#'
#' @param digest_obj Named list.
#'
#' @return Invisibly \code{TRUE}; aborts on failure.
#'
#' @importFrom jsonlite toJSON
#' @importFrom rlang abort
#' @keywords internal
#' @noRd
validate_digest_against_schema <- function(digest_obj) {

    if (!requireNamespace("jsonvalidate", quietly = TRUE)) {
        logger::log_warn(
            "jsonvalidate not installed; skipping FR-043g schema check"
        )
        return(invisible(FALSE))
    }

    schema_path <- panel_a_stats_schema_path()
    json <- jsonlite::toJSON(
        digest_obj,
        auto_unbox = TRUE,
        null = "null",
        pretty = FALSE
    )

    result <- jsonvalidate::json_validate(
        json,
        schema_path,
        engine  = "ajv",
        verbose = TRUE
    )

    if (!isTRUE(result)) {
        errs <- attr(result, "errors")
        rlang::abort(sprintf(
            "FR-043g schema validation FAILED:\n%s",
            paste(utils::capture.output(print(errs)), collapse = "\n")
        ))
    }

    invisible(TRUE)
}


#' Write \code{stats.json}
#'
#' @keywords internal
#' @noRd
write_digest_json <- function(digest_obj, path) {
    jsonlite::write_json(
        digest_obj,
        path,
        auto_unbox = TRUE,
        null       = "null",
        pretty     = TRUE
    )
}


#' Write \code{stats.csv} — one row per metric
#'
#' @keywords internal
#' @noRd
write_digest_csv <- function(digest_obj, path) {

    rows <- purrr::map(digest_obj$sections, function(sect) {
        purrr::map(sect$metrics, function(m) {
            tibble::tibble(
                section          = sect$name,
                metric_name      = m$metric_name,
                value            = m$value,
                definition_label = m$definition_label,
                state            = m$state,
                deployment       = m$deployment
            )
        })
    })

    flat <- dplyr::bind_rows(unlist(rows, recursive = FALSE))
    readr::write_csv(flat, path)
}


#' Write \code{stats.md} — human-readable Markdown summary
#'
#' @keywords internal
#' @noRd
write_digest_markdown <- function(digest_obj, path) {

    lines <- c(
        "# Figure 1 Panel A — basic database statistics",
        "",
        sprintf(
            "Snapshot id: `%s`. Deployment: `%s`. Generated: %s.",
            digest_obj$metadata$snapshot_ids[[1L]],
            digest_obj$metadata$deployments[[1L]]$name,
            digest_obj$metadata$produced_at
        ),
        ""
    )

    for (sect in digest_obj$sections) {
        lines <- c(
            lines,
            sprintf("## Section %d — %s", sect$section_id, sect$name),
            "",
            sprintf(
                "Contributing resources: **%d** (%s)",
                sect$resource_count,
                paste(sect$resources, collapse = ", ")
            ),
            "",
            "| Metric | Value | Definition | State |",
            "|--------|------:|------------|-------|"
        )
        for (m in sect$metrics) {
            lines <- c(lines, sprintf(
                "| `%s` | %d | `%s` | %s |",
                m$metric_name,
                m$value,
                m$definition_label,
                m$state
            ))
        }
        lines <- c(lines, "")
    }

    writeLines(lines, con = path, useBytes = TRUE)
}


#' Write \code{stats.pdf} — typeset summary via xelatex
#'
#' Reuses \code{tex/caption.sty} (FR-041) for visual consistency
#' with the figure caption PDFs.
#'
#' @keywords internal
#' @noRd
write_digest_pdf <- function(md_path, pdf_path, caption_sty) {

    if (!file.exists(caption_sty)) {
        logger::log_warn(paste0(
            "caption.sty not found at {caption_sty}; skipping ",
            "stats.pdf typesetting (FR-043g pdf output)"
        ))
        return(invisible(FALSE))
    }

    md_to_latex_minimal(md_path, pdf_path, caption_sty)
}


#' Minimal Markdown → LaTeX → PDF conversion using xelatex
#'
#' Limits Markdown to headings, paragraphs, code spans, tables, and
#' bold/italic emphasis — the subset used by
#' \code{\link{write_digest_markdown}}. Shells out to \code{xelatex}.
#'
#' @keywords internal
#' @noRd
md_to_latex_minimal <- function(md_path, pdf_path, caption_sty) {

    md <- readLines(md_path, warn = FALSE, encoding = "UTF-8")
    tex <- assemble_digest_tex(md, caption_sty)

    work <- tempfile(pattern = "panel-a-stats-", fileext = ".tex")
    writeLines(tex, work, useBytes = TRUE)

    file.copy(
        caption_sty,
        file.path(dirname(work), basename(caption_sty)),
        overwrite = TRUE
    )

    args <- c(
        "-interaction=nonstopmode",
        sprintf("-output-directory=%s", dirname(work)),
        work
    )
    status <- system2("xelatex", args, stdout = TRUE, stderr = TRUE)

    produced_pdf <- sub("\\.tex$", ".pdf", work)
    if (!file.exists(produced_pdf)) {
        logger::log_warn(
            "xelatex did not produce a PDF for the digest; ",
            "log tail:\n{paste(tail(status, 20), collapse = '\n')}"
        )
        return(invisible(FALSE))
    }

    file.copy(produced_pdf, pdf_path, overwrite = TRUE)
    invisible(TRUE)
}


#' Assemble the wrapper LaTeX for the digest PDF
#'
#' @keywords internal
#' @noRd
assemble_digest_tex <- function(md_lines, caption_sty) {

    body <- md_to_latex_body(md_lines)

    c(
        "\\documentclass[10pt,a4paper]{article}",
        "\\usepackage{geometry}",
        "\\geometry{margin=15mm}",
        "\\usepackage{fontspec}",
        "\\setmainfont{TeX Gyre Heros}",
        "\\setmonofont{TeX Gyre Cursor}",
        "\\usepackage{xcolor}",
        "\\usepackage{booktabs}",
        "\\usepackage{longtable}",
        sprintf("\\usepackage{%s}", tools::file_path_sans_ext(
            basename(caption_sty)
        )),
        "\\begin{document}",
        body,
        "\\end{document}"
    )
}


#' Minimal Markdown → LaTeX body conversion
#'
#' @keywords internal
#' @noRd
md_to_latex_body <- function(md_lines) {

    out <- character()
    in_table <- FALSE
    table_cols <- 0L

    for (line in md_lines) {

        if (grepl("^# ", line)) {
            in_table <- close_table_if_open(out, in_table)
            heading <- sub("^# ", "", line)
            out <- c(out, sprintf(
                "\\section*{%s}",
                latex_escape(heading)
            ))
            next
        }
        if (grepl("^## ", line)) {
            in_table <- close_table_if_open(out, in_table)
            heading <- sub("^## ", "", line)
            out <- c(out, sprintf(
                "\\subsection*{%s}",
                latex_escape(heading)
            ))
            next
        }
        if (grepl("^\\| ", line)) {
            if (grepl("^\\|[- :]+\\|", line)) {
                next   # header separator
            }
            cells <- strsplit(line, "\\|", fixed = FALSE)[[1L]]
            cells <- trimws(cells[nchar(cells) > 0L])
            if (!in_table) {
                table_cols <- length(cells)
                col_spec <- paste(rep("l", table_cols), collapse = "")
                out <- c(
                    out,
                    sprintf("\\begin{tabular}{%s}", col_spec),
                    "\\toprule"
                )
                in_table <- TRUE
            }
            out <- c(out, paste(
                purrr::map_chr(cells, md_inline_to_latex),
                collapse = " & "
            ), "\\\\")
            next
        }

        if (in_table) {
            out <- c(out, "\\bottomrule", "\\end{tabular}", "")
            in_table <- FALSE
        }

        if (nchar(trimws(line)) == 0L) {
            out <- c(out, "")
            next
        }
        out <- c(out, md_inline_to_latex(line), "")
    }

    if (in_table) {
        out <- c(out, "\\bottomrule", "\\end{tabular}", "")
    }

    paste(out, collapse = "\n")
}


#' Close an open table block (mutates \code{out} in caller via return)
#'
#' @keywords internal
#' @noRd
close_table_if_open <- function(out_caller, in_table) {
    # The caller manages out lines; here we just return the new state.
    !in_table   # placeholder; caller resets in_table to FALSE
}


#' Minimal inline Markdown → LaTeX
#'
#' @keywords internal
#' @noRd
md_inline_to_latex <- function(s) {
    s <- latex_escape(s)
    s <- gsub("`([^`]+)`", "\\\\texttt{\\1}", s)
    s <- gsub("\\*\\*([^*]+)\\*\\*", "\\\\textbf{\\1}", s)
    s <- gsub("\\*([^*]+)\\*", "\\\\textit{\\1}", s)
    s
}


#' LaTeX-escape a literal string
#'
#' @keywords internal
#' @noRd
latex_escape <- function(s) {
    s <- gsub("\\\\", "\\\\textbackslash{}", s)
    s <- gsub("([&%$#_{}])", "\\\\\\1", s)
    s <- gsub("~", "\\\\textasciitilde{}", s)
    s <- gsub("\\^", "\\\\textasciicircum{}", s)
    s
}


#' Write \code{stats.provenance.json} sidecar
#'
#' Carries each metric's verbatim SQL text, the SQL hash, the
#' deployment, and the result hash plus the digest-config.yaml
#' fingerprint and the cli arguments — everything a reader needs to
#' reproduce the digest without re-deriving the queries.
#'
#' @keywords internal
#' @noRd
write_digest_sidecar <- function(
    digest_obj,
    section_metrics,
    section_resources,
    cfg,
    path
) {

    queries_per_section <- purrr::map(names(section_metrics), function(k) {
        list(
            section_id  = as.integer(k),
            section_name = unname(panel_a_section_names()[k]),
            queries     = attr(section_metrics[[k]], "queries")
        )
    })

    resource_queries <- purrr::map(names(section_resources), function(k) {
        rows <- section_resources[[k]]
        list(
            section_id  = as.integer(k),
            sql         = attr(rows, "sql"),
            result_hash = attr(rows, "result_hash"),
            deployment  = attr(rows, "deployment"),
            row_count   = nrow(rows)
        )
    })

    sidecar <- list(
        artifact_id   = "fig01-overview/panel-a-stats",
        produced_at   = digest_obj$metadata$produced_at,
        deployments   = digest_obj$metadata$deployments,
        snapshot_ids  = digest_obj$metadata$snapshot_ids,
        script_path   = digest_obj$metadata$script_path,
        script_commit = digest_obj$metadata$script_commit,
        config = list(
            path   = cfg$path,
            sha256 = cfg$sha256
        ),
        section_queries  = queries_per_section,
        resource_queries = resource_queries
    )

    jsonlite::write_json(
        sidecar,
        path,
        auto_unbox = TRUE,
        null       = "null",
        pretty     = TRUE
    )
}
