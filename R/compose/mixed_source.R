#' Compose a mixed-source figure (R panels + TikZ + manual assets)
#'
#' Dispatches between two composition paths:
#' \itemize{
#'     \item \strong{PDF-level} (default): runs \code{xelatex} on a
#'         LaTeX template that uses \code{pdfpages} to lay out each
#'         pre-rendered PDF panel. The template path is required.
#'     \item \strong{SVG-level}: shells out to
#'         \code{python/compose/svg_assemble.py} when the spec
#'         requests SVG-level placement.
#' }
#'
#' Constitution Principle III: every step is scripted; no GUI editor
#' steps for pipeline composites.
#'
#' @param mode One of \code{"pdf"}, \code{"svg"}.
#' @param spec A list describing the composite:
#'     \itemize{
#'         \item \code{output}: path of the composite to write.
#'         \item For \code{mode = "pdf"}: \code{template} (path to a
#'             \code{.tex} file consuming \code{compose_helpers.sty}),
#'             \code{component} (log component name).
#'         \item For \code{mode = "svg"}: \code{page} +
#'             \code{panels} per \code{svg_assemble.py}.
#'     }
#' @param work_dir Working directory for intermediate LaTeX files;
#'     defaults to \code{dirname(spec$output)}.
#'
#' @return Invisibly the path of the composite.
#'
#' @importFrom rlang abort
#' @importFrom fs path path_dir path_file path_ext file_exists dir_create
#' @importFrom logger log_info
#' @export
compose_mixed_source <- function(mode, spec, work_dir = NULL) {

    out_path <- spec$output
    if (is.null(out_path)) {
        rlang::abort("compose spec must include `output`")
    }

    fs::dir_create(fs::path_dir(out_path))

    if (identical(mode, "pdf")) {
        compose_pdf(spec, work_dir %||% fs::path_dir(out_path))
    } else if (identical(mode, "svg")) {
        compose_svg(spec)
    } else {
        rlang::abort(sprintf(
            "Unknown compose mode: '%s' (expected 'pdf' or 'svg')", mode
        ))
    }

    logger::log_info("Composite ready at {out_path}")
    invisible(out_path)
}


#' Run xelatex on a pdfpages-based template
#'
#' Uses \code{lib/log.sh::log_xelatex_capture} so xelatex output lands
#' in the unified pipeline log under \code{[LaTeX][...]}.
#'
#' @keywords internal
#' @noRd
compose_pdf <- function(spec, work_dir) {

    if (is.null(spec$template)) {
        rlang::abort("PDF composite spec must include `template`")
    }
    if (!file.exists(spec$template)) {
        rlang::abort(sprintf(
            "Composite template not found: %s", spec$template
        ))
    }

    component <- spec$component %||% "compose:pdf"

    sh <- sprintf(
        paste0(
            "set -e; source lib/log.sh; ",
            "log_xelatex_capture '%s' xelatex ",
            "-interaction=nonstopmode ",
            "-output-directory='%s' '%s'"
        ),
        component,
        normalizePath(work_dir, mustWork = FALSE),
        normalizePath(spec$template, mustWork = TRUE)
    )

    rc <- system2("bash", c("-c", sh))
    if (rc != 0L) {
        rlang::abort(sprintf(
            "xelatex composite step failed (rc=%d) — see %s",
            rc, Sys.getenv("METABO_FIGURES_LOG")
        ))
    }

    # xelatex writes <template_basename>.pdf next to the .tex; move/copy
    # into the requested output path if they differ.
    template_pdf <- file.path(
        work_dir,
        sub("\\.tex$", ".pdf", basename(spec$template))
    )
    if (!identical(normalizePath(template_pdf, mustWork = FALSE),
                   normalizePath(spec$output,  mustWork = FALSE))) {
        file.copy(template_pdf, spec$output, overwrite = TRUE)
    }
}


#' Shell out to svg_assemble.py
#'
#' Passes the spec as YAML on stdin.
#'
#' @importFrom yaml as.yaml
#' @keywords internal
#' @noRd
compose_svg <- function(spec) {

    yaml_spec <- yaml::as.yaml(spec)

    python_bin <- file.path("python", ".venv", "bin", "python")
    if (!file.exists(python_bin)) python_bin <- "python3"

    component <- spec$component %||% "svg_assemble"

    rc <- system2(
        python_bin,
        c("-m", "python.compose.svg_assemble",
          "--component", component, "-"),
        input = yaml_spec
    )

    if (rc != 0L) {
        rlang::abort(sprintf(
            "svg_assemble.py failed (rc=%d) — see %s",
            rc, Sys.getenv("METABO_FIGURES_LOG")
        ))
    }
}
