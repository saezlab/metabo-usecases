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

    fs::dir_create(work_dir)
    xelatex_log <- file.path(work_dir, "compose_xelatex.log")

    # Copy the template + compose_helpers.sty into the work_dir and
    # invoke xelatex with that as the CWD, so panel paths in the
    # template are bare basenames resolved from the panel artifacts'
    # own directory (FR-005a layout: each figure folder has out/).
    template_basename <- basename(spec$template)
    file.copy(spec$template,
              file.path(work_dir, template_basename),
              overwrite = TRUE)
    sty_src <- "tex/compose_helpers.sty"
    if (file.exists(sty_src)) {
        file.copy(sty_src,
                  file.path(work_dir, "compose_helpers.sty"),
                  overwrite = TRUE)
    }

    # normalizePath() with mustWork = FALSE leaves missing files
    # un-resolved; spell out the absolute log path so it survives
    # the setwd() below.
    abs_workdir <- normalizePath(work_dir, mustWork = TRUE)
    abs_log     <- file.path(abs_workdir, basename(xelatex_log))
    old_wd <- setwd(abs_workdir)
    on.exit(setwd(old_wd), add = TRUE)

    logger::log_info("xelatex in {abs_workdir}")
    rc <- system2(
        "xelatex",
        c("-interaction=nonstopmode", template_basename),
        stdout = abs_log,
        stderr = abs_log
    )

    if (rc != 0L) {
        rlang::abort(sprintf(
            "xelatex composite step failed (rc=%d) — see %s",
            rc, xelatex_log
        ))
    }

    # xelatex writes <template_basename>.pdf in the work_dir; copy
    # into the requested output path if they differ.
    setwd(old_wd)
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
