#' Compile a \code{gt} object to a standalone LaTeX PDF + write the CSV
#'
#' Renders a \code{gt_tbl} via its LaTeX backend (\code{gt::as_latex})
#' wrapped in a minimal \code{standalone}-class document and compiled
#' with \code{xelatex}. The LaTeX path is preferred over the headless
#' Chrome path the default \code{gt::gtsave} chooses for PDF: Chrome
#' mis-measures rotated column labels (e.g. the 45° resource labels
#' used by the FR-015a checkmark table), and xelatex is already the
#' caption-pipeline dependency for figures so this introduces no new
#' tooling. The CSV is the underlying \code{data} tibble, written via
#' \code{readr::write_csv} for FR-026 (every table emits a typeset
#' PDF AND a machine-readable CSV).
#'
#' @param gt_obj A \code{gt_tbl} (the result of \code{\link[gt]{gt}}
#'     and friends).
#' @param data Tibble: the underlying wide-format data the table is
#'     built on. Written verbatim as the CSV companion.
#' @param out_dir Character: directory the artifacts are written to.
#'     Created if it does not exist.
#' @param slug Character: artifact slug — outputs are
#'     \code{<out_dir>/<slug>.pdf} and \code{<out_dir>/<slug>.csv}.
#'
#' @return A named list with elements \code{pdf} and \code{csv}
#'     pointing at the produced files. Aborts (via \code{rlang::abort})
#'     if xelatex fails or the expected PDF is not produced.
#'
#' @examples
#' \dontrun{
#' library(gt)
#' gt_obj <- gt(head(mtcars))
#' tables_save_pdf_csv(
#'     gt_obj  = gt_obj,
#'     data    = head(mtcars),
#'     out_dir = "tables/tabXX-example/out",
#'     slug    = "tabXX-example"
#' )
#' }
#'
#' @param max_width_mm Numeric: maximum natural width for the table in
#'     millimetres. The standalone wrapper uses
#'     \code{varwidth=<max_width_mm>mm} so the table sits inside its
#'     own bounding box rather than spreading across an oversized
#'     canvas — this is what controls "huge empty space between
#'     columns" in the rendered PDF.
#'
#' @importFrom fs dir_create
#' @importFrom readr write_csv
#' @importFrom gt as_latex
#' @importFrom logger log_info
#' @importFrom rlang abort
#' @export
tables_save_pdf_csv <- function(gt_obj, data, out_dir, slug,
                                max_width_mm = 180) {

    tex_body <- as.character(gt::as_latex(gt_obj))
    tables_save_latex_pdf_csv(
        tex_body     = tex_body,
        data         = data,
        out_dir      = out_dir,
        slug         = slug,
        max_width_mm = max_width_mm
    )
}


#' Compile a raw LaTeX table body to PDF + write the CSV
#'
#' Lower-level companion to \code{\link{tables_save_pdf_csv}} for
#' renderers that build their own LaTeX (e.g. the FR-015a checkmark
#' table needs \code{\\rotatebox{45}{...}} column labels that the
#' \code{gt} LaTeX backend escapes). The LaTeX body MUST be a
#' self-contained fragment renderable inside the standard table
#' wrapper from \code{\link{tables_latex_wrapper}}.
#'
#' @param tex_body Character scalar: the LaTeX body to compile.
#' @param data Tibble: the wide-format data the table is built on
#'     (written verbatim as the CSV companion).
#' @param out_dir Character: directory the artifacts are written to.
#' @param slug Character: artifact slug — outputs are
#'     \code{<out_dir>/<slug>.pdf} and \code{<out_dir>/<slug>.csv}.
#'
#' @return A named list with elements \code{pdf} and \code{csv}.
#'
#' @param max_width_mm Numeric: passed through to
#'     \code{\link{tables_latex_wrapper}} as the
#'     \code{varwidth=<n>mm} standalone option.
#'
#' @importFrom fs dir_create
#' @importFrom readr write_csv
#' @importFrom logger log_info
#' @importFrom rlang abort
#' @export
tables_save_latex_pdf_csv <- function(tex_body, data, out_dir, slug,
                                      max_width_mm = 180) {

    fs::dir_create(out_dir)

    csv_path <- file.path(out_dir, sprintf("%s.csv", slug))
    readr::write_csv(data, csv_path)

    wrapper_name <- sprintf("%s.tex", slug)
    wrapper_path <- file.path(out_dir, wrapper_name)
    writeLines(
        tables_latex_wrapper(tex_body, max_width_mm = max_width_mm),
        wrapper_path
    )

    logger::log_info(
        "Compiling table {slug} to PDF via xelatex"
    )

    log_path <- file.path(out_dir, sprintf("%s_xelatex.log", slug))
    old_wd <- setwd(out_dir)
    on.exit(setwd(old_wd), add = TRUE)
    captured <- suppressWarnings(system2(
        "xelatex",
        c("-interaction=nonstopmode", wrapper_name),
        stdout = TRUE,
        stderr = TRUE
    ))
    setwd(old_wd)
    rc <- as.integer(attr(captured, "status") %||% 0L)
    writeLines(captured, log_path)

    if (rc != 0L) {
        rlang::abort(sprintf(paste0(
            "xelatex table-compile failed (rc=%d) for %s — see %s"
        ), rc, slug, log_path))
    }

    pdf_path <- file.path(out_dir, sprintf("%s.pdf", slug))
    if (!file.exists(pdf_path)) {
        rlang::abort(sprintf(
            "Expected table PDF was not produced: %s", pdf_path
        ))
    }

    logger::log_info(
        "Wrote {pdf_path} + {csv_path}"
    )
    list(pdf = pdf_path, csv = csv_path)
}


#' Build the standalone LaTeX wrapper for a gt-emitted table body
#'
#' \code{gt::as_latex} returns a \code{\\begin{longtable}...\\end{longtable}}
#' fragment plus the package list it depends on; wrapping it in
#' \code{standalone} with \code{longtable}-friendly geometry produces
#' a self-contained PDF cropped to the table's bounding box.
#'
#' The geometry deliberately uses a wide paper size so longtable does
#' not break across pages when xelatex computes the bounding box.
#' Per-table page sizing happens in the caption-and-table wrapper, not
#' here.
#'
#' @param body Character scalar: the LaTeX source emitted by
#'     \code{\link[gt]{as_latex}}.
#'
#' @return Character vector of LaTeX lines.
#'
#' @keywords internal
#' @noRd
tables_latex_wrapper <- function(body, max_width_mm = 180) {
    c(
        "% Auto-generated by R/tables-render.R — do not edit.",
        sprintf(
            "\\documentclass[border=4mm,varwidth=%dmm]{standalone}",
            as.integer(max_width_mm)
        ),
        "\\usepackage{amsmath}",
        # amssymb supplies \checkmark for the FR-015a coverage table.
        "\\usepackage{amssymb}",
        "\\usepackage{booktabs}",
        "\\usepackage{longtable}",
        "\\usepackage{array}",
        "\\usepackage{colortbl}",
        "\\usepackage{multirow}",
        "\\usepackage{graphicx}",
        "\\usepackage{xcolor}",
        "\\usepackage{anyfontsize}",
        # Helvetica Neue LT Std (Condensed family is what fc-cache
        # picks up from ~/.local/share/fonts on beauty since only the
        # Cn / BdCn / LtCn / MdCn variants were copied over). Use the
        # canonical family name so fontconfig auto-matches Bold/Italic.
        "\\usepackage{fontspec}",
        "\\setmainfont{Helvetica Neue LT Std}",
        # fontspec's \setmainfont changes \rmfamily but leaves
        # \sffamily / \ttfamily on the default (Computer Modern Sans),
        # so any \sffamily in the hand-built LaTeX (e.g. the FR-015a
        # squared-board table's column / row labels) would render in
        # the wrong font. Pin \sffamily to the same family so the
        # whole table renders consistently.
        "\\setsansfont{Helvetica Neue LT Std}",
        # Alternating row backgrounds (light gray / white) for the
        # gt-rendered tables. gt's LaTeX backend doesn't emit
        # \\rowcolor — we apply it in the standalone wrapper so it
        # affects every row in the body uniformly.
        "\\definecolor{tablerowalt}{HTML}{F2F2F2}",
        "\\rowcolors{2}{tablerowalt}{white}",
        "\\begin{document}",
        body,
        "\\end{document}"
    )
}


#' Compile a table-and-caption PDF + emit \code{caption.txt}
#'
#' Thin shim around \code{\link{compose_caption}} that adapts it to
#' tables (FR-040..FR-041a applies to tables too per FR-042). The
#' caption.tex source describes the table layout; the typeset PDF
#' places the table-PDF above the caption in the same standard
#' manuscript layout figures use.
#'
#' \code{panel_count} is the number of \code{(a)/(b)/...} sub-letters
#' declared in \code{caption.tex}. For tables this is typically 1
#' (one main panel = the table itself); when a table is composed of
#' multiple sub-tables, pass the appropriate count.
#'
#' @param table_id Character: artifact id (e.g.
#'     \code{"tab01-id-resolving"}).
#' @param table_pdf Character: path to the bare table PDF the caption
#'     describes.
#' @param caption_source Character: path to \code{caption.tex} in the
#'     table folder (FR-040a).
#' @param out_dir Character: directory the caption artifacts are
#'     written into.
#' @param panel_count Integer: number of caption sub-letters declared
#'     in the caption source.
#'
#' @return A named list with elements \code{with_caption_pdf},
#'     \code{caption_txt}, \code{panel_letter_count},
#'     \code{composite_panel_count}, \code{caption_source} — the
#'     same shape \code{\link{compose_caption}} returns.
#'
#' @examples
#' \dontrun{
#' tables_compose_caption(
#'     table_id       = "tab01-id-resolving",
#'     table_pdf      = "tables/tab01-id-resolving/out/tab01-id-resolving.pdf",
#'     caption_source = "tables/tab01-id-resolving/caption.tex",
#'     out_dir        = "tables/tab01-id-resolving/out",
#'     panel_count    = 1L
#' )
#' }
#'
#' @export
tables_compose_caption <- function(table_id,
                                   table_pdf,
                                   caption_source,
                                   out_dir,
                                   panel_count,
                                   body_width_mm = 180) {
    compose_caption(
        figure_id          = table_id,
        composite_pdf      = table_pdf,
        caption_source     = caption_source,
        out_dir            = out_dir,
        panel_count        = panel_count,
        composite_basename = sprintf("%s.pdf", table_id),
        caption_position   = "above",
        body_width_mm      = body_width_mm
    )
}
