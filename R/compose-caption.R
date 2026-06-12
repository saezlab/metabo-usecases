#' Compile the caption-and-figure PDF and the plain-text caption
#'
#' Reads the per-figure caption source (\code{caption.tex} —
#' preferred — or \code{caption.md}), validates that the sub-letter
#' sequence matches the composite's panel count (FR-041b), generates
#' a self-contained wrapper \code{.tex} that places the composite
#' above the caption, shells out to \code{xelatex} (FR-041) to
#' produce \code{<slug>-with-caption.pdf}, and deterministically
#' strips the LaTeX markup to produce \code{caption.txt} (FR-041a).
#'
#' @param figure_id Character: figure slug, e.g.
#'     \code{"fig01-overview"}. Used to name the output
#'     \code{<figure_id>-with-caption.pdf}.
#' @param composite_pdf Character: path (relative to repo root) to
#'     the bare composite PDF the caption describes.
#' @param caption_source Character: path to \code{caption.tex} (or
#'     \code{caption.md}) — typically
#'     \code{figures/<figure_id>/caption.tex}.
#' @param out_dir Character: directory where the with-caption PDF
#'     and \code{caption.txt} land — typically
#'     \code{figures/<figure_id>/out}.
#' @param panel_count Integer: number of panels in the composite. The
#'     caption's \code{(a)/(b)/...} sub-letter count MUST equal this;
#'     a mismatch aborts with the FR-041b message.
#' @param caption_sty Character: path to the shared
#'     \code{tex/caption.sty} layout file. Copied next to the
#'     wrapper before xelatex runs so it is discoverable on the
#'     CWD's TEXINPUTS path.
#' @param composite_basename Character or \code{NULL}: filename the
#'     wrapper LaTeX uses to refer to the composite from inside
#'     \code{out_dir}. \code{NULL} (default) → use
#'     \code{basename(composite_pdf)}.
#'
#' @return A named list with elements \code{with_caption_pdf},
#'     \code{caption_txt}, \code{panel_letter_count},
#'     \code{composite_panel_count}, \code{caption_source}. The
#'     latter four are the values written into the sidecar's
#'     \code{caption} block (FR-040..FR-041b).
#'
#' @examples
#' \dontrun{
#' compose_caption(
#'     figure_id      = "fig01-overview",
#'     composite_pdf  = "figures/fig01-overview/out/fig01-overview.pdf",
#'     caption_source = "figures/fig01-overview/caption.tex",
#'     out_dir        = "figures/fig01-overview/out",
#'     panel_count    = 6L
#' )
#' }
#'
#' @importFrom logger log_info
#' @importFrom rlang abort
#' @importFrom stringr str_match_all str_squish str_replace_all
#' @importFrom fs path file_copy
#' @export
compose_caption <- function(figure_id,
                            composite_pdf,
                            caption_source,
                            out_dir,
                            panel_count,
                            caption_sty = "tex/caption.sty",
                            composite_basename = NULL,
                            caption_position = c("below", "above"),
                            body_width_mm = 180) {

    caption_position <- match.arg(caption_position)

    if (!file.exists(caption_source)) {
        rlang::abort(sprintf(
            "Caption source missing: %s. Author it per FR-040a.",
            caption_source
        ))
    }
    if (!file.exists(composite_pdf)) {
        rlang::abort(sprintf(
            "Composite PDF missing for caption-compile: %s",
            composite_pdf
        ))
    }
    if (!file.exists(caption_sty)) {
        rlang::abort(sprintf(
            "Caption style missing: %s", caption_sty
        ))
    }

    src_text <- paste(readLines(caption_source), collapse = "\n")
    panel_letters <- parse_panel_letters(src_text)

    if (length(panel_letters) != as.integer(panel_count)) {
        rlang::abort(sprintf(paste0(
            "Caption / composite panel-count mismatch (FR-041b) for %s: ",
            "caption.tex declares %d sub-letters (%s), composite has %d panels."
        ),
        figure_id,
        length(panel_letters),
        paste(panel_letters, collapse = ", "),
        as.integer(panel_count)
        ))
    }
    if (!identical(panel_letters, letters[seq_along(panel_letters)])) {
        rlang::abort(sprintf(paste0(
            "Caption sub-letter order does not match composite (FR-041b) ",
            "for %s: got %s, expected %s."
        ),
        figure_id,
        paste(panel_letters, collapse = ", "),
        paste(letters[seq_along(panel_letters)], collapse = ", ")
        ))
    }

    composite_name <-
        composite_basename %||% basename(composite_pdf)
    composite_dest <- file.path(out_dir, composite_name)
    if (!identical(normalizePath(composite_pdf, mustWork = FALSE),
                   normalizePath(composite_dest, mustWork = FALSE))) {
        file.copy(composite_pdf, composite_dest, overwrite = TRUE)
    }

    file.copy(caption_sty, file.path(out_dir, "caption.sty"),
              overwrite = TRUE)
    file.copy(caption_source, file.path(out_dir, "caption.tex"),
              overwrite = TRUE)

    wrapper_name <- sprintf("%s-with-caption.tex", figure_id)
    wrapper_path <- file.path(out_dir, wrapper_name)
    writeLines(
        caption_wrapper_tex(
            composite_name,
            caption_position,
            body_width_mm = body_width_mm
        ),
        wrapper_path
    )

    if (nchar(Sys.which("xelatex")) == 0L) {
        logger::log_warn(
            "xelatex not found — skipping caption PDF for {figure_id}"
        )
        return(invisible(list(
            figure_id             = figure_id,
            caption_source        = caption_source,
            panel_letters         = panel_letters,
            composite_panel_count = as.integer(panel_count),
            with_caption_pdf      = NULL,
            xelatex_skipped       = TRUE
        )))
    }

    logger::log_info(
        "Compiling caption-and-figure PDF for {figure_id} via xelatex"
    )
    log_path <- file.path(out_dir, "caption_xelatex.log")
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
            "xelatex caption-compile step failed (rc=%d) — see %s"
        ), rc, log_path))
    }

    with_caption_pdf <- file.path(
        out_dir, sprintf("%s-with-caption.pdf", figure_id)
    )
    if (!file.exists(with_caption_pdf)) {
        rlang::abort(sprintf(
            "Expected caption PDF was not produced: %s", with_caption_pdf
        ))
    }

    caption_txt_path <- file.path(out_dir, "caption.txt")
    writeLines(strip_latex(src_text), caption_txt_path)
    logger::log_info(
        "Wrote {with_caption_pdf} + {caption_txt_path}"
    )

    list(
        with_caption_pdf      = with_caption_pdf,
        caption_txt           = caption_txt_path,
        panel_letter_count    = length(panel_letters),
        composite_panel_count = as.integer(panel_count),
        caption_source        = caption_source
    )
}


#' Extract the (a)/(b)/... sub-letter sequence from a caption source
#'
#' Returns the lower-case letter sequence in the order it appears
#' (e.g. \code{c("a","b","c","d","e","f")}). Matches both the
#' \code{\\figpanel\{a\}\{...\}} macro form and the plain
#' \code{(a) text} prose form (so \code{caption.md} authors may use
#' the latter).
#'
#' @param text Character: caption source text.
#'
#' @return Character vector of letters in source-order.
#'
#' @keywords internal
#' @noRd
parse_panel_letters <- function(text) {

    macro_hits <- stringr::str_match_all(
        text, "\\\\figpanel\\{([a-z])\\}"
    )[[1L]][, 2L]
    if (length(macro_hits) > 0L) return(macro_hits)

    bold_hits <- stringr::str_match_all(
        text, "\\\\textbf\\{\\(([a-z])\\)\\}"
    )[[1L]][, 2L]
    if (length(bold_hits) > 0L) return(bold_hits)

    prose_hits <- stringr::str_match_all(
        text, "(?:^|\\n|\\s)\\(([a-z])\\)\\s"
    )[[1L]][, 2L]
    prose_hits
}


#' Build the wrapper LaTeX that composes the figure + caption PDF
#'
#' Self-contained standalone document; xelatex auto-crops to the
#' content's bounding box.
#'
#' @param composite_basename Character: filename of the composite
#'     PDF, relative to the xelatex CWD.
#'
#' @return Character vector of LaTeX lines.
#'
#' @keywords internal
#' @noRd
caption_wrapper_tex <- function(composite_basename,
                                caption_position = c("below", "above"),
                                body_width_mm = 180) {

    caption_position <- match.arg(caption_position)
    width_str <- sprintf("%dmm", as.integer(body_width_mm))

    composite_block <- sprintf(
        "  \\includegraphics[width=%s]{%s}\\par",
        width_str, composite_basename
    )
    caption_block <- c(
        "  \\justifying",
        "  \\input{caption.tex}"
    )

    if (identical(caption_position, "above")) {
        body <- c(caption_block, "  \\vspace{2mm}", composite_block)
    } else {
        body <- c(composite_block, "  \\vspace{2mm}", caption_block)
    }

    c(
        "% Auto-generated by R/compose-caption.R — do not edit.",
        "\\documentclass[border=4mm]{standalone}",
        "\\usepackage{graphicx}",
        "\\usepackage{caption}",
        "\\begin{document}",
        sprintf("\\begin{minipage}{%s}", width_str),
        body,
        "\\end{minipage}",
        "\\end{document}"
    )
}


#' Deterministically strip LaTeX markup to produce caption.txt
#'
#' Mapping (in order):
#' \enumerate{
#'   \item Line comments (\code{\%...} to end of line) removed.
#'   \item \code{\\figpanel\{X\}\{TEXT\}} → \code{(X) TEXT }.
#'   \item \code{\\texttt\{X\}} → \code{X}.
#'   \item \code{\\textit\{X\}} → \code{*X*}.
#'   \item \code{\\textbf\{X\}} → \code{**X**}.
#'   \item \code{\\space} → space; \code{~} → space.
#'   \item Escaped specials (\code{\\&}, \code{\\%}, \code{\\_},
#'     \code{\\$}, \code{\\#}) → bare characters.
#'   \item Whitespace squished (collapsed + trimmed).
#' }
#'
#' Pure-function determinism: identical input → identical output.
#'
#' @param text Character: raw caption.tex content.
#'
#' @return Character scalar of plain-text caption.
#'
#' @importFrom stringr str_replace_all str_squish
#' @keywords internal
#' @noRd
strip_latex <- function(text) {

    # 1. Comments — `%` to end of line (NOT preceded by backslash).
    text <- gsub("(?<!\\\\)%[^\n]*", "", text, perl = TRUE)

    # 2. \figpanel{a}{TEXT} — preserve the panel letter + text.
    text <- gsub(
        "\\\\figpanel\\{([a-z])\\}\\{((?:[^{}]|\\{[^{}]*\\})*)\\}",
        "(\\1) \\2 ",
        text, perl = TRUE
    )

    # 3-5. Inline markup.
    text <- gsub(
        "\\\\texttt\\{((?:[^{}]|\\{[^{}]*\\})*)\\}", "\\1",
        text, perl = TRUE
    )
    text <- gsub(
        "\\\\textit\\{((?:[^{}]|\\{[^{}]*\\})*)\\}", "*\\1*",
        text, perl = TRUE
    )
    text <- gsub(
        "\\\\textbf\\{((?:[^{}]|\\{[^{}]*\\})*)\\}", "**\\1**",
        text, perl = TRUE
    )

    # 6. \space and ~ → space.
    text <- gsub("\\\\space", " ", text, perl = TRUE)
    text <- gsub("~", " ", text, fixed = TRUE)

    # 7. Escaped specials.
    text <- gsub("\\\\&", "&", text, perl = TRUE)
    text <- gsub("\\\\%", "%", text, perl = TRUE)
    text <- gsub("\\\\_", "_", text, perl = TRUE)
    text <- gsub("\\\\\\$", "$", text, perl = TRUE)
    text <- gsub("\\\\#", "#", text, perl = TRUE)

    # 8. Whitespace squish.
    stringr::str_squish(text)
}
