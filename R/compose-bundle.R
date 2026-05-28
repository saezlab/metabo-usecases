#' Concatenate every pipeline PDF into the manuscript bundle
#'
#' Implements FR-029: one combined PDF containing every PDF artifact
#' (main figures, supplementary, additional plots, tables) in
#' manuscript order. Uses \code{pdftk} (research.md R-2) shelled from
#' \code{lib/log.sh::log_line} so the invocation lands in the
#' unified pipeline log under \code{[bash][INFO][pdftk:bundle]}.
#'
#' @param pdf_paths Character vector of input PDF paths in the desired
#'     bundle order.
#' @param output Path to the output bundle PDF; defaults to
#'     \code{out/manuscript-bundle.pdf}.
#'
#' @return Invisibly the bundle path.
#'
#' @examples
#' \dontrun{
#' assemble_bundle(
#'     c("figures/fig01-overview/out/fig01-overview.pdf",
#'       "figures/fig03-metalinks-versions/out/fig03-metalinks-versions.pdf"),
#'     output = "out/manuscript-bundle.pdf"
#' )
#' }
#'
#' @importFrom fs dir_create path_dir
#' @importFrom rlang abort
#' @importFrom logger log_info
#' @export
assemble_bundle <- function(pdf_paths,
                            output = "out/manuscript-bundle.pdf") {

    if (length(pdf_paths) == 0L) {
        rlang::abort("assemble_bundle() requires at least one input PDF")
    }
    missing <- pdf_paths[!file.exists(pdf_paths)]
    if (length(missing) > 0L) {
        rlang::abort(sprintf(
            "Missing input PDFs:\n  %s",
            paste(missing, collapse = "\n  ")
        ))
    }

    if (Sys.which("pdftk") == "") {
        rlang::abort(
            "pdftk not found on PATH -- install pdftk per quickstart.md section 1"
        )
    }

    fs::dir_create(fs::path_dir(output))

    sh <- sprintf(
        paste0(
            "set -e; source lib/log.sh; ",
            "log_line INFO 'pdftk:bundle' 'Concatenating %d PDFs into %s'; ",
            "pdftk %s cat output %s; ",
            "log_line INFO 'pdftk:bundle' 'Bundle ready: %s'"
        ),
        length(pdf_paths),
        shQuote(output),
        paste(shQuote(pdf_paths), collapse = " "),
        shQuote(output),
        shQuote(output)
    )

    rc <- system2("bash", c("-c", sh))
    if (rc != 0L) {
        rlang::abort(sprintf("pdftk failed (rc=%d)", rc))
    }

    logger::log_info("Wrote bundle {output}")
    invisible(output)
}
