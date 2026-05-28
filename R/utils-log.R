#' Configure the per-script pipeline logger
#'
#' Wraps the logger package so every call site writes to the unified
#' pipeline log file at \code{Sys.getenv("METABO_FIGURES_LOG")} using the
#' line format documented in
#' \code{specs/001-figures-pipeline/contracts/log-format.md}:
#'
#' \preformatted{<ISO-8601-with-tz> [R][<level>][<component>] <msg>}
#'
#' Called once at the top of every R script in the pipeline.
#'
#' @param component Character: logical owner of the lines emitted from this
#'     script. Conventionally one of \code{rebuild},
#'     \code{manifest:<build>}, \code{build:<artifact-id>},
#'     \code{compose:<artifact-id>}, \code{gate:<name>},
#'     \code{check:<name>}.
#'
#' @return Invisibly the resolved log file path.
#'
#' @examples
#' setup_pipeline_log("build:fig01-overview")
#' logger::log_info("Rendering Panel B")
#'
#' @importFrom logger log_appender log_layout layout_glue_generator appender_file
#' @importFrom logger log_threshold log_warn
#' @export
setup_pipeline_log <- function(component) {

    log_path <- Sys.getenv("METABO_FIGURES_LOG", unset = "")

    if (!nzchar(log_path)) {
        log_path <- file.path(
            "logs",
            sprintf("orphan-%d.log", Sys.getpid())
        )
        dir.create(dirname(log_path), showWarnings = FALSE, recursive = TRUE)
        Sys.setenv(METABO_FIGURES_LOG = log_path)
        warning(
            sprintf(
                "METABO_FIGURES_LOG was not set; falling back to %s",
                log_path
            ),
            call. = FALSE
        )
    }

    base_layout <- logger::layout_glue_generator(
        format = paste0(
            "{format(time, '%Y-%m-%dT%H:%M:%S%z')} ",
            "[R][{level}][", component, "] {msg}"
        )
    )

    # Truncate each formatted line to the 4 KiB POSIX-O_APPEND
    # atomicity envelope (contracts/log-format.md). Reserve one byte
    # for the trailing newline the appender writes.
    layout <- structure(
        function(level, msg, namespace, .logcall, .topcall, .topenv) {
            line <- base_layout(
                level     = level,
                msg       = msg,
                namespace = namespace,
                .logcall  = .logcall,
                .topcall  = .topcall,
                .topenv   = .topenv
            )
            # The ellipsis is 3 bytes (U+2026); reserve room for it
            # and for the trailing newline the appender writes.
            ifelse(
                nchar(line, type = "bytes") > 4095L,
                paste0(substr(line, 1L, 4093L), "..."),
                line
            )
        },
        generator = deparse(sys.call())
    )

    logger::log_appender(logger::appender_file(log_path))
    logger::log_layout(layout)
    logger::log_threshold(logger::INFO)

    invisible(log_path)
}
