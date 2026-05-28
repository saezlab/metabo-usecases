#' Verify a plot meets the FR-017a readability envelope
#'
#' Implements SC-008: at the target physical figure width (89 mm
#' single-column / 180 mm double-column), every axis label, panel
#' title, legend entry, and body text element must measure ≥ 6 pt
#' and every tick label must measure ≥ 5 pt.
#'
#' Pure pt sizes declared in \code{theme_bw_metabo()} can survive
#' \code{ggsave} faithfully when the rendered width matches the
#' design width; this gate verifies that and surfaces violations
#' before they reach reviewers.
#'
#' Phase 2B note: the check implementation walks the ggplot grob tree
#' for text elements and asserts their declared \code{font.size} (in
#' pt). It does NOT yet account for \code{ggsave}'s \code{scale}
#' argument — that hardening lands during US1 (T024 / T048) once we
#' know which call sites actually pass a non-default scale.
#'
#' @param plot A ggplot2 plot object.
#' @param width_mm Numeric: target render width in mm; \code{89} or
#'     \code{180}.
#'
#' @return A list with elements:
#'     \itemize{
#'         \item \code{pass}: Logical scalar.
#'         \item \code{violations}: A tibble of failing elements with
#'             columns \code{role}, \code{size_pt}, \code{min_pt}.
#'         \item \code{summary}: Single-line character describing
#'             min observed sizes (formatted for logging).
#'     }
#'
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(wt, mpg)) +
#'     geom_point() +
#'     theme_bw_metabo(width_mm = 89)
#' readability_check(p, width_mm = 89)$pass
#'
#' @importFrom ggplot2 ggplot_build
#' @importFrom tibble tibble
#' @importFrom dplyr filter
#' @importFrom logger log_info log_warn
#' @export
readability_check <- function(plot, width_mm) {

    # NSE workaround
    size_pt <- min_pt <- NULL

    sizes <- font_sizes()

    elements <- tibble::tibble(
        role    = c("axis_label", "panel_title", "legend", "body", "tick"),
        size_pt = c(
            theme_value_pt(plot, "axis.title"),
            theme_value_pt(plot, "plot.title"),
            theme_value_pt(plot, "legend.text"),
            theme_value_pt(plot, "text"),
            theme_value_pt(plot, "axis.text")
        ),
        min_pt  = c(
            sizes$axis_label,
            sizes$panel_title,
            sizes$legend,
            sizes$body,
            sizes$tick
        )
    )

    violations <- elements %>%
        dplyr::filter(!is.na(size_pt) & size_pt < min_pt)

    pass <- nrow(violations) == 0L

    summary_line <- sprintf(
        "min %.1fpt (label), min %.1fpt (tick)",
        min(elements$size_pt[elements$role != "tick"], na.rm = TRUE),
        elements$size_pt[elements$role == "tick"]
    )

    if (pass) {
        logger::log_info(
            "PASS readability @ {width_mm}mm: {summary_line}"
        )
    } else {
        logger::log_warn(
            "FAIL readability @ {width_mm}mm: {nrow(violations)} ",
            "elements below the FR-017a floor"
        )
    }

    list(pass = pass, violations = violations, summary = summary_line)
}


#' Resolve a theme element's effective pt size
#'
#' Walks \code{plot$theme} to find the nearest specified text element
#' for the given role and returns its declared \code{size} in pt.
#' Returns \code{NA_real_} when no value is set, which the caller
#' treats as "not violated" (the role inherits from the base theme).
#'
#' @param plot A ggplot object.
#' @param role Character: the theme element name (e.g.
#'     \code{"axis.title"}).
#' @return Numeric scalar in pt, or \code{NA_real_}.
#'
#' @keywords internal
#' @noRd
theme_value_pt <- function(plot, role) {

    element <- plot$theme[[role]]
    if (is.null(element) || is.null(element$size)) {
        return(NA_real_)
    }

    size <- element$size
    if (inherits(size, "rel")) {
        base <- plot$theme$text$size %||% 11
        return(as.numeric(size) * base)
    }

    as.numeric(size)
}


# Backport of the null-coalescing operator for the rare R version
# without it; harmless if base R already provides one.
`%||%` <- function(a, b) if (is.null(a)) b else a
