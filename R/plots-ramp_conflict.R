#' Panel G -- RaMP InChIKey-conflict counts (FR-007f)
#'
#' Vertical bar plot of conflict counts per reason category. Reasons
#' are deterministic per the cycle-001 RDKit-cartridge
#' classification (\code{stereo}, \code{specificity}, \code{tautomer},
#' \code{similar}, \code{unrelated}); the renderer registers them in
#' the \code{ramp_conflict_reasons} category-colour group on first
#' encounter so the same reason is the same colour across figures
#' (FR-022, SC-004). Uniform bar width matches Panels B-F.
#'
#' Empty-data case emits a single-bar placeholder (Edge Case:
#' zero-row queries get an explicit placeholder, not a blank file).
#'
#' @param data Tibble from \code{\link{ramp_conflict_counts}}.
#' @param width_mm Numeric: target physical width for the panel.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_identity labs
#' @importFrom ggplot2 scale_y_continuous expansion theme element_text
#' @importFrom dplyr mutate
#' @importFrom tibble tibble
#' @importFrom rlang .data
#' @export
plot_ramp_conflict <- function(data, width_mm = 89L) {

    conflict_reason <- n <- fill_hex <- NULL

    if (nrow(data) == 0L) {
        data <- tibble::tibble(
            conflict_reason = "(none)",
            n               = 0L,
            fill_hex        = "#BEBEBE"
        )
    } else {
        data <- dplyr::mutate(
            data,
            fill_hex = category_colour(
                "ramp_conflict_reasons", conflict_reason
            )
        )
    }

    ggplot2::ggplot(
        data,
        ggplot2::aes(
            x    = stats::reorder(.data$conflict_reason, -.data$n),
            y    = .data$n,
            fill = .data$fill_hex
        )
    ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::scale_fill_identity() +
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = c(0, 0.05))
        ) +
        ggplot2::labs(x = NULL, y = "RaMP id conflicts") +
        theme_bw_metabo(width_mm = width_mm) +
        ggplot2::theme(
            legend.position = "none",
            axis.text.x     = ggplot2::element_text(angle = 30, hjust = 1)
        )
}
