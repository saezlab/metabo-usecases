#' Figure 6 Panels A / B -- volcano plot for one mutation contrast
#'
#' Renders a volcano (log fold-change vs. -log10 P.Value) for the
#' Shorthouse 2022 limma output. Significance categories are coloured
#' using the lead palette: above threshold = first non-grey entry
#' (teal); below threshold = \code{"#BEBEBE"} (FR-020). X and Y axis
#' limits can be supplied via \code{xlim} / \code{ylim} to keep
#' Panels A and B aligned (FR-019); use
#' \code{\link{volcano_shared_limits}} to compute them across both
#' contrasts before rendering.
#'
#' @param diff_tibble Output of \code{\link{case_study_differential}}
#'     for the contrast being plotted; must carry \code{logFC} and
#'     \code{P.Value} columns (a \code{neg_log10_p} column is added
#'     if absent).
#' @param contrast_label Character: panel title shown above the plot
#'     (e.g. \code{"KRAS"}, \code{"EGFR"}).
#' @param pvalue_threshold Numeric: P-value cutoff for the
#'     significance colour split. Default \code{0.05}.
#' @param logfc_threshold Numeric: absolute log fold-change cutoff
#'     for the significance colour split. Default \code{0.5}.
#' @param xlim Numeric length-2 vector or \code{NULL}: shared x-axis
#'     limits. \code{NULL} → ggplot picks per-panel.
#' @param ylim Numeric length-2 vector or \code{NULL}: shared y-axis
#'     limits (on the -log10 P-value scale).
#' @param width_mm Numeric: target physical panel width in mm; passed
#'     to \code{\link{theme_bw_metabo}}.
#' @param font_scale Numeric: passed to \code{\link{theme_bw_metabo}}.
#' @param point_size Numeric: point size for \code{geom_point}.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' kras <- case_study_differential("KRAS")
#' egfr <- case_study_differential("EGFR")
#' lims <- volcano_shared_limits(list(kras, egfr))
#' volcano_panel(kras, "KRAS", xlim = lims$xlim, ylim = lims$ylim)
#' volcano_panel(egfr, "EGFR", xlim = lims$xlim, ylim = lims$ylim)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline
#' @importFrom ggplot2 scale_colour_manual scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous labs theme
#' @importFrom rlang .data abort
#' @importFrom dplyr mutate
#' @export
volcano_panel <- function(
    diff_tibble,
    contrast_label,
    pvalue_threshold = 0.05,
    logfc_threshold = 0.5,
    xlim = NULL,
    ylim = NULL,
    width_mm = 89L,
    font_scale = 1,
    point_size = 0.6
) {

    # NSE vs. R CMD check workaround
    P.Value <- logFC <- NULL

    required <- c("logFC", "P.Value")
    missing <- setdiff(required, names(diff_tibble))
    if (length(missing) > 0L) {
        rlang::abort(sprintf(
            "diff_tibble is missing required columns: %s",
            paste(missing, collapse = ", ")
        ))
    }

    if (!"neg_log10_p" %in% names(diff_tibble)) {
        diff_tibble <- dplyr::mutate(
            diff_tibble,
            neg_log10_p = -log10(P.Value)
        )
    }

    diff_tibble <- dplyr::mutate(
        diff_tibble,
        sig = ifelse(
            P.Value < pvalue_threshold &
                abs(logFC) > logfc_threshold,
            "Significant",
            "Not significant"
        )
    )

    sig_colours <- c(
        Significant       = unname(palette_lead()[["teal"]]),
        `Not significant` = unname(palette_lead()[["unknown"]])
    )

    plt <- ggplot2::ggplot(
        diff_tibble,
        ggplot2::aes(
            x      = .data$logFC,
            y      = .data$neg_log10_p,
            colour = .data$sig
        )
    ) +
        ggplot2::geom_vline(
            xintercept = c(-logfc_threshold, logfc_threshold),
            linetype   = "dashed",
            colour     = "grey70",
            linewidth  = 0.25
        ) +
        ggplot2::geom_hline(
            yintercept = -log10(pvalue_threshold),
            linetype   = "dashed",
            colour     = "grey70",
            linewidth  = 0.25
        ) +
        ggplot2::geom_point(
            size  = point_size,
            alpha = 0.7
        ) +
        ggplot2::scale_colour_manual(
            values = sig_colours,
            breaks = c("Significant", "Not significant"),
            name   = NULL
        ) +
        ggplot2::labs(
            title = contrast_label,
            x     = expression(log[2] ~ "fold-change"),
            y     = expression(-log[10] ~ italic(P))
        ) +
        theme_bw_metabo(
            width_mm   = width_mm,
            font_scale = font_scale
        ) +
        ggplot2::theme(legend.position = "top")

    if (!is.null(xlim)) {
        plt <- plt + ggplot2::scale_x_continuous(limits = xlim)
    }
    if (!is.null(ylim)) {
        plt <- plt + ggplot2::scale_y_continuous(limits = ylim)
    }

    plt
}


#' Compute shared volcano-plot axis limits across multiple contrasts
#'
#' Returns the union of x and y ranges across every input tibble so
#' that two or more volcano panels share the same coordinate system
#' (FR-019). The x range is symmetric around zero by default to keep
#' the up- and down-regulated half-planes visually balanced.
#'
#' @param diff_tibbles List of tibbles, each carrying \code{logFC}
#'     and \code{P.Value} columns.
#' @param symmetric_x Logical: when \code{TRUE} (default) the x range
#'     is widened to \code{c(-m, m)} where \code{m} is the max
#'     absolute \code{logFC} observed across all inputs.
#' @param pad_frac Numeric: fractional padding added to each end of
#'     the range. Default \code{0.05}.
#'
#' @return Named list \code{xlim} (length-2) and \code{ylim}
#'     (length-2). Y-limits are on the -log10 P-value scale.
#'
#' @examples
#' \dontrun{
#' lims <- volcano_shared_limits(
#'     list(
#'         case_study_differential("KRAS"),
#'         case_study_differential("EGFR")
#'     )
#' )
#' }
#'
#' @importFrom rlang abort
#' @export
volcano_shared_limits <- function(
    diff_tibbles,
    symmetric_x = TRUE,
    pad_frac = 0.05
) {

    if (length(diff_tibbles) == 0L) {
        rlang::abort(
            "diff_tibbles must contain at least one tibble"
        )
    }

    logfc <- unlist(lapply(diff_tibbles, function(d) d$logFC))
    pvals <- unlist(lapply(diff_tibbles, function(d) d$P.Value))
    neg_log10_p <- -log10(pvals)

    finite_logfc <- logfc[is.finite(logfc)]
    finite_y <- neg_log10_p[is.finite(neg_log10_p)]

    if (length(finite_logfc) == 0L || length(finite_y) == 0L) {
        rlang::abort(
            "No finite logFC / P.Value entries to compute limits"
        )
    }

    x_max_abs <- max(abs(finite_logfc))
    x_range <- if (isTRUE(symmetric_x)) {
        c(-x_max_abs, x_max_abs)
    } else {
        range(finite_logfc)
    }

    y_min <- 0
    y_max <- max(finite_y)
    y_range <- c(y_min, y_max)

    x_pad <- diff(x_range) * pad_frac
    y_pad <- diff(y_range) * pad_frac

    list(
        xlim = c(x_range[1L] - x_pad, x_range[2L] + x_pad),
        ylim = c(y_range[1L], y_range[2L] + y_pad)
    )
}
