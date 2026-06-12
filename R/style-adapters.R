#' Apply the metabo.figures table-style to a gt object
#'
#' Style adapter for \code{gt}-typeset tables (FR-016 / FR-017a). Sets
#' the same sans-serif family and minimum-point-size envelope used by
#' the ggplot theme, so figure + table typography is uniform across
#' the manuscript.
#'
#' @param x A \code{gt_tbl} object.
#'
#' @return The styled \code{gt_tbl}.
#'
#' @examples
#' \dontrun{
#' library(gt)
#' gt(head(mtcars)) %>% gt_metabo_style()
#' }
#'
#' @importFrom rlang check_installed
#' @export
gt_metabo_style <- function(x) {

    rlang::check_installed("gt")

    sizes <- font_sizes()

    x %>%
        gt::tab_options(
            table.font.names         = c(
                "Helvetica Neue LT Std", "Arial", "Helvetica", "sans-serif"
            ),
            table.font.size           = gt::px(sizes$body * 1.3),
            heading.title.font.size   = gt::px(sizes$panel_title * 1.3),
            heading.subtitle.font.size = gt::px(sizes$body * 1.3),
            column_labels.font.size   = gt::px(sizes$axis_label * 1.3),
            column_labels.font.weight = "bold",
            row_group.font.weight     = "bold",
            table_body.hlines.color   = "white",
            table_body.border.bottom.color = "black",
            table_body.border.top.color    = "black",
            column_labels.border.bottom.color = "black",
            data_row.padding           = gt::px(2),
            # gt-native alternating row colours: applies the stripe
            # to every body cell (stub column included) without the
            # wrapper-level \rowcolors{}{} gaps at column separators.
            row.striping.include_table_body = TRUE,
            row.striping.background_color   = "#F2F2F2",
            # Drop the stub vertical rule that previously appeared
            # as a stark | between the stub and the first data
            # column — the row striping reads as the visual divider.
            stub.border.style = "none",
            stub.border.width = gt::px(0)
        ) %>%
        gt::opt_row_striping()
}
