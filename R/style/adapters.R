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
            table.font.names         = c("Arial", "Helvetica", "sans-serif"),
            table.font.size           = gt::px(sizes$body * 1.3),
            heading.title.font.size   = gt::px(sizes$panel_title * 1.3),
            heading.subtitle.font.size = gt::px(sizes$body * 1.3),
            column_labels.font.size   = gt::px(sizes$axis_label * 1.3),
            column_labels.font.weight = "bold",
            row_group.font.weight     = "bold",
            table_body.hlines.color   = "grey90",
            table_body.border.bottom.color = "black",
            table_body.border.top.color    = "black",
            column_labels.border.bottom.color = "black"
        )
}
