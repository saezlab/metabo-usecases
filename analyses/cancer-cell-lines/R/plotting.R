# Plot theme, colour-blind-safe palettes and figure export for the notebooks.

#' ggplot2 theme used by all notebooks of the use case
theme_uc <- function(base_size = 11) {
    ggplot2::theme_bw(base_size = base_size) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            strip.background = ggplot2::element_rect(fill = "grey95"),
            legend.position = "right"
        )
}

#' Colour-blind-safe qualitative palette
#'
#' Okabe-Ito for up to 8 categories, otherwise Paul Tol's discrete rainbow
#' (up to 14). Grey is reserved for missing/unknown categories.
#'
#' @param levels Category names; the result is named by them.
#' @return Named character vector of colours.
uc_palette <- function(levels) {
    okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
    tol_rainbow <- c("#882E72", "#1965B0", "#7BAFDE", "#4EB265", "#CAE0AB", "#F7F056", "#F6C141",
                     "#F1932D", "#E8601C", "#DC050C", "#72190E", "#5289C7", "#90C987", "#D1BBD7")
    unknown <- levels %in% c("Unknown", "unknown", NA, "NA", "not reported")
    known <- levels[!unknown]
    pal <- if (length(known) <= 8) okabe_ito else tol_rainbow
    if (length(known) > length(pal)) stop("More categories than palette colours: ", length(known))
    out <- c(setNames(pal[seq_along(known)], known), setNames(rep("grey60", sum(unknown)), levels[unknown]))
    out[levels]
}

#' Scrollable preview of a table or matrix in an HTML report
#'
#' Shows the first rows and columns only, so that reports stay small.
#'
#' @param x Data frame or matrix.
#' @param rows,cols Maximum number of rows and columns shown.
#' @param caption Caption above the table.
#' @return A DT htmlwidget.
preview_table <- function(x, rows = 100, cols = 40, caption = NULL) {
    full <- dim(x)
    x <- as.data.frame(x[seq_len(min(rows, nrow(x))), seq_len(min(cols, ncol(x))), drop = FALSE],
                       check.names = FALSE)
    caption <- paste0(if (!is.null(caption)) paste0(caption, " – ") else "",
                      sprintf("showing %d of %d rows, %d of %d columns", nrow(x), full[1], ncol(x), full[2]))
    DT::datatable(
        x,
        caption = caption,
        rownames = !is.null(rownames(x)) && !identical(rownames(x), as.character(seq_len(nrow(x)))),
        options = list(scrollX = TRUE, scrollY = "300px", paging = FALSE, dom = "t")
    )
}

#' Save a ggplot as vector PDF and PNG
#'
#' @param plot ggplot object.
#' @param name File name without extension.
#' @param dir Output directory.
#' @param width,height Size in inches.
save_figure <- function(plot, name, dir, width = 7, height = 5) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(file.path(dir, paste0(name, ".pdf")), plot, width = width, height = height)
    ggplot2::ggsave(file.path(dir, paste0(name, ".png")), plot, width = width, height = height, dpi = 200)
    invisible(plot)
}
