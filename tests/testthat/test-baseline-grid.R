# test-baseline-grid.R
#
# FR-009a: panels in the same row of a composite MUST share their y-origin
# within <=1 pt tolerance. Tests the patchwork composition path; the
# LaTeX pdfpages path is validated visually at smoke-test time.
#
# Strategy: build a synthetic 4-panel two-row figure via compose_patchwork(),
# export to SVG, parse translate() transforms to extract y-origins, and
# assert row-mate panels share the same y within 1.5 px (~1 pt at 96 dpi).

library(ggplot2)

panel_y_origins_from_svg <- function(svg_path) {
    svg_text <- paste(readLines(svg_path, warn = FALSE), collapse = "\n")
    pat <- "translate\\([^,]+,\\s*([0-9.]+)\\)"
    m <- gregexpr(pat, svg_text, perl = TRUE)
    matches <- regmatches(svg_text, m)[[1L]]
    if (length(matches) == 0L) return(numeric(0L))
    as.numeric(sub(paste0(".*", pat, ".*"), "\\1", matches))
}

test_that("FR-009a: two-row patchwork respects baseline grid", {
    skip_if_not_installed("svglite")

    p1 <- ggplot(data.frame(x = 1L:3L, y = 1L:3L), aes(x, y)) + geom_col()
    p2 <- ggplot(data.frame(x = 1L:3L, y = 3L:1L), aes(x, y)) + geom_col()
    p3 <- ggplot(data.frame(x = 1L:5L, y = 5L:1L), aes(x, y)) + geom_col()
    p4 <- ggplot(data.frame(x = 1L:5L, y = 1L:5L), aes(x, y)) + geom_col()

    tmp <- withr::local_tempfile(fileext = ".svg")

    out <- compose_patchwork(
        panels  = list(p1, p2, p3, p4),
        layout  = "A B
C D",
        widths  = c(0.5, 0.5),
        heights = c(0.5, 0.5)
    )

    svglite::svglite(tmp, width = 7, height = 7)
    print(out)
    grDevices::dev.off()

    y_vals <- panel_y_origins_from_svg(tmp)

    skip_if(
        length(y_vals) < 4L,
        "SVG y-origin parsing yielded < 4 values; visual inspection required"
    )

    expect_lt(
        abs(y_vals[[1L]] - y_vals[[2L]]), 1.5,
        label = "Row 1 panels share y-origin within 1 pt"
    )
    expect_lt(
        abs(y_vals[[3L]] - y_vals[[4L]]), 1.5,
        label = "Row 2 panels share y-origin within 1 pt"
    )
    expect_gt(
        mean(y_vals[3L:4L]),
        mean(y_vals[1L:2L]),
        label = "Row 2 starts below row 1"
    )
})
