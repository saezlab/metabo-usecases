test_that("parse_panel_letters reads the \\figpanel macro form", {

    src <- "\\textbf{Figure 1.}\\space
            \\figpanel{a}{Aaa}
            \\figpanel{b}{Bbb}
            \\figpanel{c}{Ccc}"
    expect_equal(metabo.figures:::parse_panel_letters(src), c("a", "b", "c"))
})

test_that("parse_panel_letters falls back to prose `(a) ` form", {

    src <- "Figure 1: title.\n(a) First. (b) Second. (c) Third."
    expect_equal(metabo.figures:::parse_panel_letters(src), c("a", "b", "c"))
})

test_that("strip_latex handles markup and produces deterministic output", {

    src <- paste(
        "% A comment line.",
        "\\textbf{Figure 1: Title.}\\space",
        "\\figpanel{a}{A panel with \\texttt{code} and \\textit{italic}.}",
        "\\figpanel{b}{Special \\& \\% chars.}",
        sep = "\n"
    )
    plain <- metabo.figures:::strip_latex(src)

    expect_match(plain, "^\\*\\*Figure 1: Title\\.\\*\\* \\(a\\) A panel")
    expect_match(plain, "code and \\*italic\\*")
    expect_match(plain, "Special & %")

    # Determinism: same input → byte-identical output.
    expect_identical(metabo.figures:::strip_latex(src), plain)
})

test_that("compose_caption fails when (a)/(b)/... count != panel_count", {

    withr::with_tempdir({
        caption <- "\\figpanel{a}{x} \\figpanel{b}{y}"
        writeLines(caption, "caption.tex")
        writeLines("dummy composite", "fig01.pdf")
        # Provide a stub caption.sty so the pre-flight existence
        # check passes; xelatex is not invoked because the
        # FR-041b check fires first.
        writeLines("\\NeedsTeXFormat{LaTeX2e}", "caption.sty")

        expect_error(
            compose_caption(
                figure_id      = "fig01-overview",
                composite_pdf  = "fig01.pdf",
                caption_source = "caption.tex",
                out_dir        = ".",
                panel_count    = 3L,
                caption_sty    = "caption.sty"
            ),
            "FR-041b"
        )
    })
})

test_that("compose_caption errors on missing caption source", {

    withr::with_tempdir({
        writeLines("dummy", "fig01.pdf")
        writeLines("\\NeedsTeXFormat{LaTeX2e}", "caption.sty")

        expect_error(
            compose_caption(
                figure_id      = "fig01-overview",
                composite_pdf  = "fig01.pdf",
                caption_source = "missing.tex",
                out_dir        = ".",
                panel_count    = 1L,
                caption_sty    = "caption.sty"
            ),
            "Caption source missing"
        )
    })
})
