test_that("R setup_pipeline_log writes lines in the contract format", {

    tmp <- tempfile(fileext = ".log")
    withr::with_envvar(c(METABO_FIGURES_LOG = tmp), {
        setup_pipeline_log("build:fixture")
        logger::log_info("hello from R")
    })

    lines <- readLines(tmp)
    expect_length(lines, 1L)

    pattern <- paste0(
        "^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}[+-]\\d{4} ",
        "\\[R\\]\\[INFO\\]\\[build:fixture\\] hello from R$"
    )
    expect_match(lines[[1L]], pattern)
})

test_that("setup_pipeline_log falls back to logs/orphan-PID.log without env", {

    withr::with_envvar(c(METABO_FIGURES_LOG = ""), {
        withr::with_tempdir({
            expect_warning(setup_pipeline_log("test:fallback"))
            expect_true(grepl("logs/orphan-", Sys.getenv("METABO_FIGURES_LOG")))
            expect_true(dir.exists("logs"))
        })
    })
})

test_that("bash log_line in lib/log.sh writes a conforming line", {

    skip_on_os("windows")

    tmp <- tempfile(fileext = ".log")
    repo_root <- testthat::test_path("..", "..")
    log_sh <- normalizePath(file.path(repo_root, "lib", "log.sh"))

    system2(
        "bash",
        c("-c", paste0(
            "set -e; ",
            "export METABO_FIGURES_LOG=", shQuote(tmp), "; ",
            "source ", shQuote(log_sh), "; ",
            "log_line INFO 'build:fixture' 'hello from bash'"
        ))
    )

    lines <- readLines(tmp)
    expect_length(lines, 1L)

    pattern <- paste0(
        "^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}[+-]\\d{2}:\\d{2} ",
        "\\[bash\\]\\[INFO\\]\\[build:fixture\\] hello from bash$"
    )
    expect_match(lines[[1L]], pattern)
})

test_that("log lines stay under the 4 KiB atomicity envelope", {

    tmp <- tempfile(fileext = ".log")
    withr::with_envvar(c(METABO_FIGURES_LOG = tmp), {
        setup_pipeline_log("test:length")
        long_msg <- strrep("x", 5000L)
        logger::log_info("{long_msg}")
    })

    line <- readLines(tmp)
    expect_lte(nchar(line[[1L]], type = "bytes"), 4096L)
})
