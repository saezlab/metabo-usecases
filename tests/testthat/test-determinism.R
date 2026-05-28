test_that("snapshot_id is byte-stable across re-serialization", {

    m <- list(
        build         = "main",
        built_at      = "2026-05-28T10:00:00+0000",
        packages      = list(
            `omnipath-build`     = "5a3e6a2",
            `omnipath-utils`     = "7b21cc9",
            `omnipath-resources` = "9c2e571"
        ),
        resources     = list(
            list(name = "chebi",  version = "230",     record_count = 198432L,
                 expected_count = 198432L),
            list(name = "signor", version = "2024-03", record_count = 9001L,
                 expected_count = NULL)
        ),
        partial_build = FALSE
    )

    expect_equal(snapshot_id(m), snapshot_id(m))

    # Cosmetic field changes do not change the identifier.
    m2 <- m
    m2$built_at <- "2027-01-01T00:00:00+0000"
    expect_equal(snapshot_id(m), snapshot_id(m2))
})


test_that("rebuild.R --dry-run output is sorted (stable across runs)", {

    skip_if_not_installed("processx")

    withr::with_tempdir({
        dir.create("logs")
        dir.create("figures/fig01-overview", recursive = TRUE)
        dir.create("figures/fig03-metalinks-versions", recursive = TRUE)
        dir.create("tables/tab01-id-resolving", recursive = TRUE)
        writeLines("1+1", "figures/fig01-overview/build.R")
        writeLines("1+1", "figures/fig03-metalinks-versions/build.R")
        writeLines("1+1", "tables/tab01-id-resolving/build.R")

        repo_root <- testthat::test_path("..", "..")
        rebuild_R <- normalizePath(file.path(repo_root, "rebuild.R"))

        run1 <- system2(
            "Rscript",
            c(rebuild_R, "--dry-run"),
            stdout = TRUE, stderr = FALSE
        )
        run2 <- system2(
            "Rscript",
            c(rebuild_R, "--dry-run"),
            stdout = TRUE, stderr = FALSE
        )

        plan_lines <- grep("^\\[plan\\]", run1, value = TRUE)
        expect_true(length(plan_lines) >= 3L)
        expect_equal(plan_lines, sort(plan_lines))
        expect_equal(
            grep("^\\[plan\\]", run1, value = TRUE),
            grep("^\\[plan\\]", run2, value = TRUE)
        )
    })
})
