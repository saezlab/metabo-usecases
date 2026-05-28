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


test_that("rebuild.R discover_targets returns sorted, stable plan", {

    # The rebuild driver lives at the repo root, not in the package
    # namespace, so source it into a sandbox env to get
    # discover_targets() available.
    repo_root  <- testthat::test_path("..", "..")
    rebuild_R  <- file.path(repo_root, "rebuild.R")
    skip_if_not(file.exists(rebuild_R), "rebuild.R not found")

    sandbox <- new.env(parent = globalenv())
    sandbox$main <- function() invisible(NULL)  # block auto-exec
    source(rebuild_R, local = sandbox)

    withr::with_tempdir({
        dir.create("figures/fig01-overview",         recursive = TRUE)
        dir.create("figures/fig03-metalinks-versions", recursive = TRUE)
        dir.create("tables/tab01-id-resolving",      recursive = TRUE)
        writeLines("1+1", "figures/fig01-overview/build.R")
        writeLines("1+1", "figures/fig03-metalinks-versions/build.R")
        writeLines("1+1", "tables/tab01-id-resolving/build.R")

        builds_1 <- sandbox$discover_targets(character(0))
        builds_2 <- sandbox$discover_targets(character(0))

        expect_length(builds_1, 3L)
        expect_equal(builds_1, sort(builds_1))
        expect_equal(builds_1, builds_2)
    })
})
