test_that("snapshot_id reads build_id directly (post-cycle-001)", {

    # Cycle-001 build_manifest stores the snapshot identifier
    # directly as `build_id` (12-hex SHA-256). snapshot_id() is now
    # a thin accessor on that column — no in-pipeline re-derivation,
    # so the identifier is constant regardless of any other manifest
    # field changes.
    m <- list(
        build_id      = "a3f9c2e74b81",
        built_at      = "2026-05-28T10:00:00+0000",
        build         = "main",
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

    expect_equal(snapshot_id(m), "a3f9c2e74b81")
    expect_equal(snapshot_id(m), snapshot_id(m))

    # Cosmetic field changes do not change the identifier.
    m2 <- m
    m2$built_at <- "2027-01-01T00:00:00+0000"
    expect_equal(snapshot_id(m), snapshot_id(m2))
})


test_that("snapshot_id refuses a pre-cycle-001 manifest", {

    m <- list(
        build         = "main",
        packages      = list(`omnipath-build` = "aaa1111"),
        resources     = list(),
        partial_build = FALSE
    )

    expect_error(snapshot_id(m), "build_id")
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
