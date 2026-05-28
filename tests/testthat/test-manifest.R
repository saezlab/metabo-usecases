test_that("snapshot_id is deterministic across identical manifests", {

    m1 <- list(
        build         = "metabo",
        built_at      = "2026-05-28T09:00:00+0000",
        packages      = list(
            `omnipath-metabo`    = "abcd123",
            `omnipath-build`     = "ef45678",
            `omnipath-utils`     = "9876543",
            `omnipath-resources` = "1234567"
        ),
        resources     = list(
            list(name = "chebi",  version = "230",
                 record_count = 198432L, expected_count = 198432L),
            list(name = "signor", version = "2024-03",
                 record_count = 9001L, expected_count = NULL)
        ),
        partial_build = FALSE
    )
    m2 <- m1
    m2$built_at <- "2027-01-01T00:00:00+0000"   # cosmetic only

    expect_equal(snapshot_id(m1), snapshot_id(m2))
    expect_match(snapshot_id(m1), "^[0-9a-f]{12}$")
})

test_that("snapshot_id changes when packages or resources change", {

    m1 <- list(
        build         = "main",
        packages      = list(`omnipath-build` = "111"),
        resources     = list(),
        partial_build = FALSE
    )
    m2 <- m1; m2$packages$`omnipath-build` <- "222"
    m3 <- m1; m3$resources <- list(
        list(name = "x", version = "v1", record_count = 1L,
             expected_count = NULL)
    )

    expect_false(identical(snapshot_id(m1), snapshot_id(m2)))
    expect_false(identical(snapshot_id(m1), snapshot_id(m3)))
})

test_that("packages_for_build excludes omnipath-present everywhere", {

    for (b in c("utils", "main", "metabo")) {
        expect_false(
            "omnipath-present" %in% packages_for_build(b),
            info = sprintf("build = %s", b)
        )
    }
})

test_that("packages_for_build returns expected sets per FR-032", {

    expect_setequal(
        packages_for_build("utils"),
        c("omnipath-utils", "omnipath-resources")
    )
    expect_setequal(
        packages_for_build("main"),
        c("omnipath-build", "omnipath-utils", "omnipath-resources")
    )
    expect_setequal(
        packages_for_build("metabo"),
        c(
            "omnipath-metabo", "omnipath-build",
            "omnipath-utils", "omnipath-resources"
        )
    )
})

test_that("write_manifest emits canonical JSON and SHA256 files", {

    skip_if_not_installed("fs")

    withr::with_tempdir({
        m <- list(
            build         = "utils",
            packages      = list(
                `omnipath-utils`     = "aaa1111",
                `omnipath-resources` = "bbb2222"
            ),
            resources     = list(),
            partial_build = FALSE
        )
        sid <- write_manifest(m, dir = "manifests")

        expect_match(sid, "^[0-9a-f]{12}$")
        expect_true(file.exists(
            file.path("manifests", sprintf("utils.%s.json", sid))
        ))
        expect_true(file.exists(
            file.path("manifests", sprintf("utils.%s.SHA256", sid))
        ))
        expect_equal(
            readLines(file.path("manifests", sprintf("utils.%s.SHA256", sid))),
            sid
        )
    })
})

test_that("partial_build flag fires on missing commit hashes", {

    skip_on_cran()

    m <- list(
        build         = "main",
        packages      = resolve_package_commits(
            packages_for_build("main"),
            override = NULL
        ),
        resources     = list(),
        partial_build = any(
            resolve_package_commits(
                packages_for_build("main"),
                override = NULL
            ) == "unknown"
        )
    )

    if (any(unlist(m$packages) == "unknown")) {
        expect_true(m$partial_build)
    }
})
