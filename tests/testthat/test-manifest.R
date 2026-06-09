test_that("snapshot_id returns the manifest's build_id", {

    manifest <- list(
        build_id      = "a3f9c2e74b81",
        built_at      = "2026-05-28T09:00:00+0000",
        build         = "main",
        packages      = list(
            `omnipath-build`     = "aaa1111",
            `omnipath-utils`     = "bbb2222",
            `omnipath-resources` = "ccc3333"
        ),
        resources     = list(),
        partial_build = FALSE
    )

    expect_equal(snapshot_id(manifest), "a3f9c2e74b81")
    expect_match(snapshot_id(manifest), "^[0-9a-f]{12}$")
})

test_that("snapshot_id refuses a pre-cycle-001 manifest", {

    manifest <- list(
        build    = "main",
        packages = list(`omnipath-build` = "aaa1111"),
        resources = list(),
        partial_build = FALSE
    )

    expect_error(snapshot_id(manifest), "build_id")
})

test_that("infer_build_kind maps package sets to FR-032 build labels", {

    expect_equal(
        infer_build_kind(list(
            `omnipath-utils`     = "u",
            `omnipath-resources` = "r"
        )),
        "utils"
    )
    expect_equal(
        infer_build_kind(list(
            `omnipath-build`     = "b",
            `omnipath-utils`     = "u",
            `omnipath-resources` = "r"
        )),
        "main"
    )
    expect_equal(
        infer_build_kind(list(
            `omnipath-metabo`    = "m",
            `omnipath-build`     = "b",
            `omnipath-utils`     = "u",
            `omnipath-resources` = "r"
        )),
        "metabo"
    )
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

test_that("packages_for_build excludes omnipath-present everywhere", {

    for (b in c("utils", "main", "metabo")) {
        expect_false(
            "omnipath-present" %in% packages_for_build(b),
            info = sprintf("build = %s", b)
        )
    }
})

test_that("parse_jsonb handles character / list / empty / null shapes", {

    expect_equal(parse_jsonb(NULL), list())
    expect_equal(parse_jsonb(""), list())

    parsed <- parse_jsonb('{"a": 1, "b": "two"}')
    expect_equal(parsed$a, 1)
    expect_equal(parsed$b, "two")

    # When the column is already an R list (e.g. RPostgres deserialized
    # the jsonb upstream), return it unchanged.
    expect_equal(
        parse_jsonb(list(a = 1L)),
        list(a = 1L)
    )
})

test_that("write_manifest emits canonical JSON and SHA256 files", {

    skip_if_not_installed("fs")

    withr::with_tempdir({
        m <- list(
            build_id      = "abcdef012345",
            built_at      = "2026-05-28T09:00:00+0000",
            build         = "utils",
            packages      = list(
                `omnipath-utils`     = "aaa1111",
                `omnipath-resources` = "bbb2222"
            ),
            resources     = list(),
            partial_build = FALSE
        )
        sid <- write_manifest(m, deployment = "dev3", dir = "manifests")

        expect_equal(sid, "abcdef012345")
        expect_true(file.exists(
            file.path("manifests", "dev3.abcdef012345.json")
        ))
        expect_true(file.exists(
            file.path("manifests", "dev3.abcdef012345.SHA256")
        ))
        expect_equal(
            readLines(file.path("manifests", "dev3.abcdef012345.SHA256")),
            sid
        )
    })
})

test_that("build_manifest_for parses a build_manifest row natively", {

    # Synthetic DBI fixture so we don't need a live Postgres.
    fake_con <- structure(list(), class = "MockConnection")
    mock_row <- data.frame(
        build_id        = "a3f9c2e74b81",
        built_at        = as.POSIXct(
            "2026-05-28T09:00:00", tz = "UTC", format = "%Y-%m-%dT%H:%M:%S"
        ),
        package_commits = '{"omnipath-build": "aaa1111", "omnipath-utils": "bbb2222", "omnipath-resources": "ccc3333"}',
        resources       = '[{"name": "chebi", "version": "230", "record_count": 198432, "expected_count": 198432}]',
        partial_build   = FALSE,
        stringsAsFactors = FALSE
    )

    # Stub DBI::dbGetQuery via local mocking.
    withr::local_envvar(c("R_TESTS" = ""))
    with_mocked_bindings(
        dbGetQuery = function(con, sql, ...) mock_row,
        .package = "DBI",
        {
            manifest <- build_manifest_for(fake_con)
        }
    )

    expect_equal(manifest$build_id, "a3f9c2e74b81")
    expect_equal(manifest$build, "main")
    expect_equal(manifest$packages$`omnipath-build`, "aaa1111")
    expect_equal(length(manifest$resources), 1L)
    expect_equal(manifest$resources[[1L]]$name, "chebi")
    expect_false(manifest$partial_build)
})

test_that("build_manifest_for errors when the table is empty", {

    fake_con <- structure(list(), class = "MockConnection")
    empty <- data.frame(
        build_id        = character(),
        built_at        = as.POSIXct(character()),
        package_commits = character(),
        resources       = character(),
        partial_build   = logical(),
        stringsAsFactors = FALSE
    )

    with_mocked_bindings(
        dbGetQuery = function(con, sql, ...) empty,
        .package = "DBI",
        {
            expect_error(build_manifest_for(fake_con), "build_manifest table is empty")
        }
    )
})
