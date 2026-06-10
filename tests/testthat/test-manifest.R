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
        metabo.figures:::infer_build_kind(list(
            `omnipath-utils`     = "u",
            `omnipath-resources` = "r"
        )),
        "utils"
    )
    expect_equal(
        metabo.figures:::infer_build_kind(list(
            `omnipath-build`     = "b",
            `omnipath-utils`     = "u",
            `omnipath-resources` = "r"
        )),
        "main"
    )
    expect_equal(
        metabo.figures:::infer_build_kind(list(
            `omnipath-metabo`    = "m",
            `omnipath-build`     = "b",
            `omnipath-utils`     = "u",
            `omnipath-resources` = "r"
        )),
        "metabo"
    )
})

test_that("infer_build_kind accepts the cycle-001 underscore key form", {

    # The live build_manifest emits keys with underscores
    # (omnipath_build, omnipath_resources, ...) while the FR-032
    # canonical form uses hyphens. The classifier MUST recognise both.

    expect_equal(
        metabo.figures:::infer_build_kind(list(
            omnipath_build     = "b",
            omnipath_resources = "r"
        )),
        "main"
    )
    expect_equal(
        metabo.figures:::infer_build_kind(list(
            omnipath_metabo    = "m",
            omnipath_build     = "b",
            omnipath_resources = "r"
        )),
        "metabo"
    )
})

test_that("packages_for_build returns expected sets per FR-032", {

    expect_setequal(
        metabo.figures:::packages_for_build("utils"),
        c("omnipath-utils", "omnipath-resources")
    )
    expect_setequal(
        metabo.figures:::packages_for_build("main"),
        c("omnipath-build", "omnipath-utils", "omnipath-resources")
    )
    expect_setequal(
        metabo.figures:::packages_for_build("metabo"),
        c(
            "omnipath-metabo", "omnipath-build",
            "omnipath-utils", "omnipath-resources"
        )
    )
})

test_that("packages_for_build excludes omnipath-present everywhere", {

    for (b in c("utils", "main", "metabo")) {
        expect_false(
            "omnipath-present" %in% metabo.figures:::packages_for_build(b),
            info = sprintf("build = %s", b)
        )
    }
})

test_that("parse_jsonb handles character / list / empty / null shapes", {

    expect_equal(metabo.figures:::parse_jsonb(NULL), list())
    expect_equal(metabo.figures:::parse_jsonb(""), list())

    parsed <- metabo.figures:::parse_jsonb('{"a": 1, "b": "two"}')
    expect_equal(parsed$a, 1)
    expect_equal(parsed$b, "two")

    # When the column is already an R list (e.g. RPostgres deserialized
    # the jsonb upstream), return it unchanged.
    expect_equal(
        metabo.figures:::parse_jsonb(list(a = 1L)),
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

test_that("flatten_package_commits handles both rich + legacy shapes", {

    # Cycle-001 rich shape: {commit: "<sha>", dirty: <bool>}
    rich <- list(
        omnipath_build = list(
            commit = "f5ac1f4514fd5995f3a4ff9c7716e5586f73a7e2",
            dirty  = FALSE
        ),
        omnipath_resources = list(
            commit = "cad87fe4dac9a0f4217e1d3eefc197356a5204f5",
            dirty  = TRUE
        )
    )
    flat <- metabo.figures:::flatten_package_commits(rich)
    expect_equal(flat$omnipath_build,
                 "f5ac1f4514fd5995f3a4ff9c7716e5586f73a7e2")
    expect_equal(flat$omnipath_resources,
                 "cad87fe4dac9a0f4217e1d3eefc197356a5204f5")
    # No dirty flag survives.
    expect_true(is.character(flat$omnipath_build))
    expect_length(flat$omnipath_build, 1L)

    # Legacy flat shape: bare hash strings — preserved as-is.
    legacy <- list(`omnipath-build` = "abc1234", `omnipath-utils` = "def5678")
    expect_equal(metabo.figures:::flatten_package_commits(legacy), legacy)

    # Anything else falls back to "unknown" so schema validation passes.
    expect_equal(
        metabo.figures:::flatten_package_commits(list(weird = list(no_commit = "bug"))),
        list(weird = "unknown")
    )
})

test_that("build_manifest_for parses a build_manifest row natively (rich)", {

    # Synthetic DBI fixture matching the live cycle-001 build_manifest
    # row shape (rich package_commits with dirty flag).
    fake_con <- structure(list(), class = "MockConnection")
    mock_row <- data.frame(
        build_id        = "a3f9c2e74b81",
        built_at        = as.POSIXct(
            "2026-05-28T09:00:00", tz = "UTC", format = "%Y-%m-%dT%H:%M:%S"
        ),
        package_commits = '{"omnipath_build": {"commit": "aaa1111", "dirty": false}, "omnipath_resources": {"commit": "bbb2222", "dirty": true}}',
        resources       = '[{"name": "chebi", "version": "230", "record_count": 198432, "expected_count": 198432}]',
        partial_build   = FALSE,
        stringsAsFactors = FALSE
    )

    with_mocked_bindings(
        dbGetQuery = function(con, sql, ...) mock_row,
        .package = "DBI",
        {
            manifest <- build_manifest_for(fake_con)
        }
    )

    expect_equal(manifest$build_id, "a3f9c2e74b81")
    expect_equal(manifest$build, "main")
    # packages is flat (schema-compliant)
    expect_equal(manifest$packages$omnipath_build, "aaa1111")
    expect_equal(manifest$packages$omnipath_resources, "bbb2222")
    # packages_raw preserves the dirty flag
    expect_equal(manifest$packages_raw$omnipath_resources$dirty, TRUE)
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
