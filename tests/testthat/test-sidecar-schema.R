test_that("write_sidecar emits deployments: array with build_id", {

    withr::with_tempdir({

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

        deployment <- list(
            name                  = "dev3",
            url                   = "dev3.omnipathdb.org",
            db_host               = "localhost",
            db_port               = 5403L,
            db_name               = "omnipath",
            build_id              = manifest$build_id,
            build_manifest_source = "table"
        )

        artifact_path <- "out/artifact.pdf"
        dir.create("out")
        writeLines("dummy", artifact_path)

        sc <- write_sidecar(
            artifact_id   = "database-content",
            artifact_path = artifact_path,
            deployments   = list(deployment),
            manifests     = list(manifest),
            script_path   = "figures/database-content/build.R",
            queries       = list(
                list(
                    sql         = "SELECT 1",
                    row_count   = 1L,
                    result_hash = "deadbeef0123"
                )
            )
        )

        side_path <- paste0(artifact_path, ".provenance.json")
        expect_true(file.exists(side_path))

        loaded <- jsonlite::fromJSON(side_path, simplifyVector = FALSE)
        expect_equal(loaded$artifact_id, "database-content")
        expect_equal(length(loaded$deployments), 1L)
        expect_equal(loaded$deployments[[1L]]$name, "dev3")
        expect_equal(loaded$deployments[[1L]]$db_port, 5403L)
        expect_equal(loaded$deployments[[1L]]$build_id, "a3f9c2e74b81")
        expect_equal(loaded$deployments[[1L]]$build_manifest_source, "table")
        expect_equal(loaded$package_commits$`omnipath-build`, "aaa1111")
        expect_equal(length(loaded$queries), 1L)
        expect_equal(loaded$queries[[1L]]$sql, "SELECT 1")
        expect_match(loaded$deployments[[1L]]$build_id, "^[0-9a-f]{12}$")
    })
})

test_that("write_sidecar lists every deployment a panel touched", {

    withr::with_tempdir({

        m3 <- list(
            build_id  = "111111111111",
            build     = "main",
            packages  = list(`omnipath-build` = "aaa", `omnipath-utils` = "bbb",
                             `omnipath-resources` = "ccc"),
            resources = list(), partial_build = FALSE
        )
        m4 <- list(
            build_id  = "222222222222",
            build     = "main",
            packages  = list(`omnipath-build` = "ddd", `omnipath-utils` = "bbb",
                             `omnipath-resources` = "ccc"),
            resources = list(), partial_build = FALSE
        )

        dep3 <- list(name = "dev3", url = "dev3.omnipathdb.org",
                     db_host = "localhost", db_port = 5403L,
                     db_name = "omnipath", build_id = m3$build_id,
                     build_manifest_source = "table")
        dep4 <- list(name = "dev4", url = "dev4.omnipathdb.org",
                     db_host = "localhost", db_port = 5404L,
                     db_name = "omnipath", build_id = m4$build_id,
                     build_manifest_source = "table")

        dir.create("out")
        writeLines("dummy", "out/multi.pdf")

        write_sidecar(
            artifact_id   = "database-content",
            artifact_path = "out/multi.pdf",
            deployments   = list(dep3, dep4),
            manifests     = list(m3, m4),
            script_path   = "figures/database-content/build.R"
        )

        loaded <- jsonlite::fromJSON(
            "out/multi.pdf.provenance.json", simplifyVector = FALSE
        )
        expect_equal(length(loaded$deployments), 2L)
        expect_equal(
            vapply(loaded$deployments, function(d) d$name, character(1L)),
            c("dev3", "dev4")
        )
        # package_commits is the union: omnipath-build present from m3
        # AND m4 (m3's value wins for duplicates).
        expect_equal(loaded$package_commits$`omnipath-build`, "aaa")
        expect_equal(loaded$package_commits$`omnipath-resources`, "ccc")
    })
})

test_that("write_sidecar requires at least one deployment + manifest", {

    withr::with_tempdir({
        expect_error(
            write_sidecar(
                artifact_id   = "test",
                artifact_path = "out.pdf",
                deployments   = list(),
                manifests     = list(),
                script_path   = "x.R"
            ),
            "at least one deployment"
        )
        expect_error(
            write_sidecar(
                artifact_id   = "test",
                artifact_path = "out.pdf",
                deployments   = list(list(
                    name = "dev3", url = "", db_host = "",
                    db_port = 5403L, db_name = "x", build_id = "abcdef012345"
                )),
                manifests     = list(),
                script_path   = "x.R"
            ),
            "at least one manifest"
        )
    })
})

test_that("query_record carries the deployment label from pg_query_panel", {

    rows <- tibble::tibble(x = 1:3)
    attr(rows, "sql")         <- "SELECT generate_series(1,3)"
    attr(rows, "result_hash") <- "abcdef012345"
    attr(rows, "deployment")  <- "dev4"

    rec <- query_record(rows)
    expect_equal(rec$sql, "SELECT generate_series(1,3)")
    expect_equal(rec$row_count, 3L)
    expect_equal(rec$result_hash, "abcdef012345")
    expect_equal(rec$deployment, "dev4")
})

test_that("query_record omits deployment when the attribute is missing", {

    rows <- tibble::tibble(x = 1:3)
    attr(rows, "sql")         <- "SELECT 1"
    attr(rows, "result_hash") <- "deadbeef0123"

    rec <- query_record(rows)
    expect_null(rec$deployment)
})
