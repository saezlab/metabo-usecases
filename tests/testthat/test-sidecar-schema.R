test_that("write_sidecar emits a file with the contract fields", {

    withr::with_tempdir({

        manifest <- list(
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
            name    = "dev3",
            url     = "dev3.omnipathdb.org",
            db_host = "localhost",
            db_port = 5403L,
            db_name = "omnipath"
        )

        artifact_path <- "out/artifact.pdf"
        dir.create("out")
        writeLines("dummy", artifact_path)

        sc <- write_sidecar(
            artifact_id   = "fig01-overview/panelB",
            artifact_path = artifact_path,
            deployment    = deployment,
            manifests     = list(manifest),
            script_path   = "figures/fig01-overview/build.R",
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
        expect_equal(loaded$artifact_id, "fig01-overview/panelB")
        expect_equal(loaded$deployment$name, "dev3")
        expect_equal(loaded$deployment$db_port, 5403L)
        expect_equal(loaded$snapshot_ids$main, snapshot_id(manifest))
        expect_equal(loaded$package_commits$`omnipath-build`, "aaa1111")
        expect_equal(length(loaded$queries), 1L)
        expect_equal(loaded$queries[[1L]]$sql, "SELECT 1")
        expect_match(loaded$snapshot_ids$main, "^[0-9a-f]{12}$")
    })
})

test_that("write_sidecar requires at least one manifest", {

    withr::with_tempdir({
        expect_error(
            write_sidecar(
                artifact_id   = "test",
                artifact_path = "out.pdf",
                deployment    = list(name = "dev3", url = "", db_host = "",
                                     db_port = 1L, db_name = "x"),
                manifests     = list(),
                script_path   = "x.R"
            ),
            "at least one manifest"
        )
    })
})

test_that("query_record extracts pg_query attributes", {

    rows <- tibble::tibble(x = 1:3)
    attr(rows, "sql")         <- "SELECT generate_series(1,3)"
    attr(rows, "result_hash") <- "abcdef012345"

    rec <- query_record(rows)
    expect_equal(rec$sql, "SELECT generate_series(1,3)")
    expect_equal(rec$row_count, 3L)
    expect_equal(rec$result_hash, "abcdef012345")
})
