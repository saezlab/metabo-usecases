test_that("fig01 plot helpers record query attributes through pg_query", {

    # Synthesize a tibble with the attribute shape pg_query produces;
    # the sidecar writer consumes those without touching Postgres.

    rows <- tibble::tibble(
        resource   = c("signor", "chebi"),
        n_entities = c(99L, 42L)
    )
    attr(rows, "sql")         <- "SELECT ..."
    attr(rows, "result_hash") <- "0123456789ab"

    rec <- query_record(rows)
    expect_equal(rec$row_count, 2L)
    expect_equal(rec$result_hash, "0123456789ab")
    expect_equal(rec$sql, "SELECT ...")
})


test_that("category_colour fails fast when a resource is unregistered", {

    rows <- tibble::tibble(
        resource   = "brand_new_resource",
        n_entities = 1L
    )

    expect_error(
        plot_entities_by_resource(rows, width_mm = 89L),
        "missing colour entries"
    )
})


test_that("renderer with registered colours produces a valid ggplot", {

    register_category_colours(
        "resources",
        c(my_resource = "#006384")
    )
    withr::defer(metabo.figures:::init_category_colours())

    rows <- tibble::tibble(
        resource   = "my_resource",
        n_entities = 1L
    )
    expect_s3_class(
        plot_entities_by_resource(rows, width_mm = 89L),
        "ggplot"
    )
})
