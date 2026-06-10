test_that("lead palette is in canonical order with grey reserved", {

    lead <- palette_lead()

    expect_length(lead, 9L)
    expect_equal(unname(lead[[1L]]), "#006384")
    expect_equal(unname(lead[[9L]]), "#BEBEBE")
    expect_equal(names(lead)[[9L]], "unknown")
})

test_that("palette_n(1) returns the unknown grey (FR-020)", {

    expect_equal(palette_n(1L), "#BEBEBE")
})

test_that("palette_n(4) returns the first four lead colours", {

    expect_equal(
        palette_n(4L),
        c("#006384", "#9F0162", "#FEAF16", "#BBCC33")
    )
})

test_that("palette_n with unknown swaps in grey for the last slot", {

    out <- palette_n(4L, unknown = TRUE)
    expect_equal(out[[4L]], "#BEBEBE")
    expect_equal(out[seq_len(3L)], c("#006384", "#9F0162", "#FEAF16"))
})

test_that("palette_n refuses to reach the reserved unknown slot", {

    expect_error(palette_n(9L), "reserved unknown slot")
})

test_that("rwth palette parses to hex codes", {

    rwth <- palette_rwth()

    expect_true(length(rwth) > 0L)
    expect_match(unname(rwth), "^#[0-9A-F]{6}$")
})

test_that("category_colour returns the registered hex code", {

    withr::defer(metabo.figures:::init_category_colours())
    register_category_colours(
        "interaction_types",
        c(signaling = "#006384", transport = "#D03293")
    )

    expect_equal(
        category_colour("interaction_types", "signaling"),
        "#006384"
    )
    expect_equal(
        category_colour(
            "interaction_types",
            c("signaling", "transport")
        ),
        c("#006384", "#D03293")
    )
})

test_that("category_colour fails fast on unknown category", {

    expect_error(
        category_colour("nope", "x"),
        "not registered"
    )
})

test_that("category_colour fails fast on unknown value (spec Edge Case)", {

    expect_error(
        category_colour("interaction_types", "fictional_type"),
        "missing colour entries"
    )
})

test_that("register_category_colours appends but refuses overwrites", {

    withr::defer(metabo.figures:::init_category_colours())

    register_category_colours(
        "resources",
        c(signor = "#006384", chebi = "#9F0162")
    )
    expect_equal(
        category_colour("resources", "signor"),
        "#006384"
    )

    # Conflicting re-register fails.
    expect_error(
        register_category_colours(
            "resources",
            c(signor = "#FF0000")
        ),
        "Refusing to overwrite"
    )

    # Same value re-register is a no-op.
    expect_silent(
        register_category_colours(
            "resources",
            c(signor = "#006384")
        )
    )
})
