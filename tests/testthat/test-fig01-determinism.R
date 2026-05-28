test_that("fig01 quantitative panels are deterministic across re-renders", {

    skip_on_ci()
    skip_if_not_installed("ggplot2")

    # Build a fixture-driven panel deterministically twice.
    fake_data <- tibble::tibble(
        resource   = c("signor", "chebi", "omnipath_metabo"),
        n_entities = c(20L, 10L, 5L)
    )

    register_category_colours(
        "resources",
        c(signor = "#006384", chebi = "#9F0162",
          omnipath_metabo = "#FEAF16")
    )
    withr::defer(init_category_colours())

    set.seed(pipeline_seed())
    p1 <- plot_entities_by_resource(fake_data, width_mm = 89L)

    set.seed(pipeline_seed())
    p2 <- plot_entities_by_resource(fake_data, width_mm = 89L)

    out1 <- tempfile(fileext = ".svg")
    out2 <- tempfile(fileext = ".svg")
    ggplot2::ggsave(out1, p1, width = 89, height = 60, units = "mm")
    ggplot2::ggsave(out2, p2, width = 89, height = 60, units = "mm")

    expect_equal(
        digest::digest(file = out1),
        digest::digest(file = out2)
    )
})
