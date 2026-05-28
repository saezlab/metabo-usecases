test_that("compliant ggplot passes readability_check at 89 mm", {

    plot <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
        ggplot2::geom_point() +
        theme_bw_metabo(width_mm = 89)

    result <- readability_check(plot, width_mm = 89)

    expect_true(result$pass)
    expect_equal(nrow(result$violations), 0L)
})

test_that("compliant ggplot passes readability_check at 180 mm", {

    plot <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
        ggplot2::geom_point() +
        theme_bw_metabo(width_mm = 180)

    expect_true(readability_check(plot, width_mm = 180)$pass)
})

test_that("plot with too-small axis title fails the gate", {

    plot <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
        ggplot2::geom_point() +
        theme_bw_metabo(width_mm = 89) +
        ggplot2::theme(
            axis.title = ggplot2::element_text(size = 4)
        )

    result <- readability_check(plot, width_mm = 89)

    expect_false(result$pass)
    expect_true("axis_label" %in% result$violations$role)
})

test_that("plot with too-small tick labels fails the gate", {

    plot <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
        ggplot2::geom_point() +
        theme_bw_metabo(width_mm = 89) +
        ggplot2::theme(
            axis.text = ggplot2::element_text(size = 3)
        )

    result <- readability_check(plot, width_mm = 89)

    expect_false(result$pass)
    expect_true("tick" %in% result$violations$role)
})
