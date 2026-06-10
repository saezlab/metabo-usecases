test_that("fig03 and fig04 use the same category_colours registry", {
    # Both figures render interaction_type categories. The style module
    # must supply a single named vector so the same category renders
    # identically across both figures (FR-022, SC-004).

    required_categories <- c(
        "metabolic_reactions",
        "signaling",
        "transport",
        "ligand_receptor",
        "drug_target",
        "tf_target"
    )

    colours <- metabo.figures::category_colours()
    interaction_colours <- colours[["interaction_types"]]

    missing <- setdiff(required_categories, names(interaction_colours))
    expect_length(
        missing, 0L,
        label = paste(
            "interaction_type categories not registered in palette:",
            paste(missing, collapse = ", ")
        )
    )
})

test_that("no inline colour overrides in fig03 build script", {
    build_path <- here::here(
        "figures", "fig03-metalinks-versions", "build.R"
    )
    skip_if_not(file.exists(build_path), "fig03 build.R not present")

    content <- readLines(build_path)
    # Inline hex literals are a constitution violation (FR-017a)
    hex_lines <- grep('#[0-9A-Fa-f]{6}\\b', content, value = TRUE)
    expect_length(
        hex_lines, 0L,
        label = paste(
            "inline hex colours found in fig03/build.R — use style module:",
            paste(hex_lines, collapse = "; ")
        )
    )
})

test_that("no inline colour overrides in fig04 build script", {
    build_path <- here::here(
        "figures", "fig04-cosmos-pkn", "build.R"
    )
    skip_if_not(file.exists(build_path), "fig04 build.R not present")

    content <- readLines(build_path)
    hex_lines <- grep('#[0-9A-Fa-f]{6}\\b', content, value = TRUE)
    expect_length(
        hex_lines, 0L,
        label = paste(
            "inline hex colours found in fig04/build.R — use style module:",
            paste(hex_lines, collapse = "; ")
        )
    )
})
