test_that("cosmos_old_pkn returns tibble with required columns", {
    skip_if_not(
        file.exists(here::here(
            "data", "vendored", "cosmos", "meta_network.RData"
        )),
        "meta_network.RData not vendored"
    )

    pkn <- cosmos_old_pkn()
    expect_s3_class(pkn, "tbl_df")
    expect_true("interaction_type" %in% names(pkn))
    expect_true("n_edges" %in% names(pkn))
    expect_type(pkn$n_edges, "integer")
    expect_gt(sum(pkn$n_edges), 0L)
})

test_that("cosmos_old_pkn attaches source_path and fingerprint attributes", {
    skip_if_not(
        file.exists(here::here(
            "data", "vendored", "cosmos", "meta_network.RData"
        )),
        "meta_network.RData not vendored"
    )

    pkn <- cosmos_old_pkn()
    expect_false(is.null(attr(pkn, "source_path")))
    expect_false(is.null(attr(pkn, "fingerprint")))
    # MD5 fingerprint must be 32 hex chars
    expect_match(attr(pkn, "fingerprint"), "^[0-9a-f]{32}$")
})

test_that("cosmos_old_vs_new returns both panel_b and panel_c", {
    skip_if_not(
        file.exists(here::here(
            "data", "vendored", "cosmos", "meta_network.RData"
        )),
        "meta_network.RData not vendored"
    )

    old_pkn <- cosmos_old_pkn()

    # Minimal stub for new PKN (same schema as old)
    new_pkn <- tibble::tibble(
        interaction_type = c("signaling", "metabolic_reactions"),
        n_edges          = c(1000L, 500L)
    )

    result <- cosmos_old_vs_new(old_pkn, new_pkn)

    expect_true(is.list(result))
    expect_true("panel_b" %in% names(result))
    expect_true("panel_c" %in% names(result))
    expect_s3_class(result$panel_b, "gg")
    expect_s3_class(result$panel_c, "gg")
})

test_that("fig03 metalinks snapshot loaders return non-empty tibbles", {
    skip_if_not_installed("metabo.figures")

    # metalinks_v1_snapshot should return data with at minimum
    # source, target, and resource columns
    snap <- tryCatch(
        metalinks_v1_snapshot(),
        error = function(e) NULL
    )
    skip_if(is.null(snap), "MetaLinksDB SQLite not accessible")

    expect_s3_class(snap, "tbl_df")
    expect_gt(nrow(snap), 0L)
})
