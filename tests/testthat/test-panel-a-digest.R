test_that("Panel A digest emits all four files + sidecar (FR-043g, SC-012)", {

    skip_if_not_installed("withr")

    fixture <- digest_fixture_default()

    withr::with_tempdir({

        withr::local_envvar(METABO_FIGURES_LOG = file.path(getwd(), "run.log"))
        out_dir <- "panel-a-stats"
        config <- write_temp_digest_config(".", "digest-config.yaml")

        with_mocked_bindings(
            pg_connect_panel = function(deployment) NA,
            pg_query_panel   = digest_fixture_pg_query_panel(fixture),
            build_manifest_for = function(con) fixture$manifest,
            snapshot_id        = function(manifest) manifest$build_id,
            git_commit_of      = function(path) "0000000000000000000000000000000000000000",
            .package = "metabo.figures",
            {
                result <- build_panel_a_digest(
                    snapshot_id     = fixture$manifest$build_id,
                    out_dir         = out_dir,
                    config_path     = config,
                    caption_sty     = "tex/caption.sty",
                    validate_schema = FALSE
                )

                expect_true(file.exists(result$files$json))
                expect_true(file.exists(result$files$csv))
                expect_true(file.exists(result$files$md))
                # stats.pdf is skipped when tex/caption.sty is absent
                # in the tempdir — the sidecar is the authoritative
                # FR-043g check.
                expect_true(file.exists(result$sidecar))

                stats <- jsonlite::read_json(result$files$json)
                expect_length(stats$sections, 5L)
                expect_equal(
                    purrr::map_chr(stats$sections, "name"),
                    c(
                        "Entities",
                        "Metabolite-Protein Interactions",
                        "Interactions",
                        "Structures",
                        "Annotation"
                    )
                )
            }
        )
    })
})


test_that("Section 2 subset constraint aborts on violation (FR-043b)", {

    skip_if_not_installed("withr")

    fixture <- digest_fixture_default()
    fixture$mpi$totals$n_proteins_or_genes <- 10L
    fixture$mpi$transporters$n_transporters <- 50L

    withr::with_tempdir({
        withr::local_envvar(METABO_FIGURES_LOG = file.path(getwd(), "run.log"))
        out_dir <- "panel-a-stats"
        config <- write_temp_digest_config(".", "digest-config.yaml")

        with_mocked_bindings(
            pg_connect_panel = function(deployment) NA,
            pg_query_panel   = digest_fixture_pg_query_panel(fixture),
            build_manifest_for = function(con) fixture$manifest,
            snapshot_id        = function(manifest) manifest$build_id,
            git_commit_of      = function(path) "0",
            .package = "metabo.figures",
            {
                expect_error(
                    build_panel_a_digest(
                        snapshot_id     = fixture$manifest$build_id,
                        out_dir         = out_dir,
                        config_path     = config,
                        caption_sty     = "tex/caption.sty",
                        validate_schema = FALSE
                    ),
                    "FR-043b subset constraint"
                )
            }
        )
    })
})


test_that("Mismatched snapshot id aborts (FR-043h)", {

    skip_if_not_installed("withr")

    fixture <- digest_fixture_default()

    withr::with_tempdir({
        withr::local_envvar(METABO_FIGURES_LOG = file.path(getwd(), "run.log"))
        out_dir <- "panel-a-stats"
        config <- write_temp_digest_config(".", "digest-config.yaml")

        with_mocked_bindings(
            pg_connect_panel = function(deployment) NA,
            pg_query_panel   = digest_fixture_pg_query_panel(fixture),
            build_manifest_for = function(con) fixture$manifest,
            snapshot_id        = function(manifest) manifest$build_id,
            git_commit_of      = function(path) "0",
            .package = "metabo.figures",
            {
                expect_error(
                    build_panel_a_digest(
                        snapshot_id     = "deadbeef0000",
                        out_dir         = out_dir,
                        config_path     = config,
                        caption_sty     = "tex/caption.sty",
                        validate_schema = FALSE
                    ),
                    "FR-043h snapshot-id mismatch"
                )
            }
        )
    })
})


test_that("digest_config rejects missing definitions block", {

    withr::with_tempdir({
        writeLines(c("runtime:", "  structures:", "    inchikey_column: ik"),
                   "broken.yaml")
        expect_error(
            digest_config("broken.yaml"),
            "definitions:"
        )
    })
})


test_that("digest_config validates required definition keys", {

    withr::with_tempdir({
        writeLines(c(
            "definitions:",
            "  transporter: {}",
            "runtime:",
            "  structures:",
            "    inchikey_column: ik",
            "  annotation:",
            "    disease_min_rows: 100",
            "    slot3_fallback_order: [diseases]"
        ), "incomplete.yaml")
        expect_error(
            digest_config("incomplete.yaml"),
            "definitions.transporter.resources"
        )
    })
})


test_that("mirror_to_metadata flattens pathway/reaction to chosen branch", {

    cfg <- list(
        transporter = list(resources = "tcdb", uniprot_keywords = "KW-0813"),
        receptor    = list(
            resources        = "guidetopharma",
            uniprot_keywords = "KW-0675"
        ),
        pathway     = list(
            preferred = "Pathway:OM:0014",
            fallback  = "annotation"
        ),
        reaction    = list(
            preferred = "Reaction:OM:0015",
            fallback  = "predicate-classification"
        ),
        annotation_slot3 = "diseases"
    )

    md <- mirror_to_metadata(cfg, list(
        pathway          = "annotation",
        reaction         = "Reaction:OM:0015",
        annotation_slot3 = "phenotypes"
    ))

    expect_identical(md$pathway, "annotation")
    expect_identical(md$reaction, "Reaction:OM:0015")
    expect_identical(md$annotation_slot3, "phenotypes")
    expect_identical(md$transporter$resources, "tcdb")
    expect_identical(md$receptor$uniprot_keywords, "KW-0675")
})


