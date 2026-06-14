test_that("Panel A digest emits all four files + sidecar (FR-043g, SC-012)", {

    skip_if_not_installed("withr")

    fixture <- digest_fixture_default()

    withr::with_tempdir({

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


# ─── Fixture helpers ──────────────────────────────────────────────────────


#' Default canned-response fixture for a working dev5 digest
#'
#' Returns a nested list keyed by section that
#' \code{\link{digest_fixture_pg_query_panel}} dispatches on. Values
#' are tibbles shaped exactly as the corresponding query returns.
#'
#' @noRd
digest_fixture_default <- function() {

    list(
        manifest = list(
            build_id      = "abc123def456",
            build_kind    = "main",
            built_at      = "2026-06-14T12:00:00+02:00",
            partial_build = FALSE
        ),
        entities = list(
            chemical = tibble::tibble(
                chemical_class = c(
                    "metabolite", "lipid", "drug", "food", "xenobiotic"
                ),
                n_entities     = c(50000L, 12000L, 8000L, 4000L, 0L)
            ),
            pog = tibble::tibble(
                n_genes               = 75000L,
                n_unresolved_proteins = 5000L,
                n_protein_or_gene_union = 80000L
            )
        ),
        mpi = list(
            totals = tibble::tibble(
                n_mpi_relations    = 250000L,
                n_metabolites      = 8000L,
                n_proteins_or_genes = 12000L
            ),
            transporters = tibble::tibble(n_transporters = 1500L),
            receptors    = tibble::tibble(n_receptors = 2000L)
        ),
        interactions = list(
            totals = tibble::tibble(
                n_interactions     = 5000000L,
                n_proteins_or_genes = 22000L
            ),
            pathway_preferred = tibble::tibble(n_pathways = 1800L),
            pathway_fallback  = tibble::tibble(n_pathways = 1900L),
            reaction_preferred = tibble::tibble(n_reactions = 800L),
            reaction_fallback  = tibble::tibble(n_reactions = 900L)
        ),
        structures = list(
            levels = tibble::tibble(
                specificity_level = structural_specificity_levels(),
                n_structures      = c(60000L, 5000L, 3000L, 1500L, 800L, 0L)
            ),
            inchikey = tibble::tibble(
                n_skeletons     = 45000L,
                n_full_inchikey = 70000L
            ),
            ramp = tibble::tibble(n_ramp_conflicts = 320L)
        ),
        annotation = list(
            totals = tibble::tibble(
                n_annotation_records = 9000000L,
                n_organisms          = 350L
            ),
            slot3 = list(
                diseases = tibble::tibble(n_diseases = 15000L)
            ),
            localization = tibble::tibble(n_localizations = 4000L),
            pathway      = tibble::tibble(n_pathways = 1900L)
        ),
        resources = list(
            `1` = tibble::tibble(
                resource_name = c("hmdb", "chebi", "chembl"),
                n_rows = c(50000L, 100000L, 200000L)
            ),
            `2` = tibble::tibble(
                resource_name = c("guidetopharma", "tcdb"),
                n_rows = c(2000L, 1500L)
            ),
            `3` = tibble::tibble(
                resource_name = c("kegg", "reactome", "signor"),
                n_rows = c(900000L, 800000L, 100000L)
            ),
            `4` = tibble::tibble(
                resource_name = c("chembl", "lipidmaps"),
                n_rows = c(45000L, 8000L)
            ),
            `5` = tibble::tibble(
                resource_name = c("mondo", "hpo", "go"),
                n_rows = c(15000L, 25000L, 80000L)
            )
        )
    )
}


#' Build the pg_query_panel mock that dispatches on (facet, sql)
#'
#' Returns a function callable as a drop-in replacement for
#' \code{pg_query_panel}. Routes by facet → fixture branch and
#' inspects the SQL string when a facet maps to multiple queries
#' (totals vs transporter vs receptor vs pathway preferred vs
#' fallback vs reaction preferred vs fallback etc.).
#'
#' @noRd
digest_fixture_pg_query_panel <- function(fixture) {

    function(panel_id, sql, ..., facet = NULL) {

        rows <- digest_fixture_dispatch(fixture, facet, sql)

        attr(rows, "deployment")  <- "dev5"
        attr(rows, "sql")         <- trimws(sql)
        attr(rows, "result_hash") <- substr(
            digest::digest(rows, algo = "sha256"), 1L, 12L
        )
        rows
    }
}


#' Map (facet, sql) to a fixture row set
#'
#' @noRd
digest_fixture_dispatch <- function(fixture, facet, sql) {

    if (identical(facet, "panel_a_stats_entities")) {
        if (grepl("chemical_class_id", sql, fixed = TRUE)) {
            return(fixture$entities$chemical)
        }
        return(fixture$entities$pog)
    }
    if (identical(facet, "panel_a_stats_mpi")) {
        if (grepl("n_transporters", sql, fixed = TRUE)) {
            return(fixture$mpi$transporters)
        }
        if (grepl("n_receptors", sql, fixed = TRUE)) {
            return(fixture$mpi$receptors)
        }
        return(fixture$mpi$totals)
    }
    if (identical(facet, "panel_a_stats_interactions")) {
        if (grepl("n_pathways", sql, fixed = TRUE)) {
            if (grepl("entity_ontology_term", sql, fixed = TRUE)) {
                return(fixture$interactions$pathway_fallback)
            }
            return(fixture$interactions$pathway_preferred)
        }
        if (grepl("n_reactions", sql, fixed = TRUE)) {
            if (grepl("data_source ds", sql, fixed = TRUE) &&
                grepl("'brenda'", sql, fixed = TRUE)) {
                return(fixture$interactions$reaction_fallback)
            }
            return(fixture$interactions$reaction_preferred)
        }
        return(fixture$interactions$totals)
    }
    if (identical(facet, "panel_a_stats_structures")) {
        if (grepl("metabo_ramp_inchikey_conflict", sql, fixed = TRUE)) {
            return(fixture$structures$ramp)
        }
        if (grepl("SUBSTRING", sql, fixed = TRUE) ||
            grepl("n_skeletons", sql, fixed = TRUE)) {
            return(fixture$structures$inchikey)
        }
        return(fixture$structures$levels)
    }
    if (identical(facet, "panel_a_stats_annotation")) {
        if (grepl("n_diseases", sql, fixed = TRUE)) {
            return(fixture$annotation$slot3$diseases)
        }
        if (grepl("n_localizations", sql, fixed = TRUE)) {
            return(fixture$annotation$localization)
        }
        if (grepl("n_pathways", sql, fixed = TRUE)) {
            return(fixture$annotation$pathway)
        }
        return(fixture$annotation$totals)
    }

    section_key <- sub(
        "^panel_a_stats_section_([0-9])_resources$",
        "\\1",
        facet
    )
    if (!is.na(section_key) && nzchar(section_key)) {
        return(fixture$resources[[section_key]])
    }

    rlang::abort(sprintf(
        "Test fixture is missing a row set for facet='%s'", facet
    ))
}


#' Write a minimal digest-config.yaml for tests
#'
#' @noRd
write_temp_digest_config <- function(dir, name) {
    path <- file.path(dir, name)
    writeLines(c(
        "definitions:",
        "  transporter:",
        "    resources: [tcdb, slctables]",
        "    uniprot_keywords: [KW-0813]",
        "  receptor:",
        "    resources: [guidetopharma, cellphonedb]",
        "    uniprot_keywords: [KW-0675]",
        "  pathway:",
        "    preferred: Pathway:OM:0014",
        "    fallback: annotation",
        "  reaction:",
        "    preferred: Reaction:OM:0015",
        "    fallback: predicate-classification",
        "  annotation_slot3: diseases",
        "runtime:",
        "  structures:",
        "    inchikey_column: inchikey",
        "  pathway:",
        "    min_rows_for_preferred: 10",
        "  reaction:",
        "    min_rows_for_preferred: 10",
        "  annotation:",
        "    disease_min_rows: 100",
        "    slot3_fallback_order: [diseases, phenotypes, tissues]"
    ), path)
    path
}
