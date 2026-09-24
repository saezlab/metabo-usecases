# Regression tests for the 2026-06-18 Figure 5 composite merge
# (spec FR-011 family, Session 2026-06-18; tasks T056c-T056g).
#
# Guards the data-correctness fixes the merge made normative and the
# re-addition of the old-vs-COSMOS+ comparison as a composite panel.

cosmos_csvs_present <- function() {
    all(file.exists(c(
        here::here("data", "vendored", "cosmos", "cosmos_plus_human.csv"),
        here::here("data", "vendored", "cosmos", "cosmos_plus_mouse.csv")
    )))
}

test_that("comparison panel folds catalysis into metabolic_reactions", {
    # The interaction-type lookup fix: COSMOS+ `catalysis` (incl. KEGG
    # catalysis) must be aliased to metabolic_reactions so it aligns with the
    # old PKN's metabolic-reactions bar — it must not appear as its own axis
    # category.
    old_pkn <- tibble::tibble(
        interaction_type = "metabolic_reactions",
        n_edges          = 100L
    )
    new_by_type_species <- tibble::tibble(
        interaction_type = c("catalysis", "signaling", "catalysis"),
        species          = c("human", "human", "mouse"),
        n_interactions   = c(50L, 30L, 40L)
    )

    p <- fig04_cosmos_comparison_panel(old_pkn, new_by_type_species)
    expect_s3_class(p, "gg")

    types <- as.character(unique(p$data$interaction_type))
    expect_false("catalysis" %in% types)
    expect_true("metabolic_reactions" %in% types)
})

test_that("comparison panel keeps three species groups", {
    old_pkn <- tibble::tibble(
        interaction_type = "metabolic_reactions",
        n_edges          = 100L
    )
    new_by_type_species <- tibble::tibble(
        interaction_type = c("signaling", "signaling"),
        species          = c("human", "mouse"),
        n_interactions   = c(30L, 20L)
    )

    p <- fig04_cosmos_comparison_panel(old_pkn, new_by_type_species)
    groups <- as.character(unique(p$data$panel_group))
    expect_setequal(
        groups,
        c("Old COSMOS (human)", "COSMOS+ (human)", "COSMOS+ (mouse)")
    )
})

test_that("cosmos_plus_data folds child compartments into parents", {
    skip_if_not(cosmos_csvs_present(), "COSMOS+ CSVs not vendored")

    cp <- cosmos_plus_data()
    labels <- cp$by_compartment$compartment_name

    # `e`/`eg` must fold into "Cell membrane" and `i` into "Mitochondria";
    # the pre-fold child labels must not survive into the panel data.
    expect_false("Extracellular" %in% labels)
    expect_false("Extracellular (TCDB)" %in% labels)
    expect_false("Mitochondrial intermembrane space" %in% labels)
})

test_that("cosmos_plus_data canonicalises and de-duplicates resources", {
    skip_if_not(cosmos_csvs_present(), "COSMOS+ CSVs not vendored")

    cp <- cosmos_plus_data()

    # by_resource_split is the shape used by composite Panel D.
    expect_true("by_resource_split" %in% names(cp))

    # bare "Recon3D" must have been canonicalised to "GEM:Recon3D".
    expect_false(any(cp$by_resource$resource == "Recon3D"))
    expect_false(any(cp$by_resource_split$resource == "Recon3D"))

    # no semicolon-joined resource string may contain a duplicated token
    # (dedup after alias remapping).
    dup_in_string <- vapply(
        strsplit(cp$by_resource$resource, ";", fixed = TRUE),
        function(tok) any(duplicated(tok)),
        logical(1L)
    )
    expect_false(any(dup_in_string))
})

test_that("composition.yaml declares the four-panel composite", {
    skip_if_not_installed("yaml")
    cfg_path <- here::here(
        "figures", "cosmos-pkn", "composition.yaml"
    )
    skip_if_not(file.exists(cfg_path), "composition.yaml missing")

    cfg <- yaml::read_yaml(cfg_path)
    selected <- Filter(function(p) isTRUE(p$selected), cfg$panels)
    labels <- vapply(selected, function(p) p$label, character(1L))

    # Composite = schematic A + comparison B + compartments C + resources D.
    expect_setequal(labels, c("A", "B", "C", "D"))

    renderers <- vapply(
        selected,
        function(p) if (is.null(p$renderer)) NA_character_ else p$renderer,
        character(1L)
    )
    expect_true("fig04_cosmos_comparison_panel" %in% renderers)

    # MetaLinksDB panel must be standalone, not composited.
    standalone <- if (is.null(cfg$standalone)) list() else cfg$standalone
    standalone_renderers <- vapply(
        standalone,
        function(p) if (is.null(p$renderer)) NA_character_ else p$renderer,
        character(1L)
    )
    expect_true("fig04_metalinks_cosmos_panel" %in% standalone_renderers)
})
