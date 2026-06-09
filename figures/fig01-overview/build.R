# figures/fig01-overview/build.R
#
# Orchestrates the full Figure 1 build: queries the OmniPath Postgres,
# renders panels B–F via ggplot, vendors the manual architecture asset
# at the Panel A position (FR-005, FR-005a — SHA-256 fingerprint
# verified against the inst/extdata/manual/architecture/README.md pin
# before include), composes via tex/compose_fig01.tex, and writes the
# provenance sidecar.
#
# Sourced by rebuild.R; safe to source standalone too.

suppressPackageStartupMessages({
    library(metabo.figures)
    library(ggplot2)
})

setup_pipeline_log("build:fig01-overview")
set.seed(pipeline_seed())

out_dir <- "figures/fig01-overview/out"
fs::dir_create(out_dir)

# ---- Deployment + manifest -----------------------------------------------
#
# Phase-3 (US1) MVP: only dev3 is touched. T014a's per-panel resolver
# routes the Structures facet / FR-007e / FR-007f panels to dev4 — when
# those renderers land they add a `dep4 <- deployment_provenance("dev4")`
# call here and list both in the `deployments`/`manifests` args below.

logger::log_info("Resolving dev3 deployment + reading build_manifest")
dep3 <- deployment_provenance("dev3")
con  <- pg_connect_panel("dev3")

# ---- Data ------------------------------------------------------------------

logger::log_info("Querying Panel B (entities_by_resource)")
data_b <- entities_by_resource(con)

logger::log_info("Querying Panel C (interactions_by_resource)")
data_c <- interactions_by_resource(con)

logger::log_info("Querying Panel D (interactions_by_type)")
data_d <- interactions_by_type(con)

logger::log_info("Querying Panel E (annotation_classes_by_resource)")
data_e <- annotation_classes_by_resource(con)

logger::log_info("Querying Panel F (ontology_terms_by_ontology)")
data_f <- ontology_terms_by_ontology(con)

queries <- list(
    query_record(data_b),
    query_record(data_c),
    query_record(data_d),
    query_record(data_e),
    query_record(data_f)
)

# Register the categories that appear in this snapshot. The
# assertive accessor catches anything else.
all_resources <- unique(c(
    data_b$resource, data_c$resource, data_e$resource
))
if (length(all_resources) > 0L) {
    register_category_colours(
        "resources",
        setNames(
            palette_n(length(all_resources), unknown = FALSE),
            all_resources
        )
    )
}

all_interactions <- unique(data_d$interaction_type)
new_interactions <- setdiff(
    all_interactions,
    registered_category_values("interaction_types")
)
if (length(new_interactions) > 0L) {
    register_category_colours(
        "interaction_types",
        setNames(
            palette_n(
                length(new_interactions),
                unknown = FALSE
            ),
            new_interactions
        )
    )
}

# ---- Panels B–F ------------------------------------------------------------

panels <- list(
    panelB = plot_entities_by_resource(data_b, width_mm = 89L),
    panelC = plot_interactions_by_resource(data_c, width_mm = 89L),
    panelD = plot_interactions_by_type(data_d, width_mm = 89L),
    panelE = plot_annotation_classes_by_resource(data_e, width_mm = 89L),
    panelF = plot_ontology_terms_by_ontology(data_f, width_mm = 89L)
)

for (name in names(panels)) {
    out_pdf <- file.path(out_dir, sprintf("%s.pdf", name))
    out_svg <- file.path(out_dir, sprintf("%s.svg", name))

    ggsave(out_pdf, panels[[name]], width = 89, height = 60, units = "mm")
    ggsave(out_svg, panels[[name]], width = 89, height = 60, units = "mm")

    gate <- readability_check(panels[[name]], width_mm = 89L)
    logger::log_info(
        "{name}: {ifelse(gate$pass, 'PASS', 'FAIL')} ",
        "({gate$summary})"
    )
}

# ---- Panel A — vendored architecture asset (FR-005, FR-005a) ---------------

architecture_dir <- "inst/extdata/manual/architecture"
architecture_pdf <- file.path(
    architecture_dir, "omnipath-architecture-new2026.pdf"
)
architecture_readme <- file.path(architecture_dir, "README.md")

if (!file.exists(architecture_pdf)) {
    rlang::abort(sprintf(
        "Vendored architecture asset missing: %s. Re-vendor per %s.",
        architecture_pdf, architecture_readme
    ))
}

readme_pin <- stringr::str_match(
    paste(readLines(architecture_readme), collapse = "\n"),
    "pdf_sha256.*?`([0-9a-f]{64})`"
)[1L, 2L]
if (is.na(readme_pin)) {
    rlang::abort(sprintf(
        "Could not parse pdf_sha256 pin from %s", architecture_readme
    ))
}

architecture_sha256 <- digest::digest(
    file = architecture_pdf, algo = "sha256"
)
if (!identical(architecture_sha256, readme_pin)) {
    rlang::abort(sprintf(paste0(
        "Architecture asset SHA-256 mismatch (FR-005a). ",
        "Expected (README pin): %s\n",
        "Got      (on disk):    %s\n",
        "Manual asset changed — update %s and re-record fingerprint."
    ), readme_pin, architecture_sha256, architecture_readme))
}

file.copy(
    architecture_pdf,
    file.path(out_dir, "panelA.pdf"),
    overwrite = TRUE
)
logger::log_info(
    "Panel A vendored from {architecture_pdf} (sha256={substr(architecture_sha256, 1L, 12L)})"
)

# ---- Composite -------------------------------------------------------------

logger::log_info("Composing Figure 1")
compose_mixed_source(
    mode = "pdf",
    spec = list(
        template  = "tex/compose_fig01.tex",
        output    = file.path(out_dir, "fig01-overview.pdf"),
        component = "compose:fig01-overview"
    ),
    work_dir = out_dir
)

# ---- Provenance sidecar ----------------------------------------------------

write_sidecar(
    artifact_id    = "fig01-overview",
    artifact_path  = file.path(out_dir, "fig01-overview.pdf"),
    deployments    = list(dep3$deployment),
    manifests      = list(dep3$manifest),
    script_path    = "figures/fig01-overview/build.R",
    queries        = queries,
    external_inputs = list(
        list(
            kind        = "architecture-asset",
            path        = architecture_pdf,
            source      = paste0(
                "Inkscape source on the contributor's machine; ",
                "see ", architecture_readme
            ),
            fingerprint = architecture_sha256
        )
    ),
    parameters     = list(width_mm = 180L),
    seed           = pipeline_seed()
)

pg_close_all_panel()
logger::log_info("fig01-overview complete")
