# figures/fig01-overview/build.R
#
# Orchestrates the full Figure 1 build: queries the OmniPath Postgres,
# renders panels B–F via ggplot, compiles Panel A from TikZ via
# xelatex, composes everything via tex/compose_fig01.tex, and writes
# the provenance sidecar.
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

deployment <- load_connection()
con <- pg_connect()
on.exit(DBI::dbDisconnect(con), add = TRUE)

main_manifest <- build_manifest_for(con, "main")
sid_main <- write_manifest(main_manifest)
logger::log_info("Using main snapshot {sid_main}")

# ---- Data ------------------------------------------------------------------

logger::log_info("Querying database content")
data_b <- entities_by_resource(con)
data_c <- interactions_by_resource(con)
data_d <- interactions_by_type(con)
data_e <- annotation_classes_by_resource(con)
data_f <- ontology_terms_by_ontology(con)

queries <- list(
    query_record(data_b),
    query_record(data_c),
    query_record(data_d),
    query_record(data_e),
    query_record(data_f)
)

# Register the resources that appear in this snapshot. We register
# only what we see — the assertive accessor catches anything else.
all_resources <- unique(c(
    data_b$resource, data_c$resource, data_e$resource
))
register_category_colours(
    "resources",
    setNames(
        palette_n(length(all_resources), unknown = FALSE),
        all_resources
    )
)

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

# ---- Panel A (TikZ + xelatex) ----------------------------------------------

logger::log_info("Compiling Panel A (TikZ → xelatex)")
sh <- paste0(
    "source lib/log.sh && ",
    "log_xelatex_capture xelatex:fig01a xelatex ",
    "-interaction=nonstopmode ",
    "-output-directory=", shQuote(out_dir), " ",
    "tikz/fig1a_architecture.tex"
)
rc <- system2("bash", c("-c", sh))
if (rc != 0L) {
    rlang::abort(sprintf("Panel A xelatex failed (rc = %d)", rc))
}
file.copy(
    file.path(out_dir, "fig1a_architecture.pdf"),
    file.path(out_dir, "panelA.pdf"),
    overwrite = TRUE
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
    deployment     = deployment,
    manifests      = list(main_manifest),
    script_path    = "figures/fig01-overview/build.R",
    queries        = queries,
    external_inputs = list(
        list(
            kind        = "vendored-asset",
            path        = "inst/extdata/assets/postgres.svg",
            source      = "PostgreSQL trademark, see SOURCES.md",
            fingerprint = substr(
                digest::digest(file = "inst/extdata/assets/postgres.svg"),
                1L, 12L
            )
        ),
        list(
            kind        = "vendored-asset",
            path        = "inst/extdata/assets/rdkit.png",
            source      = "RDKit project (BSD)",
            fingerprint = substr(
                digest::digest(file = "inst/extdata/assets/rdkit.png"),
                1L, 12L
            )
        )
    ),
    parameters     = list(width_mm = 180L),
    seed           = pipeline_seed()
)

logger::log_info("fig01-overview complete")
