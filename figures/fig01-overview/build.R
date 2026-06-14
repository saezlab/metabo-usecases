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

# ---- Deployment + manifests -----------------------------------------------
#
# Post-2026-06-14 dev5 integrated-build promotion: every panel runs
# against dev5 (integrated build — gene-centric entities + stored
# labels including Goslin lipid names + chemical-fallback resolution
# + RDKit-derived structural specificity + RaMP-conflict tables +
# cycle-001 derived family). The per-panel dev3/dev4 split is
# retired; the registry's overrides block is empty.

logger::log_info("Resolving dev5 deployment")
dep5 <- deployment_provenance("dev5")

# ---- Data (cycle-001 derived shapes, dispatched via T014a) ----------------

logger::log_info("Querying Panel B (entities_by_resource)")
data_b <- entities_by_resource()

logger::log_info("Querying Panel C (interactions_by_resource)")
data_c <- interactions_by_resource()

logger::log_info("Querying Panel D (interactions_by_type)")
data_d <- interactions_by_type()

logger::log_info("Querying Panel E (associations_by_resource)")
data_e <- associations_by_resource()

logger::log_info("Querying Panel F (ontology_terms_by_ontology)")
data_f <- ontology_terms_by_ontology()

logger::log_info("Querying Panel G (ramp_conflict_counts on dev5)")
data_g <- ramp_conflict_counts()

# Register the RaMP-conflict reasons in the colour registry on first
# encounter (FR-022, SC-004).
new_reasons <- setdiff(
    unique(data_g$conflict_reason),
    registered_category_values("ramp_conflict_reasons")
)
if (length(new_reasons) > 0L) {
    register_category_colours(
        "ramp_conflict_reasons",
        setNames(
            palette_n(length(new_reasons), unknown = FALSE),
            new_reasons
        )
    )
}

queries <- list(
    query_record(data_b),
    query_record(data_c),
    query_record(data_d),
    query_record(data_e),
    query_record(data_f),
    query_record(data_g)
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

all_interaction_classes <- unique(data_d$interaction_class)
new_interaction_classes <- setdiff(
    all_interaction_classes,
    registered_category_values("interaction_types")
)
if (length(new_interaction_classes) > 0L) {
    register_category_colours(
        "interaction_types",
        setNames(
            palette_n(
                length(new_interaction_classes),
                unknown = FALSE
            ),
            new_interaction_classes
        )
    )
}

# ---- Panels B–G ------------------------------------------------------------

panels <- list(
    panelB = plot_entities_by_resource(data_b, width_mm = 89L),
    panelC = plot_interactions_by_resource(data_c, width_mm = 89L),
    panelD = plot_interactions_by_type(data_d, width_mm = 89L),
    panelE = plot_associations_by_resource(data_e, width_mm = 89L),
    panelF = plot_ontology_terms_by_ontology(data_f, width_mm = 89L),
    panelG = plot_ramp_conflict(data_g, width_mm = 89L)
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

# ---- FR-007a — 6-facet resource overview (stand-alone artifact) -----------
#
# Heavy SQL (~145s end-to-end against dev5); the v1 path uses
# row-scan queries. A follow-up commit will swap this for the
# facet_*_bitmap path (expected ~10s) per the cycle-001 contract.

logger::log_info("FR-007a — running 6-facet overview")
fr007a_data <- fr007a_overview()

# Append the queries to the sidecar's queries list so every value
# traces back. Use a single placeholder record at the moment — each
# per-facet query is logged via pg_query_panel and its hash is
# captured there; we don't currently combine them into one record.
queries <- c(queries, list(
    list(sql = "fr007a_overview()", row_count = nrow(fr007a_data),
         result_hash = substr(
             digest::digest(fr007a_data, algo = "sha256"), 1L, 12L
         ))
))

fr007a_plot <- plot_fr007a_overview(fr007a_data, width_mm = 320L)
# Full 6-facet × per-band resource matrix — moves to supplementary
# in the iteration plan; still emitted as a stand-alone artifact.
ggsave(file.path(out_dir, "fr007a-overview-supplementary.pdf"), fr007a_plot,
       width = 320, height = 200, units = "mm")
ggsave(file.path(out_dir, "fr007a-overview-supplementary.svg"), fr007a_plot,
       width = 320, height = 200, units = "mm")
logger::log_info(
    "FR-007a overview (supplementary) written to ",
    "{out_dir}/fr007a-overview-supplementary.{{pdf,svg}}"
)

# Tiny "Total-only" variant — 6 facets stacked vertically, one pair
# of horizontal bars each. This is the candidate for the main
# composite Figure 1.
fr007a_total_plot <- plot_fr007a_total(fr007a_data, width_mm = 180L)
ggsave(file.path(out_dir, "fr007a-total.pdf"), fr007a_total_plot,
       width = 180, height = 200, units = "mm")
ggsave(file.path(out_dir, "fr007a-total.svg"), fr007a_total_plot,
       width = 180, height = 200, units = "mm")
logger::log_info("FR-007a total written to {out_dir}/fr007a-total.{{pdf,svg}}")

# ---- FR-007e — structural specificity × chemical category (dev5) ----------
#
# Bitmap intersection of structural_specificity x chemical_class/
# metabolic_domain on dev5 (~25 ms data layer). Stand-alone artifact.

logger::log_info("FR-007e — running specificity x category")
fr007e_data <- fr007e_specificity_by_category()
queries <- c(queries, list(
    list(sql = "fr007e_specificity_by_category()",
         row_count = nrow(fr007e_data),
         result_hash = substr(
             digest::digest(fr007e_data, algo = "sha256"), 1L, 12L
         ))
))
fr007e_plot <- plot_fr007e(fr007e_data, width_mm = 180L)
ggsave(file.path(out_dir, "fr007e-specificity.pdf"), fr007e_plot,
       width = 180, height = 110, units = "mm")
ggsave(file.path(out_dir, "fr007e-specificity.svg"), fr007e_plot,
       width = 180, height = 110, units = "mm")
logger::log_info("FR-007e specificity written to {out_dir}/fr007e-specificity.{{pdf,svg}}")

# ---- FR-007b — coverage profile (Entities / Molecular / Structures) -------

logger::log_info("FR-007b — running coverage profile (all variants)")
fr007b_data <- fr007b_coverage("all")
queries <- c(queries, list(
    list(sql = "fr007b_coverage(\"all\")",
         row_count = nrow(fr007b_data),
         result_hash = substr(
             digest::digest(fr007b_data, algo = "sha256"), 1L, 12L
         ))
))
fr007b_plot <- plot_fr007b_coverage(fr007b_data, width_mm = 180L)
ggsave(file.path(out_dir, "fr007b-coverage.pdf"), fr007b_plot,
       width = 180, height = 100, units = "mm")
ggsave(file.path(out_dir, "fr007b-coverage.svg"), fr007b_plot,
       width = 180, height = 100, units = "mm")
logger::log_info("FR-007b coverage written to {out_dir}/fr007b-coverage.{{pdf,svg}}")

# ---- FR-007d — entity x interaction-type matrix (top participant types) ---

logger::log_info("FR-007d — running entity x interaction-type matrix")
fr007d_data <- fr007d_entity_x_interaction(n_types = 8L)
queries <- c(queries, list(
    list(sql = "fr007d_entity_x_interaction()",
         row_count = nrow(fr007d_data),
         result_hash = substr(
             digest::digest(fr007d_data, algo = "sha256"), 1L, 12L
         ))
))
fr007d_plot <- plot_fr007d_matrix(fr007d_data, width_mm = 180L)
ggsave(file.path(out_dir, "fr007d-matrix.pdf"), fr007d_plot,
       width = 180, height = 180, units = "mm")
ggsave(file.path(out_dir, "fr007d-matrix.svg"), fr007d_plot,
       width = 180, height = 180, units = "mm")
logger::log_info("FR-007d matrix written to {out_dir}/fr007d-matrix.{{pdf,svg}}")

# ---- FR-007c — resource-overlap networks (Molecular entities, Interactions)

logger::log_info("FR-007c — running resource overlap networks")
fr007c_data <- fr007c_overlap()
queries <- c(queries, list(
    list(sql = "fr007c_overlap()",
         row_count = nrow(fr007c_data),
         result_hash = substr(
             digest::digest(fr007c_data, algo = "sha256"), 1L, 12L
         ))
))
fr007c_plot <- plot_fr007c_networks(fr007c_data,
                                    min_overlap = 100L,
                                    width_mm    = 320L)
ggsave(file.path(out_dir, "fr007c-networks.pdf"), fr007c_plot,
       width = 320, height = 160, units = "mm")
ggsave(file.path(out_dir, "fr007c-networks.svg"), fr007c_plot,
       width = 320, height = 160, units = "mm")
logger::log_info("FR-007c networks written to {out_dir}/fr007c-networks.{{pdf,svg}}")

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

# ---- Caption (FR-040..FR-041a, SC-011) -------------------------------------
#
# Panel count is hard-coded to 6 (A=architecture + B/C/D/E/F
# quantitative panels) until composition.yaml lands; the rebuild
# fails if caption.tex declares a different number of (a)/(b)/...
# sub-letters (FR-041b).

caption_info <- compose_caption(
    figure_id      = "fig01-overview",
    composite_pdf  = file.path(out_dir, "fig01-overview.pdf"),
    caption_source = "figures/fig01-overview/caption.tex",
    out_dir        = out_dir,
    panel_count    = 7L
)

# ---- Panel A statistics digest (FR-043 family) -----------------------------
#
# Side-output of the Figure 1 build: queries dev5 for the five
# sections of headline counts that back the numbers the author
# hand-transcribes into the Inkscape-authored Panel A diagram.
# Validates against contracts/panel-a-stats.schema.json (FR-043g)
# and asserts the FR-043h snapshot binding (digest's snapshot id
# MUST equal the composite's).

logger::log_info("Building Panel A statistics digest (FR-043)")
digest_result <- build_panel_a_digest(
    snapshot_id = snapshot_id(dep5$manifest),
    out_dir     = "figures/fig01-overview/panel-a-stats",
    caption_sty = "tex/caption.sty"
)

# ---- Provenance sidecar ----------------------------------------------------

write_sidecar(
    artifact_id    = "fig01-overview",
    artifact_path  = file.path(out_dir, "fig01-overview.pdf"),
    deployments    = list(dep5$deployment),
    manifests      = list(dep5$manifest),
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
        ),
        list(
            kind        = "panel-a-stats-digest",
            path        = digest_result$sidecar,
            source      = paste0(
                "Pipeline-generated FR-043 digest; ",
                "config sha256=", digest_result$config_sha256
            ),
            fingerprint = digest_result$snapshot_id
        )
    ),
    parameters     = list(width_mm = 180L),
    seed           = pipeline_seed(),
    caption        = list(
        source_path           = caption_info$caption_source,
        source_kind           = if (endsWith(caption_info$caption_source, ".md")) "md" else "tex",
        panel_letter_count    = caption_info$panel_letter_count,
        composite_panel_count = caption_info$composite_panel_count
    )
)

pg_close_all_panel()
logger::log_info("fig01-overview complete")
