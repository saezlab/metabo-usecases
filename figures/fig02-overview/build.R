# figures/fig02-overview/build.R
#
# Orchestrates the Figure 2 build (post-2026-06-14 six-figure
# renumbering): queries dev5 for the FR-007 quantitative panel data,
# renders every FR-007a..f variant + the composite-selected subset
# (A: fr007a-total, B: fr007e, C: fr007c, D: fr007d) via ggplot,
# composes via tex/compose_fig02.tex, and writes the provenance
# sidecar. The architecture asset (formerly Panel A) and the FR-043
# statistics digest moved to Figure 1 — see
# figures/fig01-architecture/build.R.
#
# Sourced by rebuild.R; safe to source standalone too.

suppressPackageStartupMessages({
    library(metabo.figures)
    library(ggplot2)
})

setup_pipeline_log("build:fig02-overview")
set.seed(pipeline_seed())

out_dir <- "figures/fig02-overview/out"
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

# ---- Composite -------------------------------------------------------------
#
# Post-2026-06-14 six-figure renumbering: the architecture asset
# moved to Figure 1 (figures/fig01-architecture/). Figure 2's
# composite assembles the four FR-007 panels via patchwork —
# A: fr007a-total spans the full 180 mm top row (~80 mm tall, six
# horizontal facets), B/C/D fill the bottom row at ~60 mm wide each
# (~80 mm tall). Half-page composite ~180 × 160 mm.

logger::log_info("Composing Figure 2 (patchwork: A wide / B,C,D row)")

# Add capital-letter tag annotations so the composite carries the
# A/B/C/D labels FR-024 mandates.
panel_a_tagged <- fr007a_total_plot +
    patchwork::plot_annotation(tag_levels = list(c("A")))
bottom_row <- (
    (fr007e_plot +
        patchwork::plot_annotation(tag_levels = list(c("B")))) |
    (fr007c_plot +
        patchwork::plot_annotation(tag_levels = list(c("C")))) |
    (fr007d_plot +
        patchwork::plot_annotation(tag_levels = list(c("D"))))
)

composite <- panel_a_tagged / bottom_row +
    patchwork::plot_layout(heights = c(1, 1))

# 180 × 160 mm half-page composite per Session 2026-06-14.
ggsave(
    file.path(out_dir, "fig02-overview.pdf"), composite,
    width = 180, height = 160, units = "mm"
)
ggsave(
    file.path(out_dir, "fig02-overview.svg"), composite,
    width = 180, height = 160, units = "mm"
)
logger::log_info(
    "Figure 2 composite written to ",
    "{out_dir}/fig02-overview.{{pdf,svg}}"
)

# ---- Caption (FR-040..FR-041a, SC-011) -------------------------------------
#
# Figure 2 is a four-panel composite (A: fr007a-total, B: fr007e,
# C: fr007c, D: fr007d). FR-041b requires the caption's (a)/(b)/(c)/(d)
# sub-letter count to equal the composite panel count (4).

caption_info <- compose_caption(
    figure_id      = "fig02-overview",
    composite_pdf  = file.path(out_dir, "fig02-overview.pdf"),
    caption_source = "figures/fig02-overview/caption.tex",
    out_dir        = out_dir,
    panel_count    = 4L
)

# ---- Provenance sidecar ----------------------------------------------------

write_sidecar(
    artifact_id    = "fig02-overview",
    artifact_path  = file.path(out_dir, "fig02-overview.pdf"),
    deployments    = list(dep5$deployment),
    manifests      = list(dep5$manifest),
    script_path    = "figures/fig02-overview/build.R",
    queries        = queries,
    external_inputs = list(),
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
logger::log_info("fig02-overview complete")
