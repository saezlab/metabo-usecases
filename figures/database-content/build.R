# figures/database-content/build.R
#
# Orchestrates the Figure 2 build (post-2026-06-14 six-figure
# renumbering): queries dev5 for the FR-007 quantitative panel data,
# renders every FR-007a..f variant + the composite-selected subset
# (A: resource-overview-total, B: specificity, C: resource-overlap, D: entity-interaction) via ggplot,
# composes via tex/compose_fig02.tex, and writes the provenance
# sidecar. The architecture asset (formerly Panel A) and the FR-043
# statistics digest moved to Figure 1 — see
# figures/architecture/build.R.
#
# Sourced by rebuild.R; safe to source standalone too.

suppressPackageStartupMessages({
    library(metabo.figures)
    library(ggplot2)
})

setup_pipeline_log("build:database-content")
set.seed(pipeline_seed())

out_dir <- "figures/database-content/out"
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
resource_overview_data <- resource_overview_overview()

# Append the queries to the sidecar's queries list so every value
# traces back. Use a single placeholder record at the moment — each
# per-facet query is logged via pg_query_panel and its hash is
# captured there; we don't currently combine them into one record.
queries <- c(queries, list(
    list(sql = "resource_overview_overview()", row_count = nrow(resource_overview_data),
         result_hash = substr(
             digest::digest(resource_overview_data, algo = "sha256"), 1L, 12L
         ))
))

resource_overview_plot <- plot_resource_overview(resource_overview_data, width_mm = 320L)
# Full 6-facet × per-band resource matrix — moves to supplementary
# in the iteration plan; still emitted as a stand-alone artifact.
ggsave(file.path(out_dir, "resource-overview-supplementary.pdf"), resource_overview_plot,
       width = 320, height = 200, units = "mm")
ggsave(file.path(out_dir, "resource-overview-supplementary.svg"), resource_overview_plot,
       width = 320, height = 200, units = "mm")
logger::log_info(
    "FR-007a overview (supplementary) written to ",
    "{out_dir}/resource-overview-supplementary.{{pdf,svg}}"
)

# Tiny "Total-only" variant — 6 facets stacked vertically, one pair
# of horizontal bars each. This is the candidate for the main
# composite Figure 1.
resource_overview_total_plot <- plot_resource_overview_total(resource_overview_data, width_mm = 180L)
ggsave(file.path(out_dir, "resource-overview-total.pdf"), resource_overview_total_plot,
       width = 180, height = 200, units = "mm")
ggsave(file.path(out_dir, "resource-overview-total.svg"), resource_overview_total_plot,
       width = 180, height = 200, units = "mm")
logger::log_info("FR-007a total written to {out_dir}/resource-overview-total.{{pdf,svg}}")

# ---- FR-007e — structural specificity × chemical category (dev5) ----------
#
# Bitmap intersection of structural_specificity x chemical_class/
# metabolic_domain on dev5 (~25 ms data layer). Stand-alone artifact.

logger::log_info("FR-007e — running specificity x category")
specificity_data <- specificity_by_category()
queries <- c(queries, list(
    list(sql = "specificity_by_category()",
         row_count = nrow(specificity_data),
         result_hash = substr(
             digest::digest(specificity_data, algo = "sha256"), 1L, 12L
         ))
))
specificity_plot <- plot_specificity_by_category(specificity_data, width_mm = 180L)
ggsave(file.path(out_dir, "specificity-by-category.pdf"), specificity_plot,
       width = 180, height = 110, units = "mm")
ggsave(file.path(out_dir, "specificity-by-category.svg"), specificity_plot,
       width = 180, height = 110, units = "mm")
logger::log_info("FR-007e specificity written to {out_dir}/specificity-by-category.{{pdf,svg}}")

# ---- FR-007b — coverage profile (Entities / Molecular / Structures) -------

logger::log_info("FR-007b — running coverage profile (all variants)")
coverage_profile_data <- coverage_profile("all")
queries <- c(queries, list(
    list(sql = "coverage_profile(\"all\")",
         row_count = nrow(coverage_profile_data),
         result_hash = substr(
             digest::digest(coverage_profile_data, algo = "sha256"), 1L, 12L
         ))
))
coverage_profile_plot <- plot_coverage_profile(coverage_profile_data, width_mm = 180L)
ggsave(file.path(out_dir, "coverage-profile.pdf"), coverage_profile_plot,
       width = 180, height = 100, units = "mm")
ggsave(file.path(out_dir, "coverage-profile.svg"), coverage_profile_plot,
       width = 180, height = 100, units = "mm")
logger::log_info("FR-007b coverage written to {out_dir}/coverage-profile.{{pdf,svg}}")

# ---- FR-007d — entity x interaction-type matrix (top participant types) ---

logger::log_info("FR-007d — running entity x interaction-type matrix")
entity_interaction_data <- entity_by_interaction_type(n_types = 8L)
queries <- c(queries, list(
    list(sql = "entity_by_interaction_type()",
         row_count = nrow(entity_interaction_data),
         result_hash = substr(
             digest::digest(entity_interaction_data, algo = "sha256"), 1L, 12L
         ))
))
entity_interaction_plot <- plot_entity_interaction_matrix(entity_interaction_data, width_mm = 180L)
ggsave(file.path(out_dir, "entity-interaction-matrix.pdf"), entity_interaction_plot,
       width = 180, height = 180, units = "mm")
ggsave(file.path(out_dir, "entity-interaction-matrix.svg"), entity_interaction_plot,
       width = 180, height = 180, units = "mm")
logger::log_info("FR-007d matrix written to {out_dir}/entity-interaction-matrix.{{pdf,svg}}")

# ---- FR-007c — resource-overlap networks (Molecular entities, Interactions)

logger::log_info("FR-007c — running resource overlap networks")
resource_overlap_data <- resource_overlap()
queries <- c(queries, list(
    list(sql = "resource_overlap()",
         row_count = nrow(resource_overlap_data),
         result_hash = substr(
             digest::digest(resource_overlap_data, algo = "sha256"), 1L, 12L
         ))
))
resource_overlap_plot <- plot_resource_overlap(resource_overlap_data,
                                    min_overlap = 100L,
                                    width_mm    = 320L)
ggsave(file.path(out_dir, "resource-overlap-networks.pdf"), resource_overlap_plot,
       width = 320, height = 160, units = "mm")
ggsave(file.path(out_dir, "resource-overlap-networks.svg"), resource_overlap_plot,
       width = 320, height = 160, units = "mm")
logger::log_info("FR-007c networks written to {out_dir}/resource-overlap-networks.{{pdf,svg}}")

# ---- Composite -------------------------------------------------------------
#
# Post-2026-06-14 six-figure renumbering: the architecture asset
# moved to Figure 1 (figures/architecture/). Figure 2's
# composite assembles the four FR-007 panels via patchwork —
# A: resource-overview-total spans the full 180 mm top row (~80 mm tall, six
# horizontal facets), B/C/D fill the bottom row at ~60 mm wide each
# (~80 mm tall). Half-page composite ~180 × 160 mm.

logger::log_info("Composing Figure 2 (patchwork: A wide / B,C,D row)")

# Use patchwork wrap_elements() to make the inner resource-overview-total
# patchwork an atomic unit so its 6 sub-facets don't get
# auto-tagged. Outer-level plot_annotation(tag_levels = "A") then
# assigns A / B / C / D to the four top-level slots.
# Wrap every nested patchwork (resource-overview-total is a 6-facet row;
# resource-overlap is a 2-network row) in wrap_elements() so the outer
# composite treats them as atomic cells. Without this, the inner
# plot_layout / plot_annotation calls bleed up into the outer
# grid.
#
# Layout (180 × 180 mm, portrait):
#   row 1 — A: resource-overview-total           full width × ~45 mm
#   row 2 — B: resource-overlap networks        full width × ~80 mm
#   row 3 — C: specificity | D: entity-interaction     1/3 + 2/3 width × ~55 mm
#
# The top row is flat by design — each facet has just 2 bars,
# so a thicker top row would make those bars top-heavy. Row 3's
# 1:2 width split gives entity-interaction (8 entity-type facets, 2x4 grid)
# the room it needs while specificity (6 chemical-category facets,
# 2x3 grid) compresses comfortably.
bottom_row <- (
    patchwork::wrap_elements(full = specificity_plot) |
    patchwork::wrap_elements(full = entity_interaction_plot)
) +
    # C : D = 2 : 3 (≈72 mm : 108 mm at 180 mm composite width).
    # Slightly wider C than the previous 1 : 2 split so the
    # specificity facet titles + Y-axis title breathe; D drops
    # from 120 mm → 108 mm but stays large enough for the 2 x 4
    # entity-type facet grid.
    patchwork::plot_layout(widths = c(2, 3))

composite <- (
    patchwork::wrap_elements(full = resource_overview_total_plot)
    /
    patchwork::wrap_elements(full = resource_overlap_plot)
    /
    bottom_row
) +
    patchwork::plot_layout(heights = c(45, 80, 55)) +
    patchwork::plot_annotation(
        tag_levels = "A",
        theme = ggplot2::theme(
            # Larger + bolder panel letters per FR-024 (capital
            # labels) and Session 2026-06-14 feedback.
            plot.tag = ggplot2::element_text(
                size = 18, face = "bold"
            ),
            plot.tag.position = c(0.01, 0.99)
        )
    )

# 180 × 180 mm composite (rows: A 45 mm, B 80 mm, C/D 55 mm).
ggsave(
    file.path(out_dir, "database-content.pdf"), composite,
    width = 180, height = 180, units = "mm"
)
try(
    ggsave(
        file.path(out_dir, "database-content.svg"), composite,
        width = 180, height = 180, units = "mm"
    ),
    silent = FALSE
)
logger::log_info(
    "Figure 2 composite written to ",
    "{out_dir}/database-content.{{pdf,svg}}"
)

# ---- Caption (FR-040..FR-041a, SC-011) -------------------------------------
#
# Figure 2 is a four-panel composite (A: resource-overview-total, B: specificity,
# C: resource-overlap, D: entity-interaction). FR-041b requires the caption's (a)/(b)/(c)/(d)
# sub-letter count to equal the composite panel count (4).

caption_info <- compose_caption(
    figure_id      = "database-content",
    composite_pdf  = file.path(out_dir, "database-content.pdf"),
    caption_source = "figures/database-content/caption.tex",
    out_dir        = out_dir,
    panel_count    = 4L
)

# ---- Provenance sidecar ----------------------------------------------------

write_sidecar(
    artifact_id    = "database-content",
    artifact_path  = file.path(out_dir, "database-content.pdf"),
    deployments    = list(dep5$deployment),
    manifests      = list(dep5$manifest),
    script_path    = "figures/database-content/build.R",
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
logger::log_info("database-content complete")
