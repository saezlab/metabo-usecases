# figures/cancer-cell-lines/build.R
#
# Cancer cell lines use case: OmniPath Metabo evidence for azelate.
#
# Panels:
#   A — Azelate interaction-type composition by upstream source
#       (from data/derived/azelate/csvs/interaction_types_by_source.csv)
#   B — Azelate cancer associations by sample type
#       (from data/derived/azelate/csvs/cancer_assoc_by_sample_type.csv)
#
# The former KRAS / EGFR volcano panels (read from the previous Shorthouse
# pipeline) and the COSMOS PKN panels were removed in the 2026-11 cleanup;
# see RESTRUCTURING.md. The figure as of the first preprint is at the tag
# `preprint_v0`.
#
# Inputs consumed (all repo-local):
#   - data/derived/azelate/csvs/{interaction_types_by_source,cancer_assoc_by_sample_type}.csv

suppressPackageStartupMessages({
    library(metabo.figures)
    library(ggplot2)
    library(patchwork)
})

setup_pipeline_log("build:cancer-cell-lines")
set.seed(pipeline_seed())

out_dir <- "figures/cancer-cell-lines/out"
fs::dir_create(out_dir)

# ── Deployment provenance (anchor the snapshot identity) ─────────────────────

dep5 <- deployment_provenance("dev5")

# ── Composite-compact label scales ───────────────────────────────────────────
#
# font_scale: body / axis / title size multiplier. ~1.7 lands axis
# labels around 10 pt at the panel's rendered size in the composite.
# legend_scale: separately tuned so many-entry legends don't crowd the
# plot area. 1.1 lands legend text around 6.6 pt with a proportionally
# smaller key swatch.
font_scale   <- 1.7
legend_scale <- 1.1

# ── Data loading ─────────────────────────────────────────────────────────────

logger::log_info("[cancer-cell-lines] loading azelate source CSVs")
azelate_interactions <- azelate_panel_c_data()
azelate_diseases     <- azelate_panel_d_data()

# ── Category-colour registration ─────────────────────────────────────────────

register_azelate_colours()

# ── Panel rendering ──────────────────────────────────────────────────────────

logger::log_info("[cancer-cell-lines] rendering Panel A (azelate interaction types)")
panel_a <- azelate_interaction_panel(
    azelate_interactions,
    width_mm     = 89L,
    font_scale   = font_scale,
    legend_scale = legend_scale
)

logger::log_info("[cancer-cell-lines] rendering Panel B (azelate cancer associations)")
panel_b <- azelate_disease_panel(
    azelate_diseases,
    width_mm     = 89L,
    font_scale   = font_scale,
    legend_scale = legend_scale
)

# ── Save individual panels (SVG + PDF per panel) ─────────────────────────────

individual_panels <- list(
    panel_a = panel_a,
    panel_b = panel_b
)

for (name in names(individual_panels)) {
    for (ext in c("pdf", "svg")) {
        ggsave(
            filename = file.path(out_dir, paste0(name, ".", ext)),
            plot     = individual_panels[[name]],
            width    = 89,
            height   = 70,
            units    = "mm"
        )
    }
}

# ── Composite: A | B ─────────────────────────────────────────────────────────

composite <- patchwork::wrap_plots(
    list(
        panel_a + ggplot2::labs(tag = "A"),
        panel_b + ggplot2::labs(tag = "B")
    ),
    ncol = 2L
)

logger::log_info(paste0(
    "[cancer-cell-lines] assembled composite (1 row; ",
    "font_scale={font_scale}, legend_scale={legend_scale})"
))

composite_width_mm  <- 180L
composite_height_mm <- 55L

for (ext in c("pdf", "svg")) {
    ggsave(
        filename = file.path(out_dir, paste0("cancer-cell-lines.", ext)),
        plot     = composite,
        width    = composite_width_mm,
        height   = composite_height_mm,
        units    = "mm"
    )
}

# ── Caption ──────────────────────────────────────────────────────────────────

caption_info <- compose_caption(
    figure_id      = "cancer-cell-lines",
    composite_pdf  = file.path(out_dir, "cancer-cell-lines.pdf"),
    caption_source = "figures/cancer-cell-lines/caption.tex",
    out_dir        = out_dir,
    panel_count    = 2L
)

# ── Provenance sidecar ────────────────────────────────────────────────────────

external_inputs <- list(
    list(
        kind        = "azelate-interaction-types",
        path        = attr(azelate_interactions, "source_path"),
        fingerprint = attr(azelate_interactions, "fingerprint"),
        source      = "analyses/azelate/azelaic_acid_query_tables.py"
    ),
    list(
        kind        = "azelate-cancer-associations",
        path        = attr(azelate_diseases, "source_path"),
        fingerprint = attr(azelate_diseases, "fingerprint"),
        source      = "analyses/azelate/azelaic_acid_query_tables.py"
    )
)

write_sidecar(
    artifact_id     = "cancer-cell-lines",
    artifact_path   = file.path(out_dir, "cancer-cell-lines.pdf"),
    deployments     = list(dep5$deployment),
    manifests       = list(dep5$manifest),
    script_path     = "figures/cancer-cell-lines/build.R",
    queries         = list(),
    external_inputs = external_inputs,
    parameters      = list(
        font_scale        = font_scale,
        legend_scale      = legend_scale,
        composite_dims_mm = list(
            width  = composite_width_mm,
            height = composite_height_mm
        ),
        azelate_source    = paste0(
            "CSVs from data/derived/azelate/csvs/ ",
            "(generated by azelaic_acid_query_tables.py against ",
            "dev.omnipathdb.org)"
        )
    ),
    seed    = pipeline_seed(),
    caption = caption_info
)

logger::log_info("[cancer-cell-lines] build complete — outputs in {out_dir}")
