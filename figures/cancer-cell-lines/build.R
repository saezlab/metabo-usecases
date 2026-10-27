# figures/cancer-cell-lines/build.R
#
# Cancer cell lines use case (Shorthouse 2022 + azelate webapp queries).
#
# Panels:
#   A — KRAS volcano (limma logFC vs. -log10 P)
#   B — EGFR volcano (same axis limits as A, FR-019)
#   C — Azelate interaction-type composition by upstream source
#       (from data/derived/azelate/csvs/interaction_types_by_source.csv)
#   D — Azelate cancer associations by sample type
#       (from data/derived/azelate/csvs/cancer_assoc_by_sample_type.csv)
#
# The former COSMOS PKN panels (E, F) are dropped; their analysis is
# kept as an archive in analyses/cancer-cell-lines-cosmos/.
#
# Inputs consumed (all repo-local):
#   - data/derived/cancer-cell-lines/ExtendedDataTable_...xlsx
#   - data/derived/azelate/csvs/{interaction_types_by_source,cancer_assoc_by_sample_type}.csv
#
# The refactor of the analysis Rmds into the pipeline is deferred
# (tasks T058–T062).

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

# ── Composite-compact label scales (memory: composite_panel_label_sizes) ────
#
# font_scale: body / axis / title size multiplier. ~1.7 lands axis
# labels around 10 pt at the panel's rendered size in the composite.
# legend_scale: separately tuned so many-entry legends don't crowd the
# plot area. 1.1 lands legend text around 6.6 pt with a proportionally
# smaller key swatch.
font_scale   <- 1.7
legend_scale <- 1.1

# ── Data loading ─────────────────────────────────────────────────────────────

logger::log_info("[cancer-cell-lines] loading differential-analysis xlsx (KRAS, EGFR)")
kras_dem <- case_study_differential("KRAS")
egfr_dem <- case_study_differential("EGFR")

logger::log_info("[cancer-cell-lines] loading azelate Panel C / Panel D source CSVs")
azelate_c <- azelate_panel_c_data()
azelate_d <- azelate_panel_d_data()

# ── Category-colour registration ─────────────────────────────────────────────

register_azelate_colours()

# ── Volcano shared limits (FR-019) ───────────────────────────────────────────

volcano_lims <- volcano_shared_limits(list(kras_dem, egfr_dem))

# ── Panel rendering ──────────────────────────────────────────────────────────

# Highlight: azelaic acid (HMDB0000784) — Panels A / B + downstream
# Panel C / D draw on the same metabolite as the worked example.
azelate_highlight <- c(HMDB0000784 = "Azelaic acid")

logger::log_info("[cancer-cell-lines] rendering Panel A (KRAS volcano)")
panel_a <- volcano_panel(
    kras_dem,
    contrast_label = "KRAS",
    xlim           = volcano_lims$xlim,
    ylim           = volcano_lims$ylim,
    width_mm       = 89L,
    font_scale     = font_scale,
    legend_scale   = legend_scale,
    highlight      = azelate_highlight
)

logger::log_info("[cancer-cell-lines] rendering Panel B (EGFR volcano)")
panel_b <- volcano_panel(
    egfr_dem,
    contrast_label = "EGFR",
    xlim           = volcano_lims$xlim,
    ylim           = volcano_lims$ylim,
    width_mm       = 89L,
    font_scale     = font_scale,
    legend_scale   = legend_scale,
    highlight      = azelate_highlight
)

logger::log_info("[cancer-cell-lines] rendering Panel C (azelate interaction types)")
panel_c <- azelate_interaction_panel(
    azelate_c,
    width_mm     = 89L,
    font_scale   = font_scale,
    legend_scale = legend_scale
)

logger::log_info("[cancer-cell-lines] rendering Panel D (azelate cancer associations)")
panel_d <- azelate_disease_panel(
    azelate_d,
    width_mm     = 89L,
    font_scale   = font_scale,
    legend_scale = legend_scale
)

# ── Save individual panels (per FR-026: SVG + PDF per panel) ─────────────────

individual_panels <- list(
    panel_a            = panel_a,
    panel_b            = panel_b,
    panel_c            = panel_c,
    panel_d            = panel_d
)

for (name in names(individual_panels)) {
    ggsave(
        filename = file.path(out_dir, paste0(name, ".pdf")),
        plot     = individual_panels[[name]],
        width    = 89,
        height   = 70,
        units    = "mm"
    )
    ggsave(
        filename = file.path(out_dir, paste0(name, ".svg")),
        plot     = individual_panels[[name]],
        width    = 89,
        height   = 70,
        units    = "mm"
    )
}

# ── Composite (nested patchwork) ─────────────────────────────────────────────
#
# Conceptual layout:
#   Row 1: A | B                          (volcanos)
#   Row 2: C | D                          (azelate)

row_ab <- patchwork::wrap_plots(
    list(
        panel_a + ggplot2::labs(tag = "A"),
        panel_b + ggplot2::labs(tag = "B")
    ),
    ncol = 2L
)

row_cd <- patchwork::wrap_plots(
    list(
        panel_c + ggplot2::labs(tag = "C"),
        panel_d + ggplot2::labs(tag = "D")
    ),
    ncol = 2L
)

composite <- row_ab / row_cd

logger::log_info(paste0(
    "[cancer-cell-lines] assembled composite (2 rows; ",
    "font_scale={font_scale}, legend_scale={legend_scale})"
))

# 180 mm × 110 mm — two rows at the row height of the former
# three-row layout.
ggsave(
    filename = file.path(out_dir, "cancer-cell-lines.pdf"),
    plot     = composite,
    width    = 180,
    height   = 110,
    units    = "mm"
)
ggsave(
    filename = file.path(out_dir, "cancer-cell-lines.svg"),
    plot     = composite,
    width    = 180,
    height   = 110,
    units    = "mm"
)

# ── Caption ──────────────────────────────────────────────────────────────────

caption_info <- compose_caption(
    figure_id      = "cancer-cell-lines",
    composite_pdf  = file.path(
        out_dir, "cancer-cell-lines.pdf"
    ),
    caption_source = "figures/cancer-cell-lines/caption.tex",
    out_dir        = out_dir,
    panel_count    = 4L
)

# ── Provenance sidecar ────────────────────────────────────────────────────────

external_inputs <- list(
    list(
        kind        = "shorthouse-differential-analysis",
        path        = attr(kras_dem, "source_path"),
        fingerprint = attr(kras_dem, "fingerprint"),
        source      = "analyses/cancer-cell-lines/04_differential_analysis.Rmd",
        sheets      = c("KRAS_filt_limma", "EGFR_filt_limma")
    ),
    list(
        kind        = "azelate-panel-c",
        path        = attr(azelate_c, "source_path"),
        fingerprint = attr(azelate_c, "fingerprint"),
        source      = "analyses/azelate/azelaic_acid_query_tables.py"
    ),
    list(
        kind        = "azelate-panel-d",
        path        = attr(azelate_d, "source_path"),
        fingerprint = attr(azelate_d, "fingerprint"),
        source      = "analyses/azelate/azelaic_acid_query_tables.py"
    )
)

write_sidecar(
    artifact_id     = "cancer-cell-lines",
    artifact_path   = file.path(
        out_dir, "cancer-cell-lines.pdf"
    ),
    deployments     = list(dep5$deployment),
    manifests       = list(dep5$manifest),
    script_path     = "figures/cancer-cell-lines/build.R",
    queries         = list(),
    external_inputs = external_inputs,
    parameters      = list(
        volcano_pvalue_threshold = 0.05,
        volcano_logfc_threshold  = 0.5,
        font_scale               = font_scale,
        legend_scale             = legend_scale,
        composite_dims_mm        = list(width = 180L, height = 110L),
        azelate_source           = paste0(
            "CSVs from data/derived/azelate/csvs/ ",
            "(generated by azelaic_acid_query_tables.py against ",
            "dev.omnipathdb.org)"
        )
    ),
    seed    = pipeline_seed(),
    caption = caption_info
)

logger::log_info("[cancer-cell-lines] build complete — outputs in {out_dir}")
