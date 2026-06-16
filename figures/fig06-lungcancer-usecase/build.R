# figures/fig06-lungcancer-usecase/build.R
#
# Figure 6: lung-cancer use case (Shorthouse 2022 + COSMOS PKN +
# azelate webapp queries).
#
# Panels:
#   A — KRAS volcano (limma logFC vs. -log10 P)
#   B — EGFR volcano (same axis limits as A, FR-019)
#   C — Azelate interaction-type composition by upstream source
#       (from azelate_fig_tables/outputs/csvs/fig3_panel_C_data.csv)
#   D — Azelate cancer associations by sample type
#       (from azelate_fig_tables/outputs/csvs/fig3_panel_D_data.csv)
#   E — COSMOS PKN GEM/allosteric edge counts by direction × resource
#       (one sub-panel per contrast: KRAS, EGFR)
#   F — Subcellular-compartment edges (KRAS up/down, EGFR up/down)
#
# Inputs consumed (all repo-local):
#   - omnipath_metabo_case1/Results/Differential_Analysis/...xlsx
#   - omnipath_metabo_case1/data/pkn_{allosteric,enzyme_metabolite}.csv
#   - omnipath_metabo_case1/azelate_fig_tables/outputs/csvs/fig3_panel_{C,D}_data.csv
#
# The full Differential_Analysis.Rmd / 01_cosmos_pkn.py refactor is
# deferred to tasks T058–T062 and T066 (the latter will move from
# the vendored pkn_*.csv fixtures to the omnipath-client API).

suppressPackageStartupMessages({
    library(metabo.figures)
    library(ggplot2)
    library(patchwork)
})

setup_pipeline_log("build:fig06-lungcancer-usecase")
set.seed(pipeline_seed())

out_dir <- "figures/fig06-lungcancer-usecase/out"
fs::dir_create(out_dir)

# ── Deployment provenance (anchor the snapshot identity) ─────────────────────

dep5 <- deployment_provenance("dev5")

# ── Composite-compact label scales (memory: composite_panel_label_sizes) ────
#
# font_scale: body / axis / title size multiplier. ~1.7 lands axis
# labels around 10 pt at the panel's rendered size in the composite.
# legend_scale: separately tuned so the many-entry legends in Panels
# E and F (9 resources) don't crowd the plot area. 1.1 lands legend
# text around 6.6 pt with a proportionally smaller key swatch.
font_scale   <- 1.7
legend_scale <- 1.1

# ── Data loading ─────────────────────────────────────────────────────────────

logger::log_info("[fig06] loading differential-analysis xlsx (KRAS, EGFR)")
kras_dem <- case_study_differential("KRAS")
egfr_dem <- case_study_differential("EGFR")

logger::log_info("[fig06] loading COSMOS PKN fixtures")
pkn <- case_study_cosmos_pkn()

logger::log_info("[fig06] loading azelate Panel C / Panel D source CSVs")
azelate_c <- azelate_panel_c_data()
azelate_d <- azelate_panel_d_data()

# ── Category-colour registration ─────────────────────────────────────────────

register_case_study_resource_colours()
register_azelate_colours()

# ── Volcano shared limits (FR-019) ───────────────────────────────────────────

volcano_lims <- volcano_shared_limits(list(kras_dem, egfr_dem))

# ── Panel rendering ──────────────────────────────────────────────────────────

logger::log_info("[fig06] rendering Panel A (KRAS volcano)")
panel_a <- volcano_panel(
    kras_dem,
    contrast_label = "KRAS",
    xlim           = volcano_lims$xlim,
    ylim           = volcano_lims$ylim,
    width_mm       = 89L,
    font_scale     = font_scale,
    legend_scale   = legend_scale
)

logger::log_info("[fig06] rendering Panel B (EGFR volcano)")
panel_b <- volcano_panel(
    egfr_dem,
    contrast_label = "EGFR",
    xlim           = volcano_lims$xlim,
    ylim           = volcano_lims$ylim,
    width_mm       = 89L,
    font_scale     = font_scale,
    legend_scale   = legend_scale
)

logger::log_info("[fig06] rendering Panel C (azelate interaction types)")
panel_c <- azelate_interaction_panel(
    azelate_c,
    width_mm     = 89L,
    font_scale   = font_scale,
    legend_scale = legend_scale
)

logger::log_info("[fig06] rendering Panel D (azelate cancer associations)")
panel_d <- azelate_disease_panel(
    azelate_d,
    width_mm     = 89L,
    font_scale   = font_scale,
    legend_scale = legend_scale
)

logger::log_info("[fig06] building Panel E PKN summaries (KRAS, EGFR)")
kras_summary <- case_study_pkn_summary(kras_dem, pkn)
egfr_summary <- case_study_pkn_summary(egfr_dem, pkn)

panel_e_kras <- gem_allosteric_panel(
    kras_summary,
    contrast_label = "KRAS",
    width_mm       = 89L,
    font_scale     = font_scale,
    legend_scale   = legend_scale
)
panel_e_egfr <- gem_allosteric_panel(
    egfr_summary,
    contrast_label = "EGFR",
    width_mm       = 89L,
    font_scale     = font_scale,
    legend_scale   = legend_scale
)

logger::log_info("[fig06] rendering Panel F sub-panels (location × direction)")
panel_f_kras_up <- subcellular_location_panel(
    kras_summary, "KRAS", "up",
    width_mm = 89L, font_scale = font_scale,
    legend_scale = legend_scale
)
panel_f_kras_down <- subcellular_location_panel(
    kras_summary, "KRAS", "down",
    width_mm = 89L, font_scale = font_scale,
    legend_scale = legend_scale
)
panel_f_egfr_up <- subcellular_location_panel(
    egfr_summary, "EGFR", "up",
    width_mm = 89L, font_scale = font_scale,
    legend_scale = legend_scale
)
panel_f_egfr_down <- subcellular_location_panel(
    egfr_summary, "EGFR", "down",
    width_mm = 89L, font_scale = font_scale,
    legend_scale = legend_scale
)

# ── Save individual panels (per FR-026: SVG + PDF per panel) ─────────────────

individual_panels <- list(
    panel_a            = panel_a,
    panel_b            = panel_b,
    panel_c            = panel_c,
    panel_d            = panel_d,
    panel_e_kras       = panel_e_kras,
    panel_e_egfr       = panel_e_egfr,
    panel_f_kras_up    = panel_f_kras_up,
    panel_f_kras_down  = panel_f_kras_down,
    panel_f_egfr_up    = panel_f_egfr_up,
    panel_f_egfr_down  = panel_f_egfr_down
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
#   Row 1: A | B                                       (volcanos)
#   Row 2: C | D                                       (azelate)
#   Row 3: E_kras | E_egfr                             (GEM/Allos)
#   Row 4: F_kras_up | F_kras_down | F_egfr_up | F_egfr_down
#                                                      (subcellular)
#
# Panels E and F all use the same "resources" colour scale, so the
# (row_e / row_f) sub-patch collects guides into a single Resource
# legend at the bottom of the EF block. Heights weight row F a touch
# heavier because its 4-up x-axis labels need more vertical room.

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

row_e <- patchwork::wrap_plots(
    list(
        panel_e_kras + ggplot2::labs(tag = "E"),
        panel_e_egfr + ggplot2::labs(tag = "")
    ),
    ncol = 2L
)

row_f <- patchwork::wrap_plots(
    list(
        panel_f_kras_up   + ggplot2::labs(tag = "F"),
        panel_f_kras_down + ggplot2::labs(tag = ""),
        panel_f_egfr_up   + ggplot2::labs(tag = ""),
        panel_f_egfr_down + ggplot2::labs(tag = "")
    ),
    ncol = 4L
)

ef_block <- (row_e / row_f) +
    patchwork::plot_layout(
        guides  = "collect",
        heights = c(1, 1.05)
    ) &
    ggplot2::theme(legend.position = "bottom")

composite <- (row_ab / row_cd / ef_block) +
    patchwork::plot_layout(
        heights = c(1, 1, 2.25)
    )

logger::log_info(paste0(
    "[fig06] assembled compact composite (3 outer rows; EF shares ",
    "one legend; font_scale={font_scale}, legend_scale={legend_scale})"
))

# Compact: 180mm wide × 205mm tall — collecting the EF legend
# reclaims the per-panel legend gutters.
ggsave(
    filename = file.path(out_dir, "fig06-lungcancer-usecase.pdf"),
    plot     = composite,
    width    = 180,
    height   = 205,
    units    = "mm"
)
ggsave(
    filename = file.path(out_dir, "fig06-lungcancer-usecase.svg"),
    plot     = composite,
    width    = 180,
    height   = 205,
    units    = "mm"
)

# ── Caption ──────────────────────────────────────────────────────────────────

caption_info <- compose_caption(
    figure_id      = "fig06-lungcancer-usecase",
    composite_pdf  = file.path(
        out_dir, "fig06-lungcancer-usecase.pdf"
    ),
    caption_source = "figures/fig06-lungcancer-usecase/caption.tex",
    out_dir        = out_dir,
    panel_count    = 6L
)

# ── Provenance sidecar ────────────────────────────────────────────────────────

external_inputs <- list(
    list(
        kind        = "shorthouse-differential-analysis",
        path        = attr(kras_dem, "source_path"),
        fingerprint = attr(kras_dem, "fingerprint"),
        source      = "omnipath_metabo_case1 (legacy)",
        sheets      = c("KRAS_filt_limma", "EGFR_filt_limma")
    ),
    list(
        kind        = "cosmos-pkn-allosteric",
        path        = attr(pkn$allosteric, "source_path"),
        fingerprint = attr(pkn$allosteric, "fingerprint"),
        source      = "01_cosmos_pkn.py (legacy)"
    ),
    list(
        kind        = "cosmos-pkn-enzyme-metabolite",
        path        = attr(pkn$enzyme_metabolite, "source_path"),
        fingerprint = attr(pkn$enzyme_metabolite, "fingerprint"),
        source      = "01_cosmos_pkn.py (legacy)"
    ),
    list(
        kind        = "azelate-panel-c",
        path        = attr(azelate_c, "source_path"),
        fingerprint = attr(azelate_c, "fingerprint"),
        source      = "azelate_fig_tables/azelaic_acid_query_tables.py"
    ),
    list(
        kind        = "azelate-panel-d",
        path        = attr(azelate_d, "source_path"),
        fingerprint = attr(azelate_d, "fingerprint"),
        source      = "azelate_fig_tables/azelaic_acid_query_tables.py"
    )
)

write_sidecar(
    artifact_id     = "fig06-lungcancer-usecase",
    artifact_path   = file.path(
        out_dir, "fig06-lungcancer-usecase.pdf"
    ),
    deployments     = list(dep5$deployment),
    manifests       = list(dep5$manifest),
    script_path     = "figures/fig06-lungcancer-usecase/build.R",
    queries         = list(),
    external_inputs = external_inputs,
    parameters      = list(
        volcano_pvalue_threshold = 0.05,
        volcano_logfc_threshold  = 0.5,
        top_dems_per_direction   = 10L,
        font_scale               = font_scale,
        legend_scale             = legend_scale,
        composite_dims_mm        = list(width = 180L, height = 205L),
        cosmos_pkn_source        = paste0(
            "vendored fixture from omnipath_metabo_case1/data/ ",
            "(pending T066 omnipath-client refactor)"
        ),
        azelate_source           = paste0(
            "vendored CSVs from azelate_fig_tables/outputs/csvs/ ",
            "(generated by azelaic_acid_query_tables.py against ",
            "dev.omnipathdb.org)"
        )
    ),
    seed    = pipeline_seed(),
    caption = caption_info
)

logger::log_info("[fig06] build complete — outputs in {out_dir}")
