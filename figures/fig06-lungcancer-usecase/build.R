# figures/fig06-lungcancer-usecase/build.R
#
# Figure 6: lung-cancer use case (Shorthouse 2022 + COSMOS PKN).
#
# Panels:
#   A — KRAS volcano (limma logFC vs. -log10 P)
#   B — EGFR volcano (same axis limits as A, FR-019)
#   C — manual webapp screenshot (azelaic acid on metabo.omnipathdb.org)
#   D — MetaLinksDB 2.0 top-hit connections [STUB — pending DB rebuild]
#   E — GEM/allosteric edge counts by direction × resource (KRAS, EGFR)
#   F — subcellular-compartment edges (KRAS up/down, EGFR up/down)
#
# This first slice consumes pre-computed legacy artifacts:
#   - omnipath_metabo_case1/Results/Differential_Analysis/...xlsx
#   - omnipath_metabo_case1/data/pkn_{allosteric,enzyme_metabolite}.csv
# The full Differential_Analysis.Rmd / 01_cosmos_pkn.py refactor is
# tracked under tasks T058–T062 and T066 respectively.

suppressPackageStartupMessages({
    library(metabo.figures)
    library(ggplot2)
})

setup_pipeline_log("build:fig06-lungcancer-usecase")
set.seed(pipeline_seed())

out_dir <- "figures/fig06-lungcancer-usecase/out"
fs::dir_create(out_dir)

# ── Deployment provenance (anchor the snapshot identity) ─────────────────────

dep5 <- deployment_provenance("dev5")

# ── Data loading ─────────────────────────────────────────────────────────────

logger::log_info("[fig06] loading differential-analysis xlsx (KRAS, EGFR)")
kras_dem <- case_study_differential("KRAS")
egfr_dem <- case_study_differential("EGFR")

logger::log_info("[fig06] loading COSMOS PKN fixtures")
pkn <- case_study_cosmos_pkn()

# ── Category-colour registration ─────────────────────────────────────────────

register_case_study_resource_colours()

# ── Volcano shared limits (FR-019) ───────────────────────────────────────────

volcano_lims <- volcano_shared_limits(list(kras_dem, egfr_dem))

# ── Panel rendering ──────────────────────────────────────────────────────────

logger::log_info("[fig06] rendering Panel A (KRAS volcano)")
panel_a <- volcano_panel(
    kras_dem,
    contrast_label = "KRAS",
    xlim = volcano_lims$xlim,
    ylim = volcano_lims$ylim,
    width_mm = 89L
)

logger::log_info("[fig06] rendering Panel B (EGFR volcano)")
panel_b <- volcano_panel(
    egfr_dem,
    contrast_label = "EGFR",
    xlim = volcano_lims$xlim,
    ylim = volcano_lims$ylim,
    width_mm = 89L
)

# Panel D is stubbed pending the in-progress MetaLinksDB rebuild
# (task T064). Render a minimal placeholder so the composite layout
# is stable.
panel_d_stub <- function() {
    ggplot2::ggplot(
        data.frame(x = 0.5, y = 0.5,
                   label = "Panel D — pending MetaLinksDB rebuild"),
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$label)
    ) +
        ggplot2::geom_text(size = 2.5) +
        ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
        ggplot2::labs(title = NULL, x = NULL, y = NULL) +
        theme_bw_metabo(width_mm = 89L) +
        ggplot2::theme(
            axis.text   = ggplot2::element_blank(),
            axis.ticks  = ggplot2::element_blank(),
            panel.grid  = ggplot2::element_blank()
        )
}
panel_d <- panel_d_stub()

logger::log_info("[fig06] building Panel E PKN summaries (KRAS, EGFR)")
kras_summary <- case_study_pkn_summary(kras_dem, pkn)
egfr_summary <- case_study_pkn_summary(egfr_dem, pkn)

panel_e_kras <- gem_allosteric_panel(
    kras_summary, contrast_label = "KRAS", width_mm = 89L
)
panel_e_egfr <- gem_allosteric_panel(
    egfr_summary, contrast_label = "EGFR", width_mm = 89L
)

logger::log_info("[fig06] rendering Panel F sub-panels (location × direction)")
panel_f_kras_up <- subcellular_location_panel(
    kras_summary, contrast_label = "KRAS", direction = "up", width_mm = 89L
)
panel_f_kras_down <- subcellular_location_panel(
    kras_summary, contrast_label = "KRAS", direction = "down", width_mm = 89L
)
panel_f_egfr_up <- subcellular_location_panel(
    egfr_summary, contrast_label = "EGFR", direction = "up", width_mm = 89L
)
panel_f_egfr_down <- subcellular_location_panel(
    egfr_summary, contrast_label = "EGFR", direction = "down", width_mm = 89L
)

# ── Save individual panels ───────────────────────────────────────────────────

individual_panels <- list(
    panel_a            = panel_a,
    panel_b            = panel_b,
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

# ── Panel C placeholder check ────────────────────────────────────────────────

panel_c_path <- "figures/fig06-lungcancer-usecase/manual/panel_C_azelaic_acid.png"
panel_c_present <- file.exists(panel_c_path)
if (!panel_c_present) {
    logger::log_warn(paste0(
        "[fig06] Panel C placeholder absent (",
        panel_c_path, ") — composite assembled without Panel C. ",
        "Capture per manual/README.md when metabo.omnipathdb.org is ready."
    ))
}

# ── Composite (pipeline panels only) ────────────────────────────────────────
#
# Layout: 4 rows × 2 cols (matching the per-row panel groupings)
#   row 1: A   | B
#   row 2: E_kras  | E_egfr
#   row 3: F_kras_up | F_kras_down
#   row 4: F_egfr_up | F_egfr_down
# Panel D (stub) is emitted as a standalone file but not placed in
# the composite until it's a real plot. Panel C is injected by the
# mixed-source path once the manual asset lands.

# Embed manual tag labels per FR-024 capital-letter convention. Sub-
# panels of E and F share the parent letter; the placeholder tag for
# Panel D is set even though its content is a stub. compose_patchwork
# is called with tag_levels = NULL so the per-panel labs(tag=...)
# values are honoured rather than overwritten with sequential A–H.
composite_panels <- list(
    panel_a            + ggplot2::labs(tag = "A"),
    panel_b            + ggplot2::labs(tag = "B"),
    panel_e_kras       + ggplot2::labs(tag = "E"),
    panel_e_egfr       + ggplot2::labs(tag = ""),
    panel_f_kras_up    + ggplot2::labs(tag = "F"),
    panel_f_kras_down  + ggplot2::labs(tag = ""),
    panel_f_egfr_up    + ggplot2::labs(tag = ""),
    panel_f_egfr_down  + ggplot2::labs(tag = "")
)

pipeline_composite <- compose_patchwork(
    composite_panels,
    layout     = list(ncol = 2L),
    tag_levels = NULL
)

ggsave(
    filename = file.path(out_dir, "fig06-lungcancer-usecase-pipeline.pdf"),
    plot     = pipeline_composite,
    width    = 180,
    height   = 280,
    units    = "mm"
)
ggsave(
    filename = file.path(out_dir, "fig06-lungcancer-usecase-pipeline.svg"),
    plot     = pipeline_composite,
    width    = 180,
    height   = 280,
    units    = "mm"
)

# ── Caption ──────────────────────────────────────────────────────────────────

caption_info <- compose_caption(
    figure_id      = "fig06-lungcancer-usecase",
    composite_pdf  = file.path(
        out_dir, "fig06-lungcancer-usecase-pipeline.pdf"
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
        kind        = "panel-c-manual",
        path        = panel_c_path,
        fingerprint = if (panel_c_present) {
            digest::digest(file = panel_c_path, algo = "sha256")
        } else {
            NA_character_
        },
        present     = panel_c_present
    )
)

write_sidecar(
    artifact_id     = "fig06-lungcancer-usecase",
    artifact_path   = file.path(
        out_dir, "fig06-lungcancer-usecase-pipeline.pdf"
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
        panel_d_status           = "stub_pending_metalinks_rebuild",
        panel_c_status           = if (panel_c_present) {
            "present"
        } else {
            "placeholder_pending_webapp"
        },
        cosmos_pkn_source        = paste0(
            "vendored fixture from omnipath_metabo_case1/data/ ",
            "(pending T066 omnipath-client refactor)"
        )
    ),
    seed    = pipeline_seed(),
    caption = caption_info
)

logger::log_info("[fig06] build complete — outputs in {out_dir}")
