# figures/fig05-cosmos-pkn/build.R
#
# Figure 5: COSMOS+ PKN — schematic + two pipeline panels.
#
# Active composite (this version):
#   Top    : regulation-types schematic (manual asset, full width)
#   Middle : Panel B — COSMOS+ interactions per compartment
#   Bottom : Panel C — COSMOS+ entities per resource (grouped bars)
#
# Panels A (old vs. COSMOS+ comparison) and D (MetaLinksDB vs.
# COSMOS+) are still rendered as individual artifacts but are not
# part of the active composite (colleague's call; see Figure 5
# session 2026-06-16).
#
# Data sources:
#   - inst/extdata/cosmos/meta_network.RData    (old COSMOS PKN)
#   - inst/extdata/cosmos/cosmos_plus_human.csv (COSMOS+)
#   - inst/extdata/cosmos/cosmos_plus_mouse.csv (COSMOS+)
#   - custom_views.metalinksdb_relations on dev5 (MetaLinksDB 2.0)
#   - figures/fig05-cosmos-pkn/manual/regulation-types.png
#     (schematic; PDF + PNG variants both present)

suppressPackageStartupMessages({
    library(metabo.figures)
    library(ggplot2)
    library(patchwork)
})

setup_pipeline_log('build:fig05-cosmos-pkn')
set.seed(pipeline_seed())

out_dir <- 'figures/fig05-cosmos-pkn/out'
fs::dir_create(out_dir)

# ── Data loading ─────────────────────────────────────────────────────────────

logger::log_info('[fig05] loading old COSMOS PKN (vendored)')
old_pkn <- cosmos_old_pkn()

logger::log_info('[fig05] loading COSMOS+ data from vendored CSVs')
cosmos_plus <- cosmos_plus_data()

logger::log_info('[fig05] querying MetaLinksDB 2.0 by GtP protein class from dev5')
dep5 <- deployment_provenance('dev5')

# Join metalinksdb_relations with metalinksdb_protein_annotations to get
# per-interaction GtP functional class, then map to canonical panel categories.
# A protein with multiple GtP classes is counted in each matching category;
# COUNT(DISTINCT ...) prevents double-counting within a category.
metalinks_v2_sql <- paste(
    'SELECT',
    '  CASE',
    '    WHEN gtp.gtp_class LIKE \'Transporter%\'           THEN \'Transport\'',
    '    WHEN gtp.gtp_class LIKE \'Gpcr%\'',
    '      OR gtp.gtp_class LIKE \'Catalytic Receptor%\'',
    '      OR gtp.gtp_class LIKE \'Vgic%\'',
    '      OR gtp.gtp_class LIKE \'Lgic%\'                  THEN \'Ligand receptor\'',
    '    WHEN gtp.gtp_class LIKE \'Enzyme%\'                THEN \'Catalysis\'',
    '    WHEN gtp.gtp_class LIKE \'Nuclear Hormone%\'       THEN \'Gene regulation\'',
    '    ELSE \'Other\'',
    '  END AS interaction_type,',
    '  COUNT(DISTINCT r.compound_canonical_id || \'::\'',
    '       || r.protein_uniprot) AS n_interactions',
    # dev5 schema change: metalinksdb_protein_annotations.protein_uniprot
    # now holds the protein_canonical_id (which resolves to Entrez gene
    # IDs, not UniProt accessions). Join via the entity-UUID instead so
    # the relation rows still match their annotation rows.
    'FROM custom_views.metalinksdb_relations r',
    'JOIN custom_views.metalinksdb_protein_annotations a',
    '  ON r.protein_entity_id = a.protein_entity_id',
    'CROSS JOIN LATERAL unnest(a.gtp_functional_classes) AS gtp(gtp_class)',
    'WHERE r.compound_canonical_id IS NOT NULL',
    '  AND r.protein_uniprot IS NOT NULL',
    'GROUP BY 1',
    'ORDER BY n_interactions DESC'
)

metalinks_v2_types <- pg_query_panel(
    'fig05-cosmos-pkn',
    metalinks_v2_sql,
    facet = 'metalinks'
)

# ── Panel rendering ──────────────────────────────────────────────────────────

logger::log_info('[fig05] rendering Panel A (old COSMOS vs. COSMOS+)')
panel_a <- fig04_cosmos_comparison_panel(
    old_pkn_tally               = old_pkn,
    cosmos_plus_by_type_species = cosmos_plus$by_type_species,
    width_mm                    = 120L
)

logger::log_info('[fig05] rendering Panel B (COSMOS+ compartments)')
panel_b <- fig04_compartment_panel(
    cosmos_plus_by_compartment = cosmos_plus$by_compartment,
    width_mm                   = 89L
)

logger::log_info('[fig05] rendering Panel C (COSMOS+ resource contributions)')
panel_c <- fig04_resource_contribution_panel(
    cosmos_plus_by_resource = cosmos_plus$by_resource,
    width_mm                = 89L
)

logger::log_info('[fig05] rendering Panel C (grouped) — entities per resource as dodged bars')
panel_c_split <- fig04_resource_contribution_panel(
    cosmos_plus_by_resource = cosmos_plus$by_resource_split,
    width_mm                = 89L,
    position                = 'dodge'
)

logger::log_info('[fig05] rendering Panel D (MetaLinksDB 2.0 vs. COSMOS+)')
panel_d <- fig04_metalinks_cosmos_panel(
    metalinks_counts            = metalinks_v2_types,
    cosmos_plus_by_type_species = cosmos_plus$by_type_species,
    width_mm                    = 120L
)

# All individual panels including both C variants
panels_all <- list(
    panel_a       = panel_a,
    panel_b       = panel_b,
    panel_c       = panel_c,
    panel_c_split = panel_c_split,
    panel_d       = panel_d
)

# Main composite: Panel B and C (split) only; A and D saved individually
panels <- list(
    panel_b = panel_b,
    panel_c = panel_c_split
)

# ── Save individual panels ───────────────────────────────────────────────────

panel_dims <- list(
    panel_a       = c(120, 110),
    panel_b       = c(89,  140),
    panel_c       = c(89,  130),
    panel_c_split = c(89,  160),
    panel_d       = c(120, 110)
)

for (name in names(panels_all)) {
    dims <- panel_dims[[name]]
    ggsave(
        filename = file.path(out_dir, paste0(name, '.pdf')),
        plot     = panels_all[[name]],
        width    = dims[1L],
        height   = dims[2L],
        units    = 'mm'
    )
    ggsave(
        filename = file.path(out_dir, paste0(name, '.svg')),
        plot     = panels_all[[name]],
        width    = dims[1L],
        height   = dims[2L],
        units    = 'mm'
    )
}

# ── Composite (schematic on top, Panel B, Panel C stacked) ──────────────────
#
# The regulation-types schematic is embedded as a raster grob via png
# + grid + patchwork::wrap_elements. PDF is preferable for vector
# graphics but neither pdftools nor magick are in the project's R
# dependency surface, so we use the PNG variant; the source PDF is
# kept alongside as the canonical asset.

schematic_png_path <- 'figures/fig05-cosmos-pkn/manual/regulation-types.png'

schematic_panel <- if (file.exists(schematic_png_path)) {
    schematic_grob <- grid::rasterGrob(
        png::readPNG(schematic_png_path, native = TRUE),
        interpolate = TRUE
    )
    patchwork::wrap_elements(full = schematic_grob) +
        ggplot2::labs(tag = 'A')
} else {
    logger::log_warn(paste0(
        '[fig05] schematic PNG missing at ', schematic_png_path,
        ' — composite will be assembled without the top schematic panel.'
    ))
    NULL
}

if (!is.null(schematic_panel)) {
    # A: schematic at 100% width on top.
    # B | C: side-by-side bottom row.
    bottom_row <- (panel_b       + ggplot2::labs(tag = 'B')) |
                  (panel_c_split + ggplot2::labs(tag = 'C'))
    pipeline_composite <- (schematic_panel / bottom_row) +
        patchwork::plot_layout(heights = c(1.25, 1.0))
} else {
    pipeline_composite <- (
        (panel_b       + ggplot2::labs(tag = 'A')) |
        (panel_c_split + ggplot2::labs(tag = 'B'))
    )
}

# Bump panel-letter (tag) size across the composite. patchwork's `&`
# applies the theme to every nested plot, including wrap_elements
# patches like the schematic — so A, B, C all render at the same
# visually-large weight.
pipeline_composite <- pipeline_composite &
    ggplot2::theme(
        plot.tag = ggplot2::element_text(size = 16, face = 'bold')
    )

logger::log_info(paste0(
    '[fig05] assembled composite (',
    if (is.null(schematic_panel)) 'pipeline-only, schematic missing'
    else 'schematic on top, B | C bottom row',
    ')'
))

composite_width_mm  <- 180L
composite_height_mm <- if (!is.null(schematic_panel)) 220L else 110L

ggsave(
    filename = file.path(out_dir, 'fig05-cosmos-pkn.pdf'),
    plot     = pipeline_composite,
    width    = composite_width_mm,
    height   = composite_height_mm,
    units    = 'mm'
)
ggsave(
    filename = file.path(out_dir, 'fig05-cosmos-pkn.svg'),
    plot     = pipeline_composite,
    width    = composite_width_mm,
    height   = composite_height_mm,
    units    = 'mm'
)

# Retain the legacy pipeline-only artifacts for backward compatibility
# (downstream consumers reference fig05-cosmos-pkn-pipeline.pdf).
file.copy(
    file.path(out_dir, 'fig05-cosmos-pkn.pdf'),
    file.path(out_dir, 'fig05-cosmos-pkn-pipeline.pdf'),
    overwrite = TRUE
)
file.copy(
    file.path(out_dir, 'fig05-cosmos-pkn.svg'),
    file.path(out_dir, 'fig05-cosmos-pkn-pipeline.svg'),
    overwrite = TRUE
)

# ── Caption ──────────────────────────────────────────────────────────────────

caption_info <- compose_caption(
    figure_id      = 'fig05-cosmos-pkn',
    composite_pdf  = file.path(out_dir, 'fig05-cosmos-pkn.pdf'),
    caption_source = 'figures/fig05-cosmos-pkn/caption.tex',
    out_dir        = out_dir,
    panel_count    = if (!is.null(schematic_panel)) 3L else 2L
)

# ── Provenance sidecar ────────────────────────────────────────────────────────

schematic_input <- if (!is.null(schematic_panel)) {
    list(
        kind        = 'schematic-regulation-types',
        path        = schematic_png_path,
        fingerprint = digest::digest(file = schematic_png_path,
                                     algo = 'sha256'),
        source      = 'manual asset (regulation-types)'
    )
} else {
    NULL
}

write_sidecar(
    artifact_id     = 'fig05-cosmos-pkn',
    artifact_path   = file.path(out_dir, 'fig05-cosmos-pkn.pdf'),
    deployments     = list(dep5$deployment),
    manifests       = list(dep5$manifest),
    script_path     = 'figures/fig05-cosmos-pkn/build.R',
    queries         = list(query_record(metalinks_v2_types)),
    external_inputs = c(
        list(list(
            kind               = 'old-cosmos-pkn',
            path               = attr(old_pkn, 'source_path'),
            fingerprint        = attr(old_pkn, 'fingerprint'),
            species_assumption = 'human-only-or-unspecified'
        )),
        cosmos_plus$external_inputs,
        if (!is.null(schematic_input)) list(schematic_input) else list()
    ),
    parameters = list(
        old_cosmos_species_assumption = 'human-only-or-unspecified',
        panel_d_metalinks_deployment  = 'dev5',
        panel_d_metalinks_view        = 'custom_views.metalinksdb_relations',
        active_composite              = 'schematic_b_c_vertical',
        panel_c_position              = 'dodge',
        composite_dims_mm             = list(
            width  = composite_width_mm,
            height = composite_height_mm
        )
    ),
    seed    = pipeline_seed(),
    caption = caption_info
)

logger::log_info('[fig05] build complete — outputs in {out_dir}')
