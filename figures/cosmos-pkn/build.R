# figures/cosmos-pkn/build.R
#
# Figure 5: COSMOS+ PKN — schematic + three pipeline panels.
#
# Active composite (post-2026-06-18 merge):
#   Top : regulation-types schematic (manual asset, full width)   — tag A
#   Row : Panel B — old-vs-COSMOS+ comparison by interaction type
#         Panel C — COSMOS+ interactions per compartment
#         Panel D — COSMOS+ entities per resource (grouped bars)
#   (B | C | D side-by-side; set composite_layout = 'comparison_fullwidth'
#    to put B full-width under the schematic with C | D beneath.)
#
# The MetaLinksDB 2.0 vs. COSMOS+ comparison is rendered as an individual
# artifact (panel_d.{pdf,svg}) but is NOT part of the composite (FR-011d,
# spec Session 2026-06-18).
#
# Data sources:
#   - data/vendored/cosmos/meta_network.RData    (old COSMOS PKN)
#   - data/vendored/cosmos/cosmos_plus_human.csv (COSMOS+)
#   - data/vendored/cosmos/cosmos_plus_mouse.csv (COSMOS+)
#   - custom_views.metalinksdb_relations on dev5 (MetaLinksDB 2.0)
#   - figures/cosmos-pkn/manual/regulation-types.png
#     (schematic; PDF + PNG variants both present)

suppressPackageStartupMessages({
    library(metabo.figures)
    library(ggplot2)
    library(patchwork)
})

setup_pipeline_log('build:cosmos-pkn')
set.seed(pipeline_seed())

out_dir <- 'figures/cosmos-pkn/out'
fs::dir_create(out_dir)

# ── Data loading ─────────────────────────────────────────────────────────────

logger::log_info('[cosmos-pkn] loading old COSMOS PKN (vendored)')
old_pkn <- cosmos_old_pkn()

logger::log_info('[cosmos-pkn] loading COSMOS+ data from vendored CSVs')
cosmos_plus <- cosmos_plus_data()

logger::log_info('[cosmos-pkn] querying MetaLinksDB 2.0 by GtP protein class from dev5')
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
    'cosmos-pkn',
    metalinks_v2_sql,
    facet = 'metalinks'
)

# ── Panel rendering ──────────────────────────────────────────────────────────

# COSMOS+ protein-interaction PKN is sourced from OmniPath and is human-only.
# These are the authoritative edge counts from the omnipath_metabo builders:
#   from omnipath_metabo.datasets.cosmos._build import build_ppi, build_grn
#   len(pd.DataFrame(build_ppi().network))  ->  Signaling (PPI)
#   len(pd.DataFrame(build_grn().network))  ->  GRN
# Captured 2026-06-18 on beauty (~/dev/omnipath-metabo, .venv). Re-run those
# builders to refresh if the upstream OmniPath release changes.
protein_omnipath_counts <- c(
    'Signaling (PPI)' = 34367,
    'GRN'             = 45976
)

logger::log_info('[cosmos-pkn] rendering Panel A (faceted old COSMOS vs. COSMOS+)')
panel_a <- fig04_cosmos_comparison_panel(
    old_pkn_tally               = old_pkn,
    cosmos_plus_by_type_species = cosmos_plus$by_type_species,
    protein_omnipath_counts     = protein_omnipath_counts,
    width_mm                    = 120L
)

logger::log_info('[cosmos-pkn] rendering Panel B (COSMOS+ compartments)')
panel_b <- fig04_compartment_panel(
    cosmos_plus_by_compartment = cosmos_plus$by_compartment,
    width_mm                   = 89L
)

logger::log_info('[cosmos-pkn] rendering Panel C (COSMOS+ resource contributions)')
panel_c <- fig04_resource_contribution_panel(
    cosmos_plus_by_resource = cosmos_plus$by_resource,
    width_mm                = 89L
)

logger::log_info('[cosmos-pkn] rendering Panel C (grouped) — entities per resource as dodged bars')
panel_c_split <- fig04_resource_contribution_panel(
    cosmos_plus_by_resource = cosmos_plus$by_resource_split,
    width_mm                = 89L,
    position                = 'dodge'
)

logger::log_info('[cosmos-pkn] rendering Panel D (MetaLinksDB 2.0 vs. COSMOS+)')
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

# Composite (post-2026-06-18 merge): schematic (A) + comparison (B) +
# compartments (C) + resources (D). The MetaLinksDB panel (panel_d) stays a
# standalone artifact and is NOT composited. Aliases below carry the composite
# panel lettering so the assembly reads in panel order.
panel_comparison  <- panel_a        # FR-011a — composite Panel B
panel_compartment <- panel_b        # FR-011b — composite Panel C
panel_resource    <- panel_c_split  # FR-011c — composite Panel D

# ── Save individual panels ───────────────────────────────────────────────────

panel_dims <- list(
    panel_a       = c(140, 185),
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

# ── Composite (schematic A on top + B | C | D data row) ─────────────────────
#
# Post-2026-06-18 merge: four panels — the regulation-types schematic (A,
# full width, top) plus the three COSMOS+ data panels
#   B = old-vs-COSMOS+ comparison   (panel_comparison)
#   C = compartment coverage        (panel_compartment)
#   D = resource contributions      (panel_resource)
# The three data panels share a portrait coord_flip shape, so they sit
# side-by-side in one row rather than stretching the comparison across a
# full-width slot. Set composite_layout <- 'comparison_fullwidth' to instead
# place the comparison full width directly under the schematic, with C | D as
# a bottom row. The chosen layout is recorded in the sidecar (active_composite).
#
# The regulation-types schematic is embedded as a raster grob via png + grid
# + patchwork::wrap_elements (pdftools/magick are not in the dependency
# surface, so the PNG variant is used; the source PDF is kept alongside as the
# canonical asset). A missing schematic is non-fatal: the build warns and
# composes the three data panels on their own (FR-011, T056f).

composite_layout <- 'schematic_bcd_row'

schematic_png_path <- 'figures/cosmos-pkn/manual/regulation-types.png'

schematic_panel <- if (file.exists(schematic_png_path)) {
    schematic_grob <- grid::rasterGrob(
        png::readPNG(schematic_png_path, native = TRUE),
        interpolate = TRUE
    )
    patchwork::wrap_elements(full = schematic_grob) +
        ggplot2::labs(tag = 'A')
} else {
    logger::log_warn(paste0(
        '[cosmos-pkn] schematic PNG missing at ', schematic_png_path,
        ' — composite will be assembled without the top schematic panel.'
    ))
    NULL
}

# Data-row order (feedback 2026-06-18): compartment | resource | comparison,
# i.e. the old-vs-COSMOS+ comparison moves to the rightmost slot (D) and the
# compartment + resource panels shift left into B | C.
data_row <- (panel_compartment + ggplot2::labs(tag = 'B')) |
            (panel_resource    + ggplot2::labs(tag = 'C')) |
            (panel_comparison  + ggplot2::labs(tag = 'D'))

if (!is.null(schematic_panel)) {
    if (identical(composite_layout, 'comparison_fullwidth')) {
        # A: schematic on top. B: comparison full width. C | D: bottom row.
        cd_row <- (panel_compartment + ggplot2::labs(tag = 'C')) |
                  (panel_resource    + ggplot2::labs(tag = 'D'))
        pipeline_composite <-
            (schematic_panel /
             (panel_comparison + ggplot2::labs(tag = 'B')) /
             cd_row) +
            patchwork::plot_layout(heights = c(1.1, 1.0, 1.0))
    } else {
        # A: schematic on top. B | C | D: single data row underneath.
        pipeline_composite <- (schematic_panel / data_row) +
            patchwork::plot_layout(heights = c(0.9, 1.0))
    }
} else {
    # No schematic: the three data panels alone, re-tagged A | B | C in the
    # same compartment | resource | comparison order as the data row above.
    pipeline_composite <- (
        (panel_compartment + ggplot2::labs(tag = 'A')) |
        (panel_resource    + ggplot2::labs(tag = 'B')) |
        (panel_comparison  + ggplot2::labs(tag = 'C'))
    )
}

# Bump panel-letter (tag) size across the composite. patchwork's `&`
# applies the theme to every nested plot, including wrap_elements
# patches like the schematic — so A, B, C, D all render at the same
# visually-large weight.
pipeline_composite <- pipeline_composite &
    ggplot2::theme(
        plot.tag = ggplot2::element_text(size = 16, face = 'bold')
    )

logger::log_info(paste0(
    '[cosmos-pkn] assembled composite (',
    if (is.null(schematic_panel)) 'pipeline-only, schematic missing'
    else if (identical(composite_layout, 'comparison_fullwidth'))
        'schematic on top, comparison full-width, C | D bottom row'
    else 'schematic on top, B | C | D data row',
    ')'
))

composite_width_mm  <- 180L
composite_height_mm <- if (is.null(schematic_panel)) {
    110L
} else if (identical(composite_layout, 'comparison_fullwidth')) {
    260L
} else {
    200L
}

ggsave(
    filename = file.path(out_dir, 'cosmos-pkn.pdf'),
    plot     = pipeline_composite,
    width    = composite_width_mm,
    height   = composite_height_mm,
    units    = 'mm'
)
ggsave(
    filename = file.path(out_dir, 'cosmos-pkn.svg'),
    plot     = pipeline_composite,
    width    = composite_width_mm,
    height   = composite_height_mm,
    units    = 'mm'
)

# Retain the legacy pipeline-only artifacts for backward compatibility
# (downstream consumers reference cosmos-pkn-pipeline.pdf).
file.copy(
    file.path(out_dir, 'cosmos-pkn.pdf'),
    file.path(out_dir, 'cosmos-pkn-pipeline.pdf'),
    overwrite = TRUE
)
file.copy(
    file.path(out_dir, 'cosmos-pkn.svg'),
    file.path(out_dir, 'cosmos-pkn-pipeline.svg'),
    overwrite = TRUE
)

# ── Caption ──────────────────────────────────────────────────────────────────

caption_info <- compose_caption(
    figure_id      = 'cosmos-pkn',
    composite_pdf  = file.path(out_dir, 'cosmos-pkn.pdf'),
    caption_source = 'figures/cosmos-pkn/caption.tex',
    out_dir        = out_dir,
    panel_count    = if (!is.null(schematic_panel)) 4L else 3L
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
    artifact_id     = 'cosmos-pkn',
    artifact_path   = file.path(out_dir, 'cosmos-pkn.pdf'),
    deployments     = list(dep5$deployment),
    manifests       = list(dep5$manifest),
    script_path     = 'figures/cosmos-pkn/build.R',
    queries         = list(query_record(metalinks_v2_types)),
    external_inputs = c(
        list(list(
            kind               = 'old-cosmos-pkn',
            source_pkg         = attr(old_pkn, 'source_pkg'),
            fingerprint        = attr(old_pkn, 'fingerprint'),
            species_assumption = 'human-only-or-unspecified'
        )),
        cosmos_plus$external_inputs,
        if (!is.null(schematic_input)) list(schematic_input) else list()
    ),
    parameters = list(
        old_cosmos_species_assumption = 'human-only-or-unspecified',
        metalinks_deployment          = 'dev5',
        metalinks_view                = 'custom_views.metalinksdb_relations',
        metalinks_panel               = 'standalone (not composited)',
        comparison_facets             = paste(
            'Protein interactions (Old COSMOS|Signaling(PPI)|GRN)',
            'Metabolite interactions (Old COSMOS|Transporter|Allosteric|Receptor)',
            'Metabolic reactions (Old COSMOS|New COSMOS)',
            sep = '; '
        ),
        comparison_protein_source     = 'omnipath_metabo build_ppi()/build_grn() (human-only)',
        comparison_protein_counts     = list(
            `Signaling (PPI)` = unname(protein_omnipath_counts[['Signaling (PPI)']]),
            GRN               = unname(protein_omnipath_counts[['GRN']])
        ),
        active_composite              = composite_layout,
        composite_panels              = paste(
            'A:regulation-types-schematic', 'B:compartment-coverage',
            'C:resource-contributions', 'D:old-vs-COSMOS+ comparison',
            sep = '; '
        ),
        panel_c_position              = 'dodge',
        composite_dims_mm             = list(
            width  = composite_width_mm,
            height = composite_height_mm
        )
    ),
    seed    = pipeline_seed(),
    caption = caption_info
)

logger::log_info('[cosmos-pkn] build complete — outputs in {out_dir}')
