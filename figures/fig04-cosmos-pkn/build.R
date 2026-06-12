# figures/fig04-cosmos-pkn/build.R
#
# Figure 4: COSMOS+ PKN analysis — four pipeline panels.
#
# Panel A: species-aware comparison of old COSMOS PKN vs. COSMOS+ by interaction type
# Panel B: COSMOS+ interactions per annotated subcellular compartment
# Panel C: entity and interaction counts per contributing resource in COSMOS+
# Panel D: MetaLinksDB 2.0 vs. COSMOS+ comparison by interaction type
#
# Three manual schematics (pkn-to-binary-network, regulation-types,
# moon-activity-inference) must exist under manual/ before the mixed-source
# composite assembly step. Missing assets trigger a warning and skip the full
# composite (spec Edge Case).
#
# Data sources:
#   - inst/extdata/cosmos/meta_network.RData    (old COSMOS PKN, vendored)
#   - inst/extdata/cosmos/cosmos_plus_human.csv (COSMOS+, vendored via T052a)
#   - inst/extdata/cosmos/cosmos_plus_mouse.csv (COSMOS+, vendored via T052a)
#   - custom_views.metalinksdb_relations on dev4 (MetaLinksDB 2.0, Panel D)

suppressPackageStartupMessages({
    library(metabo.figures)
    library(ggplot2)
})

setup_pipeline_log('build:fig04-cosmos-pkn')
set.seed(pipeline_seed())

out_dir <- 'figures/fig04-cosmos-pkn/out'
fs::dir_create(out_dir)

# ── Data loading ─────────────────────────────────────────────────────────────

logger::log_info('[fig04] loading old COSMOS PKN (vendored)')
old_pkn <- cosmos_old_pkn()

logger::log_info('[fig04] loading COSMOS+ data from vendored CSVs')
cosmos_plus <- cosmos_plus_data()

logger::log_info('[fig04] querying MetaLinksDB 2.0 by GtP protein class from dev4')
dep4 <- deployment_provenance('dev4')

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
    'FROM custom_views.metalinksdb_relations r',
    'JOIN custom_views.metalinksdb_protein_annotations a',
    '  ON r.protein_uniprot = a.protein_uniprot',
    'CROSS JOIN LATERAL unnest(a.gtp_functional_classes) AS gtp(gtp_class)',
    'WHERE r.compound_canonical_id IS NOT NULL',
    '  AND r.protein_uniprot IS NOT NULL',
    'GROUP BY 1',
    'ORDER BY n_interactions DESC'
)

metalinks_v2_types <- pg_query_panel(
    'fig04-cosmos-pkn',
    metalinks_v2_sql,
    facet = 'metalinks'
)

# ── Panel rendering ──────────────────────────────────────────────────────────

logger::log_info('[fig04] rendering Panel A (old COSMOS vs. COSMOS+)')
panel_a <- fig04_cosmos_comparison_panel(
    old_pkn_tally               = old_pkn,
    cosmos_plus_by_type_species = cosmos_plus$by_type_species,
    width_mm                    = 120L
)

logger::log_info('[fig04] rendering Panel B (COSMOS+ compartments)')
panel_b <- fig04_compartment_panel(
    cosmos_plus_by_compartment = cosmos_plus$by_compartment,
    width_mm                   = 89L
)

logger::log_info('[fig04] rendering Panel C (COSMOS+ resource contributions)')
panel_c <- fig04_resource_contribution_panel(
    cosmos_plus_by_resource = cosmos_plus$by_resource,
    width_mm                = 89L
)

logger::log_info('[fig04] rendering Panel D (MetaLinksDB 2.0 vs. COSMOS+)')
panel_d <- fig04_metalinks_cosmos_panel(
    metalinks_counts            = metalinks_v2_types,
    cosmos_plus_by_type_species = cosmos_plus$by_type_species,
    width_mm                    = 120L
)

panels <- list(
    panel_a = panel_a,
    panel_b = panel_b,
    panel_c = panel_c,
    panel_d = panel_d
)

# ── Save individual panels ───────────────────────────────────────────────────

for (name in names(panels)) {
    ggsave(
        filename = file.path(out_dir, paste0(name, '.pdf')),
        plot     = panels[[name]],
        width    = 120,
        height   = 90,
        units    = 'mm'
    )
    ggsave(
        filename = file.path(out_dir, paste0(name, '.svg')),
        plot     = panels[[name]],
        width    = 120,
        height   = 90,
        units    = 'mm'
    )
}

# ── Pipeline panel composite ─────────────────────────────────────────────────

pipeline_composite <- compose_patchwork(panels, layout = list(ncol = 2L))

ggsave(
    filename = file.path(out_dir, 'fig04-cosmos-pkn-pipeline.pdf'),
    plot     = pipeline_composite,
    width    = 240,
    height   = 180,
    units    = 'mm'
)
ggsave(
    filename = file.path(out_dir, 'fig04-cosmos-pkn-pipeline.svg'),
    plot     = pipeline_composite,
    width    = 240,
    height   = 180,
    units    = 'mm'
)

# ── Full composite (pipeline panels + manual schematics) ─────────────────────

manual_dir <- 'figures/fig04-cosmos-pkn/manual'
schematic_slugs <- c(
    'pkn-to-binary-network',
    'regulation-types',
    'moon-activity-inference'
)
schematic_pdfs <- file.path(manual_dir, paste0(schematic_slugs, '.pdf'))
missing_schematics <- schematic_pdfs[!file.exists(schematic_pdfs)]

if (length(missing_schematics) > 0L) {
    logger::log_warn(paste0(
        '[fig04] manual schematic(s) absent — skipping full composite. ',
        'Missing: ', paste(basename(missing_schematics), collapse = ', ')
    ))
} else {
    logger::log_info('[fig04] assembling full composite (pipeline + manual schematics)')
    compose_mixed_source(
        pipeline_pdf = file.path(out_dir, 'fig04-cosmos-pkn-pipeline.pdf'),
        manual_pdfs  = schematic_pdfs,
        out_pdf      = file.path(out_dir, 'fig04-cosmos-pkn.pdf'),
        out_svg      = file.path(out_dir, 'fig04-cosmos-pkn.svg')
    )
}

# ── Caption ──────────────────────────────────────────────────────────────────

caption_info <- compose_caption(
    figure_id      = 'fig04-cosmos-pkn',
    composite_pdf  = file.path(out_dir, 'fig04-cosmos-pkn-pipeline.pdf'),
    caption_source = 'figures/fig04-cosmos-pkn/caption.tex',
    out_dir        = out_dir,
    panel_count    = 4L
)

# ── Provenance sidecar ────────────────────────────────────────────────────────

write_sidecar(
    artifact_id     = 'fig04-cosmos-pkn',
    artifact_path   = file.path(out_dir, 'fig04-cosmos-pkn-pipeline.pdf'),
    deployments     = list(dep4$deployment),
    manifests       = list(dep4$manifest),
    script_path     = 'figures/fig04-cosmos-pkn/build.R',
    queries         = list(query_record(metalinks_v2_types)),
    external_inputs = c(
        list(list(
            kind               = 'old-cosmos-pkn',
            path               = attr(old_pkn, 'source_path'),
            fingerprint        = attr(old_pkn, 'fingerprint'),
            species_assumption = 'human-only-or-unspecified'
        )),
        cosmos_plus$external_inputs
    ),
    parameters = list(
        old_cosmos_species_assumption = 'human-only-or-unspecified',
        panel_d_metalinks_deployment  = 'dev4',
        panel_d_metalinks_view        = 'custom_views.metalinksdb_relations'
    ),
    seed    = pipeline_seed(),
    caption = caption_info
)

logger::log_info('[fig04] build complete — outputs in {out_dir}')
