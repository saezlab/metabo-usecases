# figures/fig04-cosmos-pkn/build.R
#
# Figure 4: COSMOS+ prior-knowledge network — scope extension,
# compartment coverage, resource contributions, and comparison with
# MetaLinksDB 2.0.

suppressPackageStartupMessages({
    library(metabo.figures)
    library(ggplot2)
})

setup_pipeline_log('build:fig04-cosmos-pkn')
set.seed(pipeline_seed())

out_dir <- 'figures/fig04-cosmos-pkn/out'
fs::dir_create(out_dir)

# Load old COSMOS PKN
logger::log_info('Loading old COSMOS PKN')
old_pkn <- cosmos_old_pkn()

# Load COSMOS+ data (requires T052a CSVs to be vendored)
logger::log_info('Loading COSMOS+ data')
cosmos_plus <- cosmos_plus_data()

# Load MetaLinksDB v2 interaction-type counts from dev4 for Panel D
dep4 <- deployment_provenance('dev4')
metalinks_v2_sql <- paste(
    'SELECT',
    '  coalesce(rt.relation_type, \'interaction\') AS interaction_type,',
    '  COUNT(DISTINCT r.compound_canonical_id || \'::\'',
    '       || r.protein_uniprot) AS n_interactions',
    'FROM custom_views.metalinksdb_relations r',
    'LEFT JOIN LATERAL unnest(r.relation_types) AS rt(relation_type) ON TRUE',
    'WHERE r.compound_canonical_id IS NOT NULL',
    '  AND r.protein_uniprot IS NOT NULL',
    'GROUP BY coalesce(rt.relation_type, \'interaction\')',
    'ORDER BY n_interactions DESC'
)

logger::log_info('Querying MetaLinksDB v2 interaction types from dev4')
metalinks_v2_types <- pg_query_panel(
    'fig04-cosmos-pkn',
    metalinks_v2_sql,
    facet = 'metalinks'
)

# Render panels
logger::log_info('Rendering Panel A (old COSMOS vs. COSMOS+)')
panel_a <- fig04_cosmos_comparison_panel(
    old_pkn_tally            = old_pkn,
    cosmos_plus_by_type_species = cosmos_plus$by_type_species,
    width_mm                 = 120L
)

logger::log_info('Rendering Panel B (compartment coverage)')
panel_b <- fig04_compartment_panel(
    cosmos_plus_by_compartment = cosmos_plus$by_compartment,
    width_mm                   = 89L
)

logger::log_info('Rendering Panel C (resource contributions)')
panel_c <- fig04_resource_contribution_panel(
    cosmos_plus_by_resource = cosmos_plus$by_resource,
    width_mm                = 89L
)

panels <- list(
    panel_a = panel_a,
    panel_b = panel_b,
    panel_c = panel_c
)

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

# Composite of the three pipeline panels
composite <- compose_patchwork(panels, layout = list(ncol = 2L))
ggsave(
    filename = file.path(out_dir, 'fig04-cosmos-pkn.pdf'),
    plot     = composite,
    width    = 240,
    height   = 190,
    units    = 'mm'
)
ggsave(
    filename = file.path(out_dir, 'fig04-cosmos-pkn.svg'),
    plot     = composite,
    width    = 240,
    height   = 190,
    units    = 'mm'
)

caption_info <- compose_caption(
    figure_id     = 'fig04-cosmos-pkn',
    composite_pdf = file.path(out_dir, 'fig04-cosmos-pkn.pdf'),
    caption_source = 'figures/fig04-cosmos-pkn/caption.tex',
    out_dir       = out_dir,
    panel_count   = 3L
)

write_sidecar(
    artifact_id    = 'fig04-cosmos-pkn',
    artifact_path  = file.path(out_dir, 'fig04-cosmos-pkn.pdf'),
    deployments    = list(dep4$deployment),
    manifests      = list(dep4$manifest),
    script_path    = 'figures/fig04-cosmos-pkn/build.R',
    queries        = list(query_record(metalinks_v2_types)),
    external_inputs = c(
        cosmos_plus$external_inputs,
        list(list(
            kind        = 'old-cosmos-pkn',
            path        = attr(old_pkn, 'source_path'),
            fingerprint = attr(old_pkn, 'fingerprint'),
            species_assumption = 'human-only-or-unspecified'
        ))
    ),
    parameters = list(
        old_cosmos_species_assumption = 'human-only-or-unspecified'
    ),
    seed    = pipeline_seed(),
    caption = caption_info
)

logger::log_info('Figure 4 build complete')
