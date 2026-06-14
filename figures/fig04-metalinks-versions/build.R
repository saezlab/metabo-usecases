# figures/fig04-metalinks-versions/build.R
#
# First implementation slice of the Figure 3 pipeline. Uses the current
# MetaLinksDB v2 combined network contract, the OmnipathR MetaLinksDB v1
# SQLite baseline, and any vendored baseline snapshots already present.

suppressPackageStartupMessages({
    library(metabo.figures)
    library(ggplot2)
})

setup_pipeline_log('build:fig04-metalinks-versions')
set.seed(pipeline_seed())

out_dir <- 'figures/fig04-metalinks-versions/out'
fs::dir_create(out_dir)

dep5 <- deployment_provenance('dev5')

metalinks_v2_sql <- paste(
    'select',
    '  r.compound_canonical_id as hmdb,',
    '  r.protein_uniprot as uniprot,',
    '  s.source as source,',
    "  coalesce(rt.relation_type, 'interaction') as relation_type,",
    '  r.source_count,',
    '  case when r.pubmed_ids is null then null else cardinality(r.pubmed_ids) end as citation_count,',
    '  r.best_pchembl_value as affinity_value,',
    '  coalesce(ca.lipid_sub_class, ca.lipid_main_class, ca.lipid_category) as metabolite_class,',
    '  coalesce(pa.gtp_functional_classes[1], pa.uniprot_protein_families[1]) as protein_class',
    'from custom_views.metalinksdb_relations r',
    'left join lateral unnest(r.sources) as s(source) on true',
    'left join lateral unnest(r.relation_types) as rt(relation_type) on true',
    'left join custom_views.metalinksdb_compound_annotations ca',
    '  on ca.compound_entity_id = r.compound_entity_id',
    'left join custom_views.metalinksdb_protein_annotations pa',
    '  on pa.protein_entity_id = r.protein_entity_id',
    'where r.compound_canonical_id is not null and r.protein_uniprot is not null'
)

logger::log_info('Querying MetaLinksDB v2 from dev5')
metalinks_v2_rows <- pg_query_panel(
    'fig04-metalinks-versions',
    metalinks_v2_sql
)
metalinks_v2 <- normalize_mpi_resource(
    resource = 'MetaLinksDB v2.0',
    interactions = metalinks_v2_rows,
    metabolite_key = 'hmdb',
    protein_key = 'uniprot',
    source_col = 'source',
    relation_type_col = 'relation_type',
    metabolite_class = metalinks_v2_rows,
    metabolite_class_key = 'hmdb',
    metabolite_class_col = 'metabolite_class',
    protein_class = metalinks_v2_rows,
    protein_class_key = 'uniprot',
    protein_class_col = 'protein_class',
    evidence_cols = list(
        source_count = 'source_count',
        citation_count = 'citation_count',
        affinity_value = 'affinity_value',
        curation_mode = NULL
    ),
    interaction_definition = 'one HMDB-UniProt-source-relation row from custom_views.metalinksdb_relations'
)
attr(metalinks_v2, 'deployment') <- 'dev5'

logger::log_info('Loading MetaLinksDB v1 baseline')
v1 <- metalinks_v1_snapshot()

baseline_resources <- c('CellPhoneDB', 'scConnect', 'STITCH')
baselines <- list()
excluded_resources <- character(0)

for (resource in baseline_resources) {
    snapshot <- load_vendored_mpi_snapshot(resource)
    if (nrow(snapshot) == 0L) {
        excluded_resources <- c(excluded_resources, resource)
    } else {
        baselines[[resource]] <- snapshot
    }
}

all_rows <- c(
    list(metalinks_v2, v1$data),
    baselines
)
fig03_rows <- do.call(dplyr::bind_rows, all_rows)

resources <- sort(unique(fig03_rows$resource))
register_category_colours(
    'resources',
    setNames(
        palette_n(as.integer(length(resources)), unknown = FALSE),
        resources
    )
)

relation_types <- sort(unique(fig03_rows$relation_type))
register_category_colours(
    'interaction_types',
    setNames(
        palette_n(as.integer(length(relation_types)), unknown = FALSE),
        relation_types
    )
)

panels <- list(
    coverage = fig03_coverage_panel(fig03_rows, width_mm = 120L),
    metabolite_classes = fig03_metabolite_class_panel(fig03_rows, width_mm = 120L),
    protein_classes = fig03_protein_class_panel(fig03_rows, width_mm = 120L),
    evidence = fig03_evidence_confidence_panel(fig03_rows, width_mm = 120L),
    source_relationship = fig03_source_relationship_panel(fig03_rows, width_mm = 120L)
)

for (name in names(panels)) {
    ggsave(
        filename = file.path(out_dir, paste0(name, '.pdf')),
        plot = panels[[name]],
        width = 120,
        height = 90,
        units = 'mm'
    )
    ggsave(
        filename = file.path(out_dir, paste0(name, '.svg')),
        plot = panels[[name]],
        width = 120,
        height = 90,
        units = 'mm'
    )
}

composite <- compose_patchwork(
    panels,
    layout = list(ncol = 2)
)

ggsave(
    filename = file.path(out_dir, 'fig04-metalinks-versions.pdf'),
    plot = composite,
    width = 240,
    height = 220,
    units = 'mm'
)
ggsave(
    filename = file.path(out_dir, 'fig04-metalinks-versions.svg'),
    plot = composite,
    width = 240,
    height = 220,
    units = 'mm'
)

caption_info <- compose_caption(
    figure_id = 'fig04-metalinks-versions',
    composite_pdf = file.path(out_dir, 'fig04-metalinks-versions.pdf'),
    caption_source = 'figures/fig04-metalinks-versions/caption.tex',
    out_dir = out_dir,
    panel_count = 5L
)

queries <- list(query_record(metalinks_v2_rows))
write_sidecar(
    artifact_id = 'fig04-metalinks-versions',
    artifact_path = file.path(out_dir, 'fig04-metalinks-versions.pdf'),
    deployments = list(dep5$deployment),
    manifests = list(dep5$manifest),
    script_path = 'figures/fig04-metalinks-versions/build.R',
    queries = queries,
    external_inputs = c(list(v1$external_input), unname(lapply(names(baselines), function(resource) {
        list(
            kind = 'vendored-mpi-baseline',
            path = file.path('inst/extdata/fig04-mpi-baselines', paste0(snapshot_slug(resource), '.csv')),
            source = resource,
            fingerprint = NA_character_
        )
    }))),
    parameters = list(
        counting_basis = list(
            metabolites = 'HMDB',
            proteins = 'UniProt'
        ),
        included_resources = resources,
        excluded_resources = excluded_resources,
        optional_artifacts = character(0),
        class_system = list(
            metabolites = 'ChEBI-oriented baseline when available; current baseline snapshots may still carry source-native labels',
            proteins = 'UniProt-derived / Guide to Pharmacology'
        ),
        evidence_availability = lapply(resources, function(resource) {
            resource_rows <- fig03_rows[fig03_rows$resource == resource, ]
            list(
                resource = resource,
                source_count = TRUE,
                citation_count = any(!is.na(resource_rows$citation_count)),
                affinity_value = any(!is.na(resource_rows$affinity_value)),
                curation_mode = any(!is.na(resource_rows$curation_mode) & resource_rows$curation_mode != '')
            )
        })
    ),
    seed = pipeline_seed(),
    caption = caption_info
)

logger::log_info('Figure 4 build complete')
