# tables/tab01-id-resolving/build.R
#
# Orchestrates the FR-014 Methods table: queries the OmniPath Postgres
# (dev3 via the panel_deployment registry), pivots to wide, renders via
# gt + the LaTeX backend, saves PDF + CSV, composes the caption-and-
# table PDF + plain-text caption, and writes the provenance sidecar.
#
# Sourced by rebuild.R; safe to source standalone too.

suppressPackageStartupMessages({
    library(metabo.figures)
})

setup_pipeline_log("build:tab01-id-resolving")
set.seed(pipeline_seed())

out_dir <- "tables/tab01-id-resolving/out"
fs::dir_create(out_dir)

# ---- Deployment + manifest -----------------------------------------------

logger::log_info("Resolving dev3 deployment for tab01-id-resolving")
dep3 <- deployment_provenance("dev3")

# ---- Data layer ----------------------------------------------------------

logger::log_info(
    "Querying identifier counts (15 hand-picked buckets); ",
    "expect ~30 minutes until identifier_source_count lands"
)
data_long <- tbl_id_resolving_counts()
resource_labels <- resources_label_map("tab01-id-resolving")
data_wide <- tbl_id_resolving_wide(
    data_long, resource_labels = resource_labels
)

queries <- list(query_record(data_long))

# ---- Render + save -------------------------------------------------------

gt_obj <- tbl_id_resolving_gt(data_wide)

artifacts <- tables_save_pdf_csv(
    gt_obj  = gt_obj,
    data    = data_wide,
    out_dir = out_dir,
    slug    = "tab01-id-resolving"
)

# ---- Caption (FR-040..FR-041a, SC-011) -----------------------------------

caption_info <- tables_compose_caption(
    table_id       = "tab01-id-resolving",
    table_pdf      = artifacts$pdf,
    caption_source = "tables/tab01-id-resolving/caption.tex",
    out_dir        = out_dir,
    panel_count    = 1L
)

# ---- Provenance sidecar --------------------------------------------------

write_sidecar(
    artifact_id    = "tab01-id-resolving",
    artifact_path  = artifacts$pdf,
    deployments    = list(dep3$deployment),
    manifests      = list(dep3$manifest),
    script_path    = "tables/tab01-id-resolving/build.R",
    queries        = queries,
    parameters     = list(buckets = vapply(
        tbl_id_resolving_buckets(),
        function(b) b$label,
        character(1L)
    )),
    seed           = pipeline_seed(),
    caption        = list(
        source_path           = caption_info$caption_source,
        source_kind           = if (endsWith(caption_info$caption_source, ".md")) "md" else "tex",
        panel_letter_count    = caption_info$panel_letter_count,
        composite_panel_count = caption_info$composite_panel_count
    )
)

pg_close_all_panel()
logger::log_info("tab01-id-resolving complete")
