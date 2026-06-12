# tables/tab02-ramp-comparison/build.R
#
# Orchestrates the FR-015 Methods table: queries the OmniPath Postgres
# (dev4 via the panel_deployment registry — metabo_ramp_inchikey_conflict
# is dev4-only), renders via gt + the LaTeX backend, saves PDF + CSV,
# composes the caption-and-table PDF + plain-text caption, and writes
# the provenance sidecar.
#
# Sourced by rebuild.R; safe to source standalone too.

suppressPackageStartupMessages({
    library(metabo.figures)
})

setup_pipeline_log("build:tab02-ramp-comparison")
set.seed(pipeline_seed())

out_dir <- "tables/tab02-ramp-comparison/out"
fs::dir_create(out_dir)

# ---- Deployment + manifest -----------------------------------------------

logger::log_info(
    "Resolving dev4 deployment for tab02-ramp-comparison"
)
dep4 <- deployment_provenance("dev4")

# ---- Data layer ----------------------------------------------------------

example_count <- 2L

logger::log_info(
    "Querying RaMP conflict summary (example_count={example_count})"
)
data_summary <- tbl_ramp_comparison_summary(
    example_count = example_count
)

queries <- list(query_record(data_summary))

# ---- Render + save -------------------------------------------------------

gt_obj <- tbl_ramp_comparison_gt(data_summary)

table_width_mm <- 180L

artifacts <- tables_save_pdf_csv(
    gt_obj       = gt_obj,
    data         = data_summary,
    out_dir      = out_dir,
    slug         = "tab02-ramp-comparison",
    max_width_mm = table_width_mm
)

# ---- Caption (FR-040..FR-041a, SC-011) -----------------------------------

caption_info <- tables_compose_caption(
    table_id       = "tab02-ramp-comparison",
    table_pdf      = artifacts$pdf,
    caption_source = "tables/tab02-ramp-comparison/caption.tex",
    out_dir        = out_dir,
    panel_count    = 1L,
    body_width_mm  = table_width_mm
)

# ---- Provenance sidecar --------------------------------------------------

write_sidecar(
    artifact_id    = "tab02-ramp-comparison",
    artifact_path  = artifacts$pdf,
    deployments    = list(dep4$deployment),
    manifests      = list(dep4$manifest),
    script_path    = "tables/tab02-ramp-comparison/build.R",
    queries        = queries,
    parameters     = list(example_count = example_count),
    seed           = pipeline_seed(),
    caption        = list(
        source_path           = caption_info$caption_source,
        source_kind           = if (endsWith(caption_info$caption_source, ".md")) "md" else "tex",
        panel_letter_count    = caption_info$panel_letter_count,
        composite_panel_count = caption_info$composite_panel_count
    )
)

pg_close_all_panel()
logger::log_info("tab02-ramp-comparison complete")
