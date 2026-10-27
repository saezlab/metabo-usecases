# tables/record-coverage/build.R
#
# Orchestrates the FR-015a resource × record-type checkmark table.
# Queries the OmniPath Postgres (dev5 integrated build — post-
# 2026-06-14 every row, including the Structures row, resolves on
# dev5), pivots to wide, builds the LaTeX body with rotated column
# labels, saves PDF + CSV, composes the caption-and-table PDF +
# plain-text caption, and writes the provenance sidecar.
#
# Sourced by rebuild.R; safe to source standalone too.

suppressPackageStartupMessages({
    library(metabo.figures)
})

setup_pipeline_log("build:record-coverage")
set.seed(pipeline_seed())

out_dir <- "tables/record-coverage/out"
fs::dir_create(out_dir)

# ---- Deployments + manifests --------------------------------------------

logger::log_info("Resolving dev5 for record-coverage")
dep5 <- deployment_provenance("dev5")

# ---- Data layer ---------------------------------------------------------

threshold <- 1L

logger::log_info(
    "Running per-row bitmap-intersection queries (threshold={threshold})"
)
data_long <- record_coverage_long()
queries   <- record_coverage_queries(data_long)
overrides <- attr(data_long, "expert_overrides")

resource_labels <- resources_label_map("record-coverage")
data_wide <- record_coverage_wide(
    data_long,
    threshold       = threshold,
    resource_labels = resource_labels
)

# ---- Render + save ------------------------------------------------------

tex_body <- record_coverage_latex(data_wide)

table_width_mm <- 360L

artifacts <- tables_save_latex_pdf_csv(
    tex_body     = tex_body,
    data         = data_wide,
    out_dir      = out_dir,
    slug         = "record-coverage",
    max_width_mm = table_width_mm
)

# ---- Caption (FR-040..FR-041a, SC-011) ----------------------------------

caption_info <- tables_compose_caption(
    table_id       = "record-coverage",
    table_pdf      = artifacts$pdf,
    caption_source = "tables/record-coverage/caption.tex",
    out_dir        = out_dir,
    panel_count    = 1L,
    body_width_mm  = table_width_mm
)

# ---- Provenance sidecar -------------------------------------------------

write_sidecar(
    artifact_id    = "record-coverage",
    artifact_path  = artifacts$pdf,
    deployments    = list(dep5$deployment),
    manifests      = list(dep5$manifest),
    script_path    = "tables/record-coverage/build.R",
    queries        = queries,
    parameters     = list(
        threshold        = threshold,
        expert_overrides = overrides
    ),
    seed           = pipeline_seed(),
    caption        = list(
        source_path           = caption_info$caption_source,
        source_kind           = if (endsWith(caption_info$caption_source, ".md")) "md" else "tex",
        panel_letter_count    = caption_info$panel_letter_count,
        composite_panel_count = caption_info$composite_panel_count
    )
)

pg_close_all_panel()
logger::log_info("record-coverage complete")
