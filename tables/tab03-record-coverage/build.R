# tables/tab03-record-coverage/build.R
#
# Orchestrates the FR-015a resource × record-type checkmark table.
# Queries the OmniPath Postgres (dev3 by default; the Structures row
# is routed to dev4 via the panel_deployment registry), pivots to
# wide, builds the LaTeX body with rotated column labels, saves PDF
# + CSV, composes the caption-and-table PDF + plain-text caption, and
# writes the provenance sidecar carrying both deployment build_ids.
#
# Sourced by rebuild.R; safe to source standalone too.

suppressPackageStartupMessages({
    library(metabo.figures)
})

setup_pipeline_log("build:tab03-record-coverage")
set.seed(pipeline_seed())

out_dir <- "tables/tab03-record-coverage/out"
fs::dir_create(out_dir)

# ---- Deployments + manifests --------------------------------------------

logger::log_info("Resolving dev3 + dev4 for tab03-record-coverage")
dep3 <- deployment_provenance("dev3")
dep4 <- deployment_provenance("dev4")

# ---- Data layer ---------------------------------------------------------

threshold <- 1L

logger::log_info(
    "Running per-row bitmap-intersection queries (threshold={threshold})"
)
data_long <- record_coverage_long()
queries   <- record_coverage_queries(data_long)
overrides <- attr(data_long, "expert_overrides")

resource_labels <- resources_label_map("tab03-record-coverage")
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
    slug         = "tab03-record-coverage",
    max_width_mm = table_width_mm
)

# ---- Caption (FR-040..FR-041a, SC-011) ----------------------------------

caption_info <- tables_compose_caption(
    table_id       = "tab03-record-coverage",
    table_pdf      = artifacts$pdf,
    caption_source = "tables/tab03-record-coverage/caption.tex",
    out_dir        = out_dir,
    panel_count    = 1L,
    body_width_mm  = table_width_mm
)

# ---- Provenance sidecar -------------------------------------------------

write_sidecar(
    artifact_id    = "tab03-record-coverage",
    artifact_path  = artifacts$pdf,
    deployments    = list(dep3$deployment, dep4$deployment),
    manifests      = list(dep3$manifest,   dep4$manifest),
    script_path    = "tables/tab03-record-coverage/build.R",
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
logger::log_info("tab03-record-coverage complete")
