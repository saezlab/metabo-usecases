# figures/fig01-architecture/build.R
#
# Orchestrates the new (post-2026-06-14 six-figure renumbering) Figure 1:
# the workflow-architecture diagram (vendored Inkscape PDF/SVG, FR-005)
# plus the FR-043 statistics digest. No quantitative panels — those moved
# to figures/fig02-overview/. This script:
#
# - verifies the vendored architecture asset's SHA-256 against
#   inst/extdata/manual/architecture/README.md (FR-005a),
# - copies the vendored PDF + SVG to out/fig01-architecture.{pdf,svg},
# - invokes build_panel_a_digest() to refresh
#   panel-a-stats/{stats.json, stats.csv, stats.md, stats.pdf,
#   stats.provenance.json} (FR-043),
# - compiles the Figure 1 caption-and-figure PDF (FR-041) + caption.txt
#   (FR-041a) from caption.tex,
# - writes the provenance sidecar carrying the dev5 build_id, the
#   architecture asset fingerprint, and a pointer at the digest sidecar.
#
# Sourced by rebuild.R; safe to source standalone too.

suppressPackageStartupMessages({
    library(metabo.figures)
})

setup_pipeline_log("build:fig01-architecture")
set.seed(pipeline_seed())

out_dir <- "figures/fig01-architecture/out"
fs::dir_create(out_dir)

# ---- Deployment + manifests -----------------------------------------------

logger::log_info("Resolving dev5 deployment")
dep5 <- deployment_provenance("dev5")

# ---- Vendored architecture asset (FR-005, FR-005a) ------------------------

architecture_dir    <- "inst/extdata/manual/architecture"
architecture_pdf    <- file.path(
    architecture_dir, "omnipath-architecture-new2026.pdf"
)
architecture_svg    <- file.path(
    architecture_dir, "omnipath-architecture-new2026.svg"
)
architecture_readme <- file.path(architecture_dir, "README.md")

if (!file.exists(architecture_pdf)) {
    rlang::abort(sprintf(
        "Vendored architecture asset missing: %s. Re-vendor per %s.",
        architecture_pdf, architecture_readme
    ))
}

readme_pin <- stringr::str_match(
    paste(readLines(architecture_readme), collapse = "\n"),
    "pdf_sha256.*?`([0-9a-f]{64})`"
)[1L, 2L]
if (is.na(readme_pin)) {
    rlang::abort(sprintf(
        "Could not parse pdf_sha256 pin from %s", architecture_readme
    ))
}

architecture_sha256 <- digest::digest(
    file = architecture_pdf, algo = "sha256"
)
if (!identical(architecture_sha256, readme_pin)) {
    rlang::abort(sprintf(paste0(
        "Architecture asset SHA-256 mismatch (FR-005a). ",
        "Expected (README pin): %s\n",
        "Got      (on disk):    %s\n",
        "Manual asset changed — update %s and re-record fingerprint."
    ), readme_pin, architecture_sha256, architecture_readme))
}

file.copy(
    architecture_pdf,
    file.path(out_dir, "fig01-architecture.pdf"),
    overwrite = TRUE
)
if (file.exists(architecture_svg)) {
    file.copy(
        architecture_svg,
        file.path(out_dir, "fig01-architecture.svg"),
        overwrite = TRUE
    )
}
logger::log_info(paste0(
    "Figure 1 architecture vendored from {architecture_pdf} ",
    "(sha256={substr(architecture_sha256, 1L, 12L)})"
))

# ---- Panel A statistics digest (FR-043 family) ----------------------------

logger::log_info("Building Panel A statistics digest (FR-043)")
digest_result <- build_panel_a_digest(
    snapshot_id = snapshot_id(dep5$manifest),
    out_dir     = "figures/fig01-architecture/panel-a-stats",
    caption_sty = "tex/caption.sty"
)

# ---- Caption (FR-040..FR-041a, SC-011) ------------------------------------
#
# Figure 1 is a single-panel figure (the architecture diagram itself). The
# caption has no `(a)/(b)/...` sub-letters; FR-041b's count match is the
# trivial 1 == 1.

caption_info <- compose_caption(
    figure_id      = "fig01-architecture",
    composite_pdf  = file.path(out_dir, "fig01-architecture.pdf"),
    caption_source = "figures/fig01-architecture/caption.tex",
    out_dir        = out_dir,
    panel_count    = 1L
)

# ---- Provenance sidecar ---------------------------------------------------

write_sidecar(
    artifact_id    = "fig01-architecture",
    artifact_path  = file.path(out_dir, "fig01-architecture.pdf"),
    deployments    = list(dep5$deployment),
    manifests      = list(dep5$manifest),
    script_path    = "figures/fig01-architecture/build.R",
    queries        = list(),
    external_inputs = list(
        list(
            kind        = "architecture-asset",
            path        = architecture_pdf,
            source      = paste0(
                "Inkscape source on the contributor's machine; ",
                "see ", architecture_readme
            ),
            fingerprint = architecture_sha256
        ),
        list(
            kind        = "panel-a-stats-digest",
            path        = digest_result$sidecar,
            source      = paste0(
                "Pipeline-generated FR-043 digest; ",
                "config sha256=", digest_result$config_sha256
            ),
            fingerprint = digest_result$snapshot_id
        )
    ),
    parameters     = list(),
    seed           = pipeline_seed(),
    caption        = list(
        source_path           = caption_info$caption_source,
        source_kind           = if (endsWith(
            caption_info$caption_source, ".md"
        )) "md" else "tex",
        panel_letter_count    = caption_info$panel_letter_count,
        composite_panel_count = caption_info$composite_panel_count
    )
)

pg_close_all_panel()
logger::log_info("fig01-architecture complete")
