#!/usr/bin/env Rscript
# Entry point of the cancer cell lines use case.
#
# Renders the steps in order, each in a fresh R session, so that every step
# depends only on the files written by earlier steps.
#
# Usage (from the repository root):
#   Rscript analyses/cancer-cell-lines/run_analysis.R            # all steps
#   Rscript analyses/cancer-cell-lines/run_analysis.R 02 03      # selected steps
#   Rscript analyses/cancer-cell-lines/run_analysis.R --list     # show steps
#
# Rendering needs pandoc on the PATH (on beauty: nix-shell -p pandoc).

steps <- c(
    "01" = "preprocessing/01_import_shorthouse.Rmd",
    "02" = "preprocessing/02_cell_lines_mutations.Rmd",
    "03" = "preprocessing/03_feature_ids.Rmd",
    "04" = "preprocessing/04_depmap_status.Rmd",
    "05" = "preprocessing/05_analysis_dataset.Rmd",
    # analysis steps (Phase 3 and 4)
    "m02" = "analysis/02_mutation_landscape.Rmd",
    "m03" = "analysis/03_contrast_selection.Rmd",
    "m04" = "analysis/04_metabolomics_qc.Rmd",
    "m05" = "analysis/05_positive_control_nrf2.Rmd",
    # audit notebooks: analyses that support decisions (Phase 1b)
    "a01" = "audit/01_scale_transform.Rmd",
    "a02" = "audit/02_structure_batch_tissue.Rmd",
    "a03" = "audit/03_tissue_and_overlaps.Rmd",
    "a04" = "audit/04_mutation_sources.Rmd",
    "a05" = "audit/05_metabolite_annotation.Rmd",
    "a06" = "audit/06_ortmayr_descriptive.Rmd",
    "a07" = "audit/07_functional_groups.Rmd",
    "a08" = "audit/08_genotype_reliability.Rmd"
)

use_case_dir <- here::here("analyses", "cancer-cell-lines")
source(file.path(use_case_dir, "R", "use_case.R"))
report_dir <- uc_dir(load_config(), "output", "reports")

args <- commandArgs(trailingOnly = TRUE)
if ("--list" %in% args) {
    cat(paste(names(steps), steps), sep = "\n")
    quit(status = 0)
}
selected <- if (length(args)) args else names(steps)
unknown <- setdiff(selected, names(steps))
if (length(unknown)) stop("Unknown step(s): ", paste(unknown, collapse = ", "))
if (!rmarkdown::pandoc_available()) stop("pandoc not found; on beauty run inside nix-shell -p pandoc")

for (id in selected) {
    rmd <- file.path(use_case_dir, steps[[id]])
    message(sprintf("[%s] %s: rendering %s", format(Sys.time(), "%H:%M:%S"), id, steps[[id]]))
    started <- Sys.time()
    callr::r(
        function(input, output_dir) {
            rmarkdown::render(input, output_dir = output_dir, envir = new.env(), quiet = TRUE)
        },
        args = list(input = rmd, output_dir = report_dir),
        show = FALSE
    )
    message(sprintf("[%s] %s: done in %.0f s", format(Sys.time(), "%H:%M:%S"), id,
                    as.numeric(difftime(Sys.time(), started, units = "secs"))))
}
