# cancer-cell-lines

Metabolomics of cancer cell lines from Shorthouse et al. 2022
(<https://pubmed.ncbi.nlm.nih.gov/36321551/>, data of Cherkaoui et al.
2022): pre-processing, mutation status, mutation-based contrasts,
differential analysis and enrichment on MetSigDB sets.

**Status (2026-11-24):** metabolomics phase complete (Checkpoint 3,
`docs/checkpoint_3.md`); the multi-omics expansion is planned in
`docs/handoff_multiomics.md`. The figure (`figures/cancer-cell-lines/`) is
not built yet. OmniPath Metabo is read from dev5 Postgres for now
(`R/metsigdb.R`); to be moved to the Parquet release and its clients
(`docs/DB_ACCESS.md`).

## Layout

| Path | Content |
|---|---|
| `run_analysis.R` | Entry point: renders the steps in fresh R sessions (`--list` shows them) |
| `config/config.yml` | Paths, excluded and merged lines, contrasts, manual accessions, releases |
| `R/` | Helper functions sourced by the Rmds |
| `preprocessing/` | Steps 01–06: from the raw data to the analysis-ready data set |
| `analysis/` | Steps m02–m08: mutation landscape, contrasts, QC, differential analysis, enrichment |
| `audit/` | a01–a13: notebooks that support individual decisions |
| `docs/` | Use case overview, analysis plan, decision log, caveats, checkpoints, open questions |
| `output/` | Rendered reports, QC tables and plots (gitignored) |

Start with `docs/use_case_overview.md`; decisions are numbered (D1, D2, …)
in `docs/decision_log.md`.

## Running

From the repository root, with pandoc available (on beauty:
`nix-shell -p pandoc`):

```bash
Rscript analyses/cancer-cell-lines/run_analysis.R            # all steps
Rscript analyses/cancer-cell-lines/run_analysis.R 01 02      # selected steps
Rscript analyses/cancer-cell-lines/run_analysis.R --list     # list the steps
```

## Pre-processing

Raw inputs and their sources: `data/raw/cancer-cell-lines/README.md`.
Results go to `data/derived/cancer-cell-lines/shorthouse_v1/`; each step
also writes a `0x_manifest.json` with input checksums, git commit and
package versions.

| Step | Writes |
|---|---|
| `01_import_shorthouse.Rmd` | `01_injections.rds` (injection × ion), `01_lines.rds` (line × ion) |
| `02_cell_lines_mutations.Rmd` | `02_cell_lines.tsv`, `02_sequence_variations.tsv`, `02_gene_status.tsv` (Cellosaurus) |
| `03_feature_ids.Rmd` | `03_features.rds` (metabolite identifiers, MetaProViz) |
| `04_depmap_status.Rmd` | `04_mutation_status.tsv`, `04_copy_number.rds`, `04_homozygous_deletions.tsv`, `04_depmap_models.tsv` (DepMap 24Q4) |
| `05_analysis_dataset.Rmd` | `shorthouse_analysis_v1.rds`, the input of `analysis/` |
| `06_metsigdb_sets.Rmd` | `06_sets.tsv`, `06_ion_sets.tsv` (ions linked to MetSigDB sets, from OmniPath Metabo) |

`preprocessing/verify_against_previous.R` compares the results with the
previous pipeline, read from git history (commit `03c9674`, see
`config.yml`); `docs/preprocessing_changes.md` lists every change.
