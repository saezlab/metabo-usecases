# Cancer cell lines use case

Metabolomics of cancer cell lines from Shorthouse et al. 2022
(<https://pubmed.ncbi.nlm.nih.gov/36321551/>), who reanalysed the data of
Cherkaoui et al. 2022. Raw inputs and their provenance:
`data/raw/cancer-cell-lines/README.md`.

## Layout

| Path | Content |
|---|---|
| `run_analysis.R` | Entry point: renders the steps in fresh R sessions (`--list` shows them) |
| `config/config.yml` | Paths, excluded and merged lines, manual accessions, Cellosaurus release |
| `R/` | Helper functions sourced by the Rmds |
| `preprocessing/` | Refactored pre-processing (steps 01–03) and `verify_against_previous.R` |
| `docs/preprocessing_changes.md` | Every change to the previous pipeline, with its verified effect |
| `output/` | QC tables, plots and rendered reports (gitignored) |
| `01_preprocess.Rmd` … `04_differential_analysis.Rmd` | Previous pipeline; kept because Figure 6 reads its differential analysis until that is replaced |

## Running

From the repository root, with pandoc available (on beauty:
`nix-shell -p pandoc`):

```bash
Rscript analyses/cancer-cell-lines/run_analysis.R          # all steps
Rscript analyses/cancer-cell-lines/run_analysis.R 02 03    # selected steps
```

## Pre-processing steps

| Step | Reads | Writes to `data/derived/cancer-cell-lines/shorthouse_v1/` |
|---|---|---|
| `01_import_shorthouse.Rmd` | Datasets EV1–EV3 | `01_injections.rds` (injection × ion), `01_lines.rds` (line × ion, mean over injections) |
| `02_cell_lines_mutations.Rmd` | `01_lines.rds`, `cellosaurus.txt` | `02_cell_lines.tsv`, `02_sequence_variations.tsv`, `02_gene_status.tsv` |
| `03_feature_ids.Rmd` | `01_lines.rds`, MetaProViz prior knowledge | `03_features.rds` |

Each step also writes a `0x_manifest.json` with input checksums, git commit
and package versions.

## Previous pipeline

`01_preprocess.Rmd` → `02_add_metadata.Rmd` → `03_feature_processing.Rmd` →
`04_differential_analysis.Rmd`, writing to `data/derived/cancer-cell-lines/`.
`04` produces the workbook read by Figure 6. Known problems of the previous
steps are listed in `docs/preprocessing_changes.md`.
