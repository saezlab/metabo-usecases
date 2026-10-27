# Cancer cell lines use case

Metabolomics of cancer cell lines from Shorthouse et al. 2022
(<https://pubmed.ncbi.nlm.nih.gov/36321551/>). The current contrasts are
KRAS- and EGFR-mutant vs. other lung-cancer cell lines; further contrasts
within the cancer cell lines may follow. The differential-analysis output
feeds Figure 6.

## Run order

| Step | Reads | Writes to `data/derived/cancer-cell-lines/` |
|---|---|---|
| `01_preprocess.Rmd` | Shorthouse supplementary xlsx (`data/raw/cancer-cell-lines/`) | `Shorthouse_2022_preprocessed.RData` |
| `02_add_metadata.Rmd` | + `cellosaurus.txt` (`data/raw/cancer-cell-lines/`) | `Shorthouse_2022.RData`, `cellosaurus_entry_lookup.RData` |
| `03_feature_processing.Rmd` | `Shorthouse_2022.RData` | `Shorthouse_2022_extended.RData` |
| `04_differential_analysis.Rmd` | `Shorthouse_2022_extended.RData` | `ExtendedDataTable_DifferentialAnalysis_Shorthouse_LungMutations.xlsx` |

QC tables and plots of each step go to `output/<step>/` (gitignored).

## Paths

All paths are resolved from the repository root with `here::here()`, so the
Rmds can be knitted from any working directory. Do not put an `.Rproj` or
`.here` file inside `analyses/`, since `here` would then treat that folder as
the root.

## Known issue

`03_feature_processing.Rmd` loads
`features_Shorthouse_compatibility_check.RData` and
`features_Shorthouse_traverse_ids.RData`, but the code that computes them is
commented out. A fresh run therefore depends on these cached files.
