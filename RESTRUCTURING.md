# Restructuring after the first preprint

This document records how the repository was reorganised after the first
OmniPath Metabo preprint, on the `restructure` branch. It lists what moved or
was renamed and why, pins each change to its commit, and describes the state
to continue from (e.g. with a new Spec Kit specification).

- **Baseline:** tag [`preprint_v0`](https://github.com/saezlab/metabo-usecases/tree/preprint_v0)
  = commit `25136b9`, the state used for the first preprint. Anything removed
  below can be recovered from there, e.g.
  `git checkout preprint_v0 -- backup/`.
- **Branch:** `restructure`, short-lived. It is merged into `main` after review
  by the collaborators; development then continues trunk-based on `main`.
- **Status:** all changes are syntax-checked only. Nothing has been run yet,
  see [Before merging](#before-merging).

## Goals

1. The repository should read at a glance: separate top-level folders for
   analyses, data, figures, tables and docs.
2. Names describe content. Figure and table numbers change with every
   manuscript revision, so they live in one file (`manuscript.yaml`) and
   nowhere else.
3. Every file in `data/derived/` is written by a script in the repository.
   Outputs whose origin is unknown are not trusted.
4. Everything must stay reproducible for publication (environment work
   follows after the restructure, see [Planned work](#planned-work)).

## Layout after the restructure

```
data/raw/<name>/        external inputs, README with the source
      │
analyses/<name>/        use-case analyses (Rmd, Python, notebooks)
      │
data/derived/<name>/    analysis results, written only by analyses/
      │                                   OmniPath Postgres (live)
figures/<id>/, tables/<id>/  build.R ◄────────────┘   data/vendored/ (pinned snapshots)
      │
figures/<id>/out/, tables/<id>/out/, build/   (all gitignored)
```

| Path | Contents |
|---|---|
| `analyses/` | `cancer-cell-lines/`, `cancer-cell-lines-cosmos/`, `azelate/`, `ramp-ambiguity/` |
| `data/raw/`, `data/derived/`, `data/vendored/` | External inputs, analysis results, pinned third-party snapshots |
| `figures/`, `tables/` | One folder per artifact, named by content; `manual/` for hand-made assets |
| `manuscript.yaml` | Figure and table numbering; the only place numbers are defined |
| `R/`, `man/`, `tests/`, `DESCRIPTION`, `NAMESPACE` | R package `metabo.figures` (unchanged in structure) |
| `inst/extdata/` | Package resources only: palettes, logos, JSON schema, connection template |
| `python/`, `tex/`, `lib/` | Pipeline helpers |
| `docs/` | `CONFIGURATION.md`, `DB_ACCESS.md` |
| `build/` | Pipeline run output: logs, manifests, manuscript bundle (gitignored) |
| `AGENTS.md`, `README.md` | Conventions for agents; overview for humans |

## Changes by commit

In order on `restructure`, starting from `preprint_v0` (`25136b9`):

| Commit | Change | Why |
|---|---|---|
| `461c6b7` | Add `.gitattributes` (LF line endings, binary data types) | Same checkouts on Windows/WSL, macOS and NixOS. No file had CRLF, so no content changed. |
| `746cdb3` | Remove `backup/`, `prompts-archive/`, D3 dashboard in `figures/` | Stale copy of old layout; LLM prompts; standalone dashboard unrelated to the pipeline whose README shadowed `figures/`. |
| `0219cb2` | Move the cancer cell lines analysis to `analyses/` + `data/{raw,derived}/` | Separate analysis code, external inputs and results; the Figure 6 input no longer lives in a gitignored `Results/` folder. |
| `25a13ef` | Rmds resolve paths from the repo root (`here::here()`) | `../data` paths broke after the move and only worked from one working directory. |
| `5658aac` | `lung-cancer` → `cancer-cell-lines` | The use case may get other contrasts within cancer cell lines. |
| `c66b7f4` | RaMP notebook → `analyses/ramp-ambiguity/` | An analysis like the others; one top-level folder fewer. |
| `56ee69e` | Remove `cosmos-figure/COSMOS_schematic_figure_V2-CS_YB.svg` | Not referenced anywhere. |
| `7d80c04` | `CONFIGURATION.md`, `DB_ACCESS.md` → `docs/`; `omnipath_metabo_case1/agents.md` → root `AGENTS.md` | Docs in one place; agent notes rewritten for the new layout, credentials not duplicated. |
| `c3fceb6` | COSMOS-PKN extension → `analyses/cancer-cell-lines-cosmos/` | Next to the analysis it builds on; hard-coded `/Users/priscilla/...` paths replaced; untrusted copy of the DA xlsx removed. |
| `069aef7` | Azelate workflow → `analyses/azelate/`, results → `data/derived/azelate/` | Same layout as the other analyses; `fig3_panel_*` names were an outdated figure number. |
| `ff17b6e` | Rewrite `README.md` | The old README had outdated numbering and relied on a private spec repository. |
| `247ea7c` | `.Rbuildignore`: `analyses/`, `data/`, `docs/`, … | `data/` has a special meaning in R packages (bundled datasets); `R CMD check` would have treated the files there as package data. |
| `3611cde` | Split `inst/extdata/` by purpose | Pinned data belongs under `data/`; the loaders mixed `system.file()` with relative fallbacks, which left a dead path behind. |
| `80128de` | Figure/table folders named by content; `manuscript.yaml` | The June renumbering left stale ids everywhere; the next one would have done the same. |
| `bf64b32` | Functions, R files and output files named by content | Names carried old figure numbers (`fig03_*` drew Figure 4) or private spec ids (`fr007a`). |
| `f6b6778` | `logs/`, `manifests/`, bundle → `build/` | One gitignored place for pipeline output instead of three root folders. |

## Rename maps

### Figures and tables (`80128de`)

| Old folder / id | New folder / id | Number (`manuscript.yaml`) |
|---|---|---|
| `figures/fig01-architecture` | `figures/architecture` | Figure 1 |
| `figures/fig02-overview` | `figures/database-content` | Figure 2 |
| — | — (webapp screenshots, manual, not built here) | Figure 3 |
| `figures/fig04-metalinks-versions` | `figures/metalinks-versions` | Figure 4 |
| `figures/fig05-cosmos-pkn` | `figures/cosmos-pkn` | Figure 5 |
| `figures/fig06-lungcancer-usecase` | `figures/cancer-cell-lines` | Figure 6 |
| `tables/tab01-id-resolving` | `tables/id-resolving` | Table 1 |
| `tables/tab02-ramp-comparison` | `tables/ramp-comparison` | Table 2 |
| `tables/tab03-record-coverage` | `tables/record-coverage` | Table 3 |

Ids left over from before the June renumbering were mapped to their current
folder: `fig01-overview` → `database-content` (or `architecture` where it
meant the architecture asset), `fig03-metalinks-versions` →
`metalinks-versions`, `fig04-cosmos-pkn` → `cosmos-pkn`. Artifact ids, output
files (`fig05-cosmos-pkn.pdf` → `cosmos-pkn.pdf`, …), `panel_id` defaults and
log labels (`[fig06]` → `[cancer-cell-lines]`) follow the folder names.
`rebuild.R` now orders builds and the manuscript bundle by `manuscript.yaml`
instead of by the number prefix.

Test files: `test-fig01-determinism.R` → `test-database-content-determinism.R`,
`test-fig01-provenance.R` → `test-database-content-provenance.R`,
`test-fig03-comparison-basis.R` → `test-comparison-basis.R`,
`test-fig03-fig04-palette-consistency.R` → `test-palette-consistency.R`,
`test-fig05-composite-merge.R` → `test-cosmos-pkn-composite-merge.R`.

### Functions (`bf64b32`)

| Old | New |
|---|---|
| `fr007a_*` (`fr007a_entities`, `_interactions`, `_associations`, `_identifiers_major`, `_identifiers_authoritative`, `_structures`, `_literature`, `_overview`, …) | `resource_overview_*` |
| `plot_fr007a_overview`, `plot_fr007a_total` | `plot_resource_overview`, `plot_resource_overview_total` |
| `fr007b_coverage`, `plot_fr007b_coverage` | `coverage_profile`, `plot_coverage_profile` |
| `fr007c_overlap`, `fr007c_node_sizes`, `plot_fr007c_networks` | `resource_overlap`, `resource_overlap_node_sizes`, `plot_resource_overlap` |
| `fr007d_entity_x_interaction`, `plot_fr007d_matrix` | `entity_by_interaction_type`, `plot_entity_interaction_matrix` |
| `fr007e_specificity_by_category`, `plot_fr007e` | `specificity_by_category`, `plot_specificity_by_category` |
| `fig03_coverage_panel`, `_metabolite_class_panel`, `_protein_class_panel`, `_metalinks_overview_panel`, `_relationship_types_panel` | `metalinks_coverage_panel`, `metalinks_metabolite_class_panel`, `metalinks_protein_class_panel`, `metalinks_overview_panel`, `metalinks_relationship_types_panel` |
| `empty_fig03_panel`, `fig03_baseline_dir` | `empty_metalinks_panel`, `mpi_baseline_dir` |
| `fig04_cosmos_comparison_panel`, `fig04_compartment_panel`, `fig04_resource_contribution_panel`, `fig04_metalinks_cosmos_panel` | `cosmos_comparison_panel`, `cosmos_compartment_panel`, `cosmos_resource_contribution_panel`, `metalinks_cosmos_panel` |

R files: `R/{data,plots}-fr007a.R` → `-resource_overview.R`, `fr007b` →
`coverage_profile`, `fr007c` → `resource_overlap`, `fr007d` →
`entity_interaction`, `fr007e` → `specificity`.

Output files of `figures/database-content/`: `fr007a-total` →
`resource-overview-total`, `fr007a-overview-supplementary` →
`resource-overview-supplementary`, `fr007b-coverage` → `coverage-profile`,
`fr007c-networks` → `resource-overlap-networks`, `fr007d-matrix` →
`entity-interaction-matrix`, `fr007e-specificity` → `specificity-by-category`.

`NAMESPACE` and `man/` were updated by hand with the same mapping (roxygen2
was not available), keeping roxygen's sort order.

### Analyses and data

| Old | New | Commit |
|---|---|---|
| `omnipath_metabo_case1/scripts/Preprocesss_Shorthouse_2022_data.Rmd` | `analyses/cancer-cell-lines/01_preprocess.Rmd` | `0219cb2`, `5658aac` |
| `…/scripts/Metadata_addition.Rmd` | `analyses/cancer-cell-lines/02_add_metadata.Rmd` | 〃 |
| `…/scripts/Feature_processing.Rmd` | `analyses/cancer-cell-lines/03_feature_processing.Rmd` | 〃 |
| `…/scripts/Differential_Analysis.Rmd` | `analyses/cancer-cell-lines/04_differential_analysis.Rmd` | 〃 |
| `…/data/msb202211006-sup-000{2,3,4}-datasetev{1,2,3}.xlsx`, `…/data/README.md` | `data/raw/cancer-cell-lines/` | 〃 |
| `…/data/*.RData` (Shorthouse intermediates, Cellosaurus lookup) | `data/derived/cancer-cell-lines/` | 〃 |
| `…/Results/Differential_Analysis/ExtendedDataTable_DifferentialAnalysis_Shorthouse_LungMutations.xlsx` | `data/derived/cancer-cell-lines/` (same name) | 〃 |
| `…/scripts/01_cosmos_pkn.py`, `02_connect_dem_pkn.R`, `03_visualization.R` | `analyses/cancer-cell-lines-cosmos/` | `c3fceb6` |
| `…/session_2026-05-20_DEM_COSMOS_PKN.md` | `analyses/cancer-cell-lines-cosmos/NOTES_2026-05-20.md` | 〃 |
| `…/data/pkn_allosteric.csv`, `pkn_enzyme_metabolite.csv` | `data/derived/cancer-cell-lines-cosmos/` | 〃 |
| `…/azelate_fig_tables/notebooks/*`, `README.md` | `analyses/azelate/` | `069aef7` |
| `…/azelate_fig_tables/outputs/{csvs/,azelaic_acid_tables.xlsx,manifest.json}` | `data/derived/azelate/` | 〃 |
| `fig3_panel_C_data` (table, CSV, sheet) | `interaction_types_by_source` | 〃 |
| `fig3_panel_D_data` (table, CSV, sheet) | `cancer_assoc_by_sample_type` (Excel sheet names are limited to 31 characters) | 〃 |
| `RaMP_Exploration/2026-06-ramp-ambiguity-investigation.ipynb` | `analyses/ramp-ambiguity/ramp_ambiguity_investigation.ipynb` | `c66b7f4` |
| `inst/extdata/cosmos/` | `data/vendored/cosmos/` | `3611cde` |
| `inst/extdata/fig04-mpi-baselines/` | `data/vendored/mpi-baselines/` | 〃 |
| `inst/extdata/manual/architecture/` | `figures/architecture/manual/` | 〃 |
| `CONFIGURATION.md`, `DB_ACCESS.md` | `docs/` | `7d80c04` |
| `omnipath_metabo_case1/agents.md` | `AGENTS.md` (rewritten) | 〃 |
| `logs/`, `manifests/`, `out/manuscript-bundle.pdf` | `build/logs/`, `build/manifests/`, `build/manuscript-bundle.pdf` | `f6b6778` |

`omnipath_metabo_case1/` has no tracked files left. Folder-specific QC
output of each analysis now goes to `analyses/<name>/output/` (gitignored).

### Removed

| What | Commit | Reason |
|---|---|---|
| `backup/01/` | `746cdb3` | Stale copy of an earlier `figures/` + `tables/` layout |
| `prompts-archive/` | `746cdb3` | LLM prompts, not project content |
| D3 dashboard in `figures/` (`index.html`, `package.json`, `src/`, `scripts/`, `sql/`, `data/`, `.env.example`, `README.md`) | `746cdb3` | Unrelated to the pipeline, referenced only by itself |
| `omnipath_metabo_case1/data/CCLE_proteomics_metadata_extended.csv` | `0219cb2` | Its generating script had been deleted in `9a1f8d6` |
| `omnipath_metabo_case1/data/ExtendedDataTable_…xlsx` | `c3fceb6` | Manual copy that differed from the output of `04_differential_analysis.Rmd` |
| Azelate preview PNGs + figure index | `069aef7` | Regenerated into `analyses/azelate/output/`; Figure 6 draws those panels in R |
| `cosmos-figure/` | `56ee69e` | Not referenced |
| tracked `.DS_Store`, `.gitkeep` in `logs/`, `manifests/`, `inst/extdata/{cosmos,manual}/` | various | OS noise; folders are created on demand |

## Conventions

These also appear in `AGENTS.md` and `README.md`.

- Analyses read from `data/raw/` (and earlier `data/derived/` steps) and write
  to `data/derived/<name>/`. Figures and tables read only from `data/derived/`,
  `data/vendored/` and the live database.
- Every file in `data/derived/` is written by a script in `analyses/`.
- No personal absolute paths. R uses `here::here()` in analyses and paths
  relative to the repository root in the pipeline; Python uses
  `Path(__file__)`.
- Folders, ids, functions and file names are named by content, never by
  figure or table number. Numbers live only in `manuscript.yaml`.
- Paper deliverables come out of `figures/` or `tables/`; `data/derived/`
  holds intermediate results.

## Before merging

None of the commits above was executed, only syntax-checked. The local
machine used for the restructure has R 4.1.2 (the pipeline needs R ≥ 4.3) and
no access to the databases.

- [ ] Full pipeline on beauty: `./rebuild.sh --dry-run`, lint,
      `devtools::test()`, then `./rebuild.sh --check --bundle`. Confirm every
      figure and table builds and the bundle is in `manuscript.yaml` order.
- [ ] `test-palette-consistency.R` pointed at build paths that no longer
      existed, so its build-script checks were skipped silently since June.
      They run again now and may fail.
- [ ] Regenerate `man/` and `NAMESPACE` with `devtools::document()`; expect
      no or only minor changes.
- [ ] Re-run `analyses/cancer-cell-lines/01`→`04` with the new paths.
- [ ] Collaborators review their moved parts (open questions are in each
      README): Yunfan `analyses/cancer-cell-lines-cosmos/`, Daniele
      `analyses/azelate/`, Christina/Jonathan `analyses/ramp-ambiguity/`,
      Denes the pipeline changes.
- [ ] After the re-run: delete the local, untracked `omnipath_metabo_case1/`
      and its entries in `.gitignore`.

## Known issues (not fixed in the restructure)

- **Old COSMOS PKN is not pinned.** `cosmos_old_pkn()` installs cosmosR from
  GitHub HEAD at build time and uses its `meta_network`. It ignores the
  vendored `data/vendored/cosmos/meta_network.RData` (cosmosR 1.18.1), although
  its documentation says it reads it. `test-comparison-basis.R` expects
  attributes the function does not set. Fixing it may change Figure 5, so this
  needs Denes.
- **Cached intermediates in `03_feature_processing.Rmd`.** It loads
  `features_Shorthouse_compatibility_check.RData` and
  `features_Shorthouse_traverse_ids.RData`, but the code computing them is
  commented out, so a fresh run depends on the cached files.
- **Azelate workbook out of date.** The tracked `azelaic_acid_tables.xlsx` and
  `manifest.json` still use the `fig3_panel_*` names until the next run.
- **Dead LaTeX template.** `tex/compose_fig01.tex` is not used by any code;
  `figures/database-content/` mentions `compose_fig01.tex` / `compose_fig02.tex`
  in comments.
- **Private spec references.** Comments and READMEs still cite `FR-0xx` and
  `T0xx` ids from the saezverse spec; see the mapping below.

## Planned work

- **Environment (next step).** A Docker image as the reference environment:
  `rocker/r-ver` with a dated Posit Package Manager CRAN snapshot, a fixed
  Bioconductor release, GitHub packages at fixed commits, Python via
  `uv.lock`, plus LaTeX and pdftk. The image acts as the lock file. Results
  for the paper are produced in the container (wrapper script or Dev
  Container); quick native runs stay possible. CI builds the image. No renv
  during development; an `renv.lock` may be exported at the end as a record.
  On beauty (NixOS), Docker or Podman has to be enabled. The analyses'
  dependencies (limma, MetaProViz, SummarizedExperiment, openxlsx, httr2,
  readxl, here, `omnipath_metabo`, pandas, requests, duckdb, …) are not yet
  declared anywhere.
- **Database versioning.** Queries stay live during development. At
  publication, the database is deposited as parquet files (a few GB). The
  provenance should record the database version, and the pipeline should be
  able to rebuild from the parquet deposit. RDKit structure queries run inside
  Postgres, so their results have to be exported as well.
- **Azelate supplementary workbook** becomes a `tables/` artifact built by
  `rebuild.sh` from `data/derived/azelate/` (with Daniele).
- **Large files in git** (COSMOS+ CSVs 49 MB, PKN CSV 13 MB, xlsx): left as
  they are; Git LFS or Zenodo to be decided.
- **R package `metabo.figures`: left as it is.** Options for later (a team
  decision):
  1. move it to `code/metabo.figures/` and load it with
     `pkgload::load_all()` in `rebuild.R`;
  2. slim it down (drop generated `man/` and `R CMD check`, keep lint and
     tests);
  3. replace it with plain sourced R files.
  `lib/log.sh` stays until then.

## Picking up with Spec Kit

The first preprint was developed from the spec
`saezverse/ai/specifications/omnipath-metabo-figures/specs/001-figures-pipeline/`.
Its requirement and task ids, figure numbers and paths refer to the old layout.
A new specification (e.g. for the next version of the cancer cell lines use
case) should:

- start from this branch (or `main` after the merge), not from `preprint_v0`;
- refer to folders and functions by their content names, and to numbers only
  through `manuscript.yaml`;
- follow the flow `data/raw` → `analyses/<name>` → `data/derived/<name>` →
  `figures/<id>`;
- decide where specs live now: in this repository (`specs/`, `.specify/`) or
  in saezverse as before.

Old spec ids as referenced in the code and READMEs, mapped to the new names:

| Old spec id | What | Now |
|---|---|---|
| FR-005, FR-005a | Architecture diagram as fingerprinted manual asset | `figures/architecture/`, `figures/architecture/manual/` |
| FR-043 | Statistics digest for the architecture panel | `figures/architecture/panel-a-stats/`, `R/stats-*.R` |
| FR-007a | Faceted resource overview | `resource_overview_*`, `plot_resource_overview*` (`figures/database-content/`) |
| FR-007b | Coverage profile | `coverage_profile`, `plot_coverage_profile` |
| FR-007c | Resource-overlap networks | `resource_overlap*`, `plot_resource_overlap` |
| FR-007d | Entity × interaction-type matrix | `entity_by_interaction_type`, `plot_entity_interaction_matrix` |
| FR-007e | Structural specificity × chemical category | `specificity_by_category`, `plot_specificity_by_category` |
| FR-007f | RaMP conflicts | `ramp_conflict_counts`, `plot_ramp_conflict` |
| FR-010a–e | MetaLinksDB comparison panels | `metalinks_*_panel` (`figures/metalinks-versions/`) |
| FR-011, FR-011a | Old COSMOS PKN vs COSMOS+ | `cosmos_*_panel`, `cosmos_old_pkn`, `cosmos_plus_data` (`figures/cosmos-pkn/`) |
| FR-012 | Lung-cancer use case | `figures/cancer-cell-lines/`, `analyses/cancer-cell-lines*/`, `analyses/azelate/` |
| FR-015 | RaMP comparison table | `tables/ramp-comparison/`, `tbl_ramp_comparison_*` |
| FR-015a | Record-coverage table | `tables/record-coverage/`, `record_coverage_*` |
| FR-029 | Manuscript bundle | `assemble_bundle()`, `build/manuscript-bundle.pdf`, ordered by `manuscript.yaml` |
| T058–T062, T069 | Refactor of the use-case Rmds into the pipeline | open; the Rmds now live in `analyses/cancer-cell-lines/` |
| T066 | COSMOS PKN from the OmniPath Metabo API instead of vendored CSVs | open; see `analyses/cancer-cell-lines-cosmos/` |

Open items that a new spec for the cancer cell lines use case will likely
cover:

- re-running `analyses/cancer-cell-lines/` in the pinned environment, with new
  contrasts;
- replacing the cached intermediates in `03_feature_processing.Rmd`;
- deciding with Yunfan which COSMOS-PKN scripts are still needed and which DEM
  selection rule applies;
- T066 (PKN from the API) and the azelate workbook as a `tables/` artifact.
