# cancer-cell-lines — Cancer cell lines use case (FR-012)

**Mode**: pipeline (six pipeline panels A–F, all rendered from
vendored CSV / xlsx sources).

## Purpose

Refactor of the legacy `omnipath_metabo_case1/` Rmd + script chain (now
`analyses/cancer-cell-lines*/` and `analyses/azelate/`) into
the figures-pipeline R-package layout. See spec FR-012 family in
`specs/001-figures-pipeline/spec.md`.

## Panels

- **Panel A** (pipeline): KRAS volcano plot — log fold-change vs.
  −log₁₀ p-value for the KRAS-mutated vs. wild-type comparison on
  the Shorthouse 2022 lung-cancer cohort.
- **Panel B** (pipeline): EGFR volcano plot — same on the EGFR
  contrast. X and Y axis limits aligned with Panel A (FR-019).
- **Panel C** (pipeline): Azelate interaction-type composition by
  upstream source — stacked bar of relation counts per
  interaction type (associated_with, interacts_with, …) coloured
  by source (ChEMBL, DrugCentral, FooDB, HMDB, MACDB, MetAtlas,
  PFOCR, RHEA, STITCH, ChEBI, SwissLipids). Reads
  `data/derived/azelate/csvs/interaction_types_by_source.csv`.
- **Panel D** (pipeline): Azelate cancer associations by sample
  type — stacked bar of MACDB evidence counts per disease type
  coloured by sample (tissue) type (Blood, Serum, Plasma, …).
  Reads `data/derived/azelate/csvs/cancer_assoc_by_sample_type.csv`.
- **Panel E** (pipeline): GEM / allosteric edge counts by direction
  (up- / down-regulated) × resource for KRAS and EGFR contrasts.
  Stacked bars coloured by COSMOS PKN resource.
- **Panel F** (pipeline): subcellular-location edge counts —
  four faceted stacked bars (KRAS up, KRAS down, EGFR up, EGFR
  down) by single-letter compartment code × resource.

## Artifact IDs

| ID | Kind | Mode |
|----|------|------|
| `cancer-cell-lines` | figure-composite | pipeline |
| `cancer-cell-lines/panel_a` | figure-panel | pipeline |
| `cancer-cell-lines/panel_b` | figure-panel | pipeline |
| `cancer-cell-lines/panel_c` | figure-panel | pipeline |
| `cancer-cell-lines/panel_d` | figure-panel | pipeline |
| `cancer-cell-lines/panel_e` | figure-panel | pipeline |
| `cancer-cell-lines/panel_f` | figure-panel | pipeline |

## Deployment

Currently DB-independent: every panel reads a CSV / xlsx produced
by an analysis (`data/derived/cancer-cell-lines/`,
`data/derived/cancer-cell-lines-cosmos/`, `data/derived/azelate/`). The build still resolves `dev5`
solely to anchor the snapshot identifier in the provenance sidecar
so Figure 6 binds to the same `build_id` as the rest of the
pipeline. Panels C / D will move to live `metabo.omnipathdb.org` /
omnipath-client queries once that wiring lands; Panels E / F will
move off the vendored `pkn_*.csv` fixtures to the OmniPath Metabo
API as part of the T066 refactor.

## Inputs

| Path | Role |
|------|------|
| `data/derived/cancer-cell-lines/ExtendedDataTable_DifferentialAnalysis_Shorthouse_LungMutations.xlsx` (written by `analyses/cancer-cell-lines/04_differential_analysis.Rmd`) | Differential analysis output — sheets `KRAS_filt_limma`, `EGFR_filt_limma` drive Panels A, B; both drive top-DEM ChEBI lookup for Panels E, F. |
| `data/derived/cancer-cell-lines-cosmos/pkn_allosteric.csv` | COSMOS allosteric PKN edges (vendored fixture). Drives Panels E, F. |
| `data/derived/cancer-cell-lines-cosmos/pkn_enzyme_metabolite.csv` | COSMOS enzyme-metabolite PKN edges (vendored fixture). Drives Panels E, F. |
| `data/derived/azelate/csvs/interaction_types_by_source.csv` | Resolved-only interaction-type × source × relation-count table for Azelate; drives Panel C. |
| `data/derived/azelate/csvs/cancer_assoc_by_sample_type.csv` | Resolved-only disease-type × sample-type × evidence-count table for Azelate (MACDB); drives Panel D. |

## Outputs

| File | Description |
|------|-------------|
| `out/cancer-cell-lines.{pdf,svg}` | Compact composite — nested patchwork over rows AB / CD / E (×2) / F (×4), 180×220 mm. |
| `out/panel_{a,b,c,d,e_kras,e_egfr,f_kras_up,f_kras_down,f_egfr_up,f_egfr_down}.{pdf,svg}` | Individual pipeline panels. |
| `out/cancer-cell-lines-with-caption.pdf` | Composite + typeset caption (xelatex + `tex/caption.sty`, FR-041). |
| `out/caption.txt` | Plain-text caption (FR-041a). |
| `out/cancer-cell-lines.pdf.provenance.json` | Provenance sidecar — input fingerprints + composite font scale. |

## Dependencies

- `R/data-case_study.R` — `case_study_differential()` loader for the
  DA xlsx and `case_study_cosmos_pkn()` loader for the PKN CSVs.
- `R/data-azelate.R` — `azelate_panel_c_data()` / `azelate_panel_d_data()`
  loaders for the azelate query-table CSVs.
- `R/plots-volcano.R` — `volcano_panel()` (A / B).
- `R/plots-azelate.R` — `azelate_interaction_panel()` (C) and
  `azelate_disease_panel()` (D).
- `R/plots-stacked_bars.R` — `gem_allosteric_panel()` (E) and
  `subcellular_location_panel()` (F).

## Files

- `build.R` — orchestration.
- `caption.tex` — FR-040a caption source.
- `manual/` — Panel C placeholder + capture README.
- `out/` — build outputs.

## Open items

- **T066** — refactor `01_cosmos_pkn.py` into a CLI helper that
  fetches the COSMOS PKN from the OmniPath Metabo API (via
  `omnipath-client`) rather than the legacy `pypath` path. Until
  this lands, Panels E/F consume the vendored CSV fixtures.
- **T058–T062, T069** — full Rmd-source refactor and parity gate.
  Out of scope for the current slice; Panels A/B currently read the
  pre-computed limma output directly.
- Panels C/D currently read pre-computed CSVs generated by
  `analyses/azelate/azelaic_acid_query_tables.py`
  against the OmniPath FastAPI. A future slice will fold that
  notebook's queries into the unified pipeline through
  `omnipath-client`.
