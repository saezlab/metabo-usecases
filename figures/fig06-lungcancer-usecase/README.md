# fig06-lungcancer-usecase — Lung-cancer use case (FR-012)

**Mode**: mixed (four pipeline panels A/B/E/F + one manual webapp
screenshot at Panel C + one DB-dependent panel D currently stubbed).

## Purpose

Refactor of the legacy `omnipath_metabo_case1/` Rmd + script chain into
the figures-pipeline R-package layout. See spec FR-012 family in
`specs/001-figures-pipeline/spec.md`.

## Panels

- **Panel A** (pipeline): KRAS volcano plot — log fold-change vs.
  −log₁₀ p-value for the KRAS-mutated vs. wild-type comparison on
  the Shorthouse 2022 lung-cancer cohort.
- **Panel B** (pipeline): EGFR volcano plot — same on the EGFR
  contrast. X and Y axis limits aligned with Panel A (FR-019).
- **Panel C** (manual): screenshot of the `metabo.omnipathdb.org`
  azelaic-acid entry view. Currently a placeholder PNG with
  re-capture instructions — see `manual/README.md` (Assumption 3).
- **Panel D** (pipeline, **STUBBED**): differential top-hit
  connections in MetaLinksDB 2.0. Build currently emits a
  TODO-placeholder until the in-progress MetaLinksDB rebuild
  settles; the SQL against `custom_views.metalinksdb_relations`
  will be wired in once the schema is finalised.
- **Panel E** (pipeline): GEM / allosteric edge counts by direction
  (up- / down-regulated) × resource for KRAS and EGFR contrasts.
  Stacked bars coloured by COSMOS PKN resource.
- **Panel F** (pipeline): subcellular-location edge counts —
  four faceted stacked bars (KRAS up, KRAS down, EGFR up, EGFR
  down) by single-letter compartment code × resource.

## Artifact IDs

| ID | Kind | Mode |
|----|------|------|
| `fig06-lungcancer-usecase` | figure-composite | mixed |
| `fig06-lungcancer-usecase/panel_a` | figure-panel | pipeline |
| `fig06-lungcancer-usecase/panel_b` | figure-panel | pipeline |
| `fig06-lungcancer-usecase/panel_c` | figure-panel | manual |
| `fig06-lungcancer-usecase/panel_d` | figure-panel | pipeline (stub) |
| `fig06-lungcancer-usecase/panel_e` | figure-panel | pipeline |
| `fig06-lungcancer-usecase/panel_f` | figure-panel | pipeline |

## Deployment

Currently DB-independent for the actively-built panels. Panel D will
route to `dev5` once wired in (post-MetaLinksDB rebuild). Panels E/F
consume the vendored COSMOS PKN CSVs in
`omnipath_metabo_case1/data/pkn_*.csv` as a fixture until the
`build_pkn.py` refactor (T066) ports the COSMOS PKN construction to
the OmniPath Metabo API.

## Inputs

| Path | Role |
|------|------|
| `omnipath_metabo_case1/Results/Differential_Analysis/ExtendedDataTable_DifferentialAnalysis_Shorthouse_LungMutations.xlsx` | Differential analysis output — sheets `KRAS_filt_limma`, `EGFR_filt_limma` drive Panels A, B; both drive top-DEM ChEBI lookup for Panels E, F. |
| `omnipath_metabo_case1/data/pkn_allosteric.csv` | COSMOS allosteric PKN edges (vendored fixture). Drives Panels E, F. |
| `omnipath_metabo_case1/data/pkn_enzyme_metabolite.csv` | COSMOS enzyme-metabolite PKN edges (vendored fixture). Drives Panels E, F. |
| `figures/fig06-lungcancer-usecase/manual/panel_C_azelaic_acid.png` | Manual webapp screenshot — placeholder. |

## Outputs

| File | Description |
|------|-------------|
| `out/fig06-lungcancer-usecase.{pdf,svg}` | Full composite (patchwork over Panels A/B/D/E/F + injected Panel C). |
| `out/panel_{a,b,d,e,f}.{pdf,svg}` | Individual pipeline panels. |
| `out/fig06-lungcancer-usecase-with-caption.pdf` | Composite + typeset caption (xelatex + `tex/caption.sty`, FR-041). |
| `out/caption.txt` | Plain-text caption (FR-041a). |
| `out/fig06-lungcancer-usecase.pdf.provenance.json` | Provenance sidecar — input fingerprints, panel-D stub flag. |

## Dependencies

- `R/data-case_study.R` — `case_study_differential()` loader for the
  DA xlsx and `case_study_cosmos_pkn()` loader for the PKN CSVs.
- `R/plots-volcano.R` — `volcano_panel(diff_tibble, contrast_label)`.
- `R/plots-stacked_bars.R` — `gem_allosteric_panel()` (E) and
  `subcellular_location_panel()` (F).

## Files

- `build.R` — orchestration.
- `caption.tex` — FR-040a caption source.
- `manual/` — Panel C placeholder + capture README.
- `out/` — build outputs.

## Open items

- **T064** — Panel D top-hit MetaLinksDB connections; deferred until
  the in-progress MetaLinksDB rebuild settles.
- **T066** — refactor `01_cosmos_pkn.py` into a CLI helper that
  fetches the COSMOS PKN from the OmniPath Metabo API (via
  `omnipath-client`) rather than the legacy `pypath` path. Until
  this lands, Panels E/F consume the vendored CSV fixtures.
- **T058–T062, T069** — full Rmd-source refactor and parity gate.
  Out of scope for this first slice; Panels A/B currently read the
  pre-computed limma output directly.
