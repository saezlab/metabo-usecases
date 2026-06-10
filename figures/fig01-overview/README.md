# Figure 1 — OmniPath Metabo architecture and database content

**Mode**: mixed (vendored manual architecture asset at Panel A
+ pipeline quantitative panels at Panel B+)

## Composition

| Panel | Source | Notes |
|------:|:-------|:------|
| A | Vendored manual asset (`inst/extdata/manual/architecture/omnipath-architecture-new2026.pdf`) | FR-005 / FR-005a — SHA-256 fingerprint pinned in the asset folder's `README.md` and verified by `build.R` before include; mismatch fails the rebuild. |
| B+ | Pipeline (R/plots/db_content + R/plots/ramp_conflict) | One PDF per FR-007a–FR-007f variant; the composite-selection config (FR-007g, FR-039) picks the subset that appears in the composite. |

## Caption draft (pre–FR-040 source)

> Architecture and content of the OmniPath Metabo database. **(A)**
> Workflow architecture of the OmniPath ecosystem — upstream
> resources flow through download / processing in
> `omnipath-resources`, three parallel database builds
> (`omnipath-utils`, `omnipath-build`, `omnipath-metabo`), three web
> services (`utils.omnipathdb.org`, `omnipathdb.org`,
> `metabo.omnipathdb.org`), and the R + Python clients. The metabo
> column features MetalinksDB, the COSMOS PKN, and the server-side
> RDKit cheminformatics layer; MetalinksDB is built as a materialised
> view inside the main OmniPath Postgres database during the main
> build and exposed by the Metabo web service. **(B+)** Quantitative
> characterisation of the database contents — see the per-panel
> sub-captions in `caption.tex` once that source is added (FR-040,
> FR-040a). Counts come from the per-panel deployment matrix recorded
> in the figure's provenance sidecar (FR-030, FR-031); panels with
> structural-specificity / RaMP-conflict data come from `dev4`
> (cycle-001 protein-centric build), the rest from `dev3` (cycle-002
> gene-centric build).

The canonical caption consumed by the build will live alongside this
README as `caption.tex` (FR-040, FR-040a). The draft above is for the
folder's human reader; once the caption pipeline (T047a / T047b /
T047c) lands, the typeset `<slug>-with-caption.pdf` and `caption.txt`
artifacts will be derived from `caption.tex`.

## FR-007 panel inventory

Each panel is emitted as a standalone artifact under `out/`. The
final composite picks a subset via `composition.yaml` (FR-007g,
FR-039) — the rest stay as supplementary artifacts.

| FR | Artifact (`out/...`) | Size | Source / deployment | What it shows |
|---|---|---|---|---|
| 007a | `fr007a-overview.{pdf,svg}` | 320×200 mm | bitmap path on dev3 + dev4 (~85 s) | 6-facet resource overview (Entities, Associations, Interactions, Identifiers, Structures, Literature) × 3 magnitude bands; each (facet, resource) shows shared/unique + major-class stacked bars; 7 titled legends |
| 007b | `fr007b-coverage.{pdf,svg}` | 180×100 mm | `entity_source_count` on dev3 (+ `metabo_entity_structural_specificity` on dev4) | Coverage line graph: items present in ≥ N resources, colour-coded by variant (Entities, Molecular entities, Structures); log y |
| 007c | `fr007c-networks.{pdf,svg}` | 320×160 mm | `resource_overlap_summary` on dev3 (537 + 90 edges) | Two networks (Molecular entities, Interactions); shared Kamada-Kawai layout so the same resource sits at the same position in both; edges ≥ 100 only, log-scaled thickness |
| 007d | `fr007d-matrix.{pdf,svg}` | 180×110 mm | `facet_relation_bitmap` on dev3 (~150 ms) | Top-8 `participant_type` × 3 `interaction_class` bar plot; log y |
| 007e | `fr007e-specificity.{pdf,svg}` | 180×110 mm | `facet_entity_bitmap` on dev4 (~25 ms) | 6 chemical categories × 6 specificity levels; log y |
| 007f | `panelG.{pdf,svg}` (legacy slug) | 89×60 mm | `metabo_ramp_inchikey_conflict` on dev4 | RaMP InChIKey-conflict reasons (stereo / specificity / tautomer / similar / unrelated) |

Composite + caption pipeline:

- `out/fig01-overview.{pdf,svg}` — bare composite (Panel A vendored
  + panels B–G whose layout is in `tex/compose_fig01.tex`; the
  composite-selection still references the prior milestone's
  single-bar `panel{B..F}.pdf` until `composition.yaml` is
  updated).
- `out/fig01-overview-with-caption.pdf` — composite + typeset
  caption (xelatex + `tex/caption.sty`, FR-041).
- `out/caption.txt` — plain-text caption (deterministic strip,
  FR-041a).
- `out/fig01-overview.pdf.provenance.json` — sidecar.

## Files

- `build.R` — orchestration. Queries the per-panel deployments,
  runs the FR-007a..f renderers, verifies the architecture-asset
  fingerprint, copies the vendored PDF into `out/panelA.pdf`,
  invokes the composite assembler, writes the caption-pipeline
  artifacts, and writes the provenance sidecar.
- `caption.tex` — FR-040a source for the caption text (main title +
  per-panel sub-captions in `(a)` / `(b)` / ... form).
- `composition.yaml` — FR-039 panel-selection config consumed by
  the composite assembler.
- `out/` — every artifact emitted by the build (see inventory above).
- `manual/` — placeholder for any figure-local manual sub-assets.
