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

## Files

- `build.R` — orchestration. Queries the per-panel deployments, runs
  the FR-007a..f renderers, verifies the architecture-asset
  fingerprint, copies the vendored PDF into `out/panelA.pdf`, invokes
  `R/compose/mixed_source.R` (which dispatches to
  `tex/compose_fig01.tex`), and writes the provenance sidecar.
- `out/` — emitted artifacts: `fig01-overview.{pdf,svg}` (composite),
  `panel{A,B,C,D,E,F}.pdf` (constituents), and the
  `fig01-overview.pdf.provenance.json` sidecar.
- `manual/` — placeholder for any figure-local manual sub-assets
  (currently empty; the architecture asset lives under
  `inst/extdata/manual/architecture/` so it is reachable via
  `system.file()` once the package is installed).
