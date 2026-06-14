# fig05-cosmos-pkn — COSMOS+ Prior-Knowledge Network (FR-011)

**Mode**: mixed (four pipeline panels A–D + three manual schematic
graphics S1–S3 included in the composite).

## Purpose

Characterises the COSMOS+ prior-knowledge network produced by the
OmniPath Metabo build and compares it to the legacy cosmosR PKN
(v1.18.1) and to MetaLinksDB 2.0. See FR-011 family in
`specs/001-figures-pipeline/spec.md`.

## Panels

- **Panel A** (FR-011a, pipeline): species-aware grouped-bar
  comparison of old COSMOS PKN vs. COSMOS+ (human and mouse
  separately) by interaction-type category.
- **Panel B** (FR-011b, pipeline): COSMOS+ interactions per
  annotated subcellular compartment.
- **Panel C** (FR-011c, pipeline): unique metabolite (ChEBI) and
  protein (UniProt) entity counts per contributing resource in
  COSMOS+ (top-15 + Other).
- **Panel D** (FR-011d, pipeline): MetaLinksDB 2.0 vs. COSMOS+
  comparison by canonical interaction category (Transport, Ligand
  receptor, Catalysis, Gene regulation, Allosteric regulation).
- **Schematic S1** (manual): PKN-to-binary-network construction.
- **Schematic S2** (manual): regulation-types legend.
- **Schematic S3** (manual): MOON recursive-activity-inference
  diagram.

## Artifact IDs

| ID | Kind | Mode |
|----|------|------|
| `fig05-cosmos-pkn` | figure-composite | mixed |
| `fig05-cosmos-pkn/panel_a` | figure-panel | pipeline |
| `fig05-cosmos-pkn/panel_b` | figure-panel | pipeline |
| `fig05-cosmos-pkn/panel_c` | figure-panel | pipeline |
| `fig05-cosmos-pkn/panel_d` | figure-panel | pipeline |

## Deployment

`dev5` — the dev5 integrated build provides
`custom_views.metalinksdb_relations` for Panel D and is the only
deployment in the post-2026-06-14 active rotation.

## Outputs

| File | Description |
|------|-------------|
| `out/fig05-cosmos-pkn.pdf` | Full composite (pipeline + schematics, xelatex-assembled). Skipped while `manual/*.gitkeep` placeholders are in place. |
| `out/fig05-cosmos-pkn-pipeline.{pdf,svg}` | Pipeline-only composite (4 panels via patchwork) — the artifact actually built today. |
| `out/panel_{a,b,c,d}.{pdf,svg}` | Individual pipeline panels. |
| `out/fig05-cosmos-pkn-pipeline.pdf.provenance.json` | Provenance sidecar — dev5 build_id, vendored COSMOS+ CSV fingerprints, old-COSMOS species assumption. |

## Dependencies

- `R/data-cosmos_old.R` — `cosmos_old_pkn()` loader for the
  vendored `meta_network.RData`.
- `R/data-cosmos_plus.R` — `cosmos_plus_data()` for the vendored
  `cosmos_plus_{human,mouse}.csv`.
- `R/plots-version_comparison.R` — the four `fig04_*` renderer
  functions (legacy naming — kept until a deliberate package-wide
  function rename).
- `inst/extdata/cosmos/{meta_network.RData,cosmos_plus_*.csv}` —
  vendored input data.
- `manual/pkn-to-binary-network.pdf`,
  `manual/regulation-types.pdf`,
  `manual/moon-activity-inference.pdf` — hand-authored schematics
  (currently PLACEHOLDERS; the full-composite step is skipped with
  a warning until they land).
