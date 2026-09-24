# fig05-cosmos-pkn — COSMOS+ Prior-Knowledge Network (FR-011)

**Mode**: mixed (one manual schematic at the top + three pipeline
data panels in a row beneath).

## Purpose

Characterises the OmniPath Metabo COSMOS+ prior-knowledge network:
the regulation-type repertoire it encodes, its coverage across
subcellular compartments, and the relative contribution of each
upstream resource. See FR-011 family in
`specs/001-figures-pipeline/spec.md`.

## Active composite (post-2026-06-18 merge)

Four panels: the regulation-types schematic (A, full width, top) plus
the three COSMOS+ data panels (B comparison, C compartments, D
resources) in a side-by-side row beneath. The three data panels share a
portrait `coord_flip` shape, so a single 3-across row fits them better
than stretching the comparison across a full-width slot. Set
`composite_layout <- 'comparison_fullwidth'` in `build.R` to instead put
the comparison full-width under the schematic with C | D as a bottom row.

| Panel | Source | Description |
|------:|:-------|:------------|
| A | Manual (`manual/regulation-types.png`) | Regulation-type schematic — visual key for the signed interaction categories the PKN carries. **Top, full width.** |
| B | Pipeline (`fig04_cosmos_comparison_panel`) | Species-aware old COSMOS vs. COSMOS+ comparison by interaction type (3 bars/group). **Data row.** |
| C | Pipeline (`fig04_compartment_panel`) | COSMOS+ interactions per annotated subcellular compartment, stacked by interaction type. **Data row.** |
| D | Pipeline (`fig04_resource_contribution_panel`, `position = "dodge"`) | Metabolite (teal) + protein (magenta) entity counts per top-15 resources, **grouped (side-by-side) bars**. **Data row.** |

One pipeline panel stays rendered as a standalone artifact, not in the
composite:

- `panel_d.{pdf,svg}` — FR-011d MetaLinksDB 2.0 vs. COSMOS+ comparison
  by interaction type (the comparison panel, formerly standalone
  `panel_a`, is now composite Panel B).

## Artifact IDs

| ID | Kind | Mode |
|----|------|------|
| `fig05-cosmos-pkn` | figure-composite | mixed |
| `fig05-cosmos-pkn/panel_a` | figure-panel | pipeline (comparison; composite Panel B) |
| `fig05-cosmos-pkn/panel_b` | figure-panel | pipeline (compartments; composite Panel C) |
| `fig05-cosmos-pkn/panel_c` | figure-panel | pipeline (resources, stacked; standalone variant) |
| `fig05-cosmos-pkn/panel_c_split` | figure-panel | pipeline (resources, semicolon-split + dodged; composite Panel D) |
| `fig05-cosmos-pkn/panel_d` | figure-panel | pipeline (MetaLinksDB; standalone) |
| `fig05-cosmos-pkn/regulation-types` | schematic | manual (composite Panel A) |

## Deployment

`dev5` — provides `custom_views.metalinksdb_relations` (used by
the still-rendered Panel D standalone) and is the sole anchor for
the snapshot identifier recorded in the provenance sidecar.

## Outputs

| File | Description |
|------|-------------|
| `out/fig05-cosmos-pkn.{pdf,svg}` | **Active composite** — schematic (A) full width on top; Panel B (comparison), C (compartments), D (resources, grouped) side by side below. 180×200 mm. |
| `out/fig05-cosmos-pkn-pipeline.{pdf,svg}` | Identical copy of the composite kept for backward compatibility with downstream consumers. |
| `out/panel_{a,b,c,c_split,d}.{pdf,svg}` | Individual pipeline panels. |
| `out/fig05-cosmos-pkn-with-caption.pdf` | Composite + typeset caption (xelatex + `tex/caption.sty`, FR-041). |
| `out/caption.txt` | Plain-text caption (FR-041a). |
| `out/fig05-cosmos-pkn.pdf.provenance.json` | Provenance sidecar — dev5 build_id, vendored COSMOS+ CSV fingerprints, old-COSMOS species assumption, schematic PNG fingerprint, Panel C position parameter. |

## Dependencies

- `R/data-cosmos_old.R` — `cosmos_old_pkn()` loader for the
  vendored `meta_network.RData`.
- `R/data-cosmos_plus.R` — `cosmos_plus_data()` for the vendored
  `cosmos_plus_{human,mouse}.csv`.
- `R/plots-version_comparison.R` — the four `fig04_*` renderer
  functions (legacy naming — kept until a deliberate package-wide
  function rename).
- `data/vendored/cosmos/{meta_network.RData,cosmos_plus_*.csv}` —
  vendored input data.
- `manual/regulation-types.{pdf,png}` — hand-authored schematic
  (the PNG variant is embedded in the composite via
  `png::readPNG` + `grid::rasterGrob` + `patchwork::wrap_elements`).

## Open items

- The other two schematics from the original FR-011 plan
  (`pkn-to-binary-network`, `moon-activity-inference`) remain as
  `.gitkeep` placeholders. They are NOT included in the active
  composite. When authored, the build.R schematic slot can be
  extended to a multi-element top row.
