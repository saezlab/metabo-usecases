# fig04-cosmos-pkn — COSMOS Prior-Knowledge Network Comparison

**Mode**: mixed (pipeline panels B–C; manual schematics A, D–F)

## Purpose

Characterises the OmniPath Metabo prior-knowledge network (PKN) used by COSMOS
and compares it to the legacy cosmosR PKN included with cosmosR v1.18.1.  The
figure is organised as a two-row composite:

- **Panel A** (manual): Architecture schematic — how the full OmniPath Metabo
  network is collapsed to a signed binary network for COSMOS.
- **Panel B** (pipeline): Grouped-bar chart — edge counts by interaction-type
  category, old PKN vs new OmniPath Metabo PKN side-by-side.
- **Panel C** (pipeline): Proportional stacked bar — relative composition of
  each PKN by category, enabling "what fraction is metabolic vs signaling" at a
  glance.
- **Panel D** (manual): Regulation-type legend schematic.
- **Panel E** (manual): MOON activity-inference diagram.

## Artifact IDs

| ID | Kind | Mode |
|----|------|------|
| `fig04-cosmos-pkn` | figure-composite | mixed |
| `fig04-cosmos-pkn/panelB` | figure-panel | pipeline |
| `fig04-cosmos-pkn/panelC` | figure-panel | pipeline |

## Deployment

Default: `dev3` (only edge counts queried; no RDKit dependency).

## Outputs

| File | Description |
|------|-------------|
| `out/fig04-cosmos-pkn.pdf` | Full composite (xelatex assembled) |
| `out/fig04-cosmos-pkn.svg` | SVG composite |
| `out/fig04-cosmos-pkn-with-caption.pdf` | Composite + typeset caption |
| `out/panelB.pdf` / `panelB.svg` | Pipeline panel B |
| `out/panelC.pdf` / `panelC.svg` | Pipeline panel C |
| `out/fig04-cosmos-pkn.pdf.provenance.json` | Provenance sidecar |

## Dependencies

- `R/data-cosmos_old.R` — `cosmos_old_pkn()` loader for vendored meta_network.RData
- `R/plots-version_comparison.R` — `cosmos_old_vs_new()` renderer
- `inst/extdata/cosmos/meta_network.RData` — vendored old PKN
- `manual/pkn-to-binary-network.pdf`, `manual/regulation-types.pdf`,
  `manual/moon-activity-inference.pdf` — hand-authored schematics (PLACEHOLDERS)
