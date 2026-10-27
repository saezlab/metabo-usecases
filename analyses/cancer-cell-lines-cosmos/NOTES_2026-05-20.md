# Session Notes: DEM Filtering & COSMOS PKN Connection
**Date:** 2026-05-20 (updated 2026-05-21)  
**Project:** CCLE Case Study — EGFR/KRAS lung mutation metabolomics

---

## Data

- **Input file:** `data/ExtendedDataTable_DifferentialAnalysis_Shorthouse_LungMutations.xlsx`
- **Sheets used:** `EGFR_filt_limma`, `KRAS_filt_limma`
- **Original correction:** BH (too strict; 0 significant hits in both)
- **Decision:** Use **top 10 up-regulated + top 10 down-regulated by t-statistic** (20 DEMs per comparison, no p-value threshold)

---

## Scripts

### `script/01_cosmos_pkn.py`
Builds the COSMOS PKN using `omnipath_metabo`:

```python
from omnipath_metabo.datasets import cosmos

allosteric = cosmos.build_allosteric()
df_a = pd.DataFrame(allosteric.network)

enzyme_met = cosmos.build_enzyme_metabolite(gem={'include_orphans': False})
df_e = pd.DataFrame(enzyme_met.network)
```

- Must be run from `~/biotools/omnipath-metabo/` with `uv run python`
- PubChem name lookup (~4124 entries) takes ~40 min on first run
- Outputs saved to `data/`:
  - `pkn_allosteric.csv` — 12,978 edges (BRENDA: 3,786 + STITCH: 9,192)
  - `pkn_enzyme_metabolite.csv` — includes Human-GEM, Recon3D, and **KEGG** (added 2026-05-21); KEGG edges have `resource == 'KEGG'` and pass the `build_enzyme_metabolite` row filter alongside GEM edges; overlapping edges across sources are merged with `;`-joined resource strings (e.g. `GEM:Human-GEM;KEGG`)

### `script/02_connect_dem_pkn.R`
Loads DEMs, explodes semicolon-separated ChEBI IDs, and screens PKN for matching edges.

**Key logic:**
- DEM selection: `slice_max(t, n=10)` (up) + `slice_min(t, n=10)` (down) → 20 DEMs per comparison; arranged by `desc(t)`
- ChEBI column contains multiple IDs separated by `;` — must `strsplit` before matching
- Match: `source %in% chebi_ids | target %in% chebi_ids`
- Separate tables for allosteric and enzyme-metabolite PKN (enzyme-metabolite now includes KEGG)
- Separate sections for EGFR and KRAS (no loop)
- No CSV output; display only via `knitr::kable()`

### `script/03_visualization.R`
Two figures using `ggplot2` + `patchwork`.

---

## Figure Formatting Guidelines

- **Y-axis length:** 8 cm (fixed plot area)
- **X-axis length:** 0.5 cm × bars
- **Grid:** none (`panel.grid.major/minor = element_blank()`)
- **Axis tick labels:** Arial 9
- **Axis labels:** Arial 10
- **Legend title:** Arial 10, bold (Fig 1) / Arial 7, bold (Fig 2)
- **Legend entries:** Arial 10 (Fig 1) / Arial 7 (Fig 2)
- **Legend key size:** 0.4 cm (Fig 1) / 0.25 cm (Fig 2)
- **PDF device:** `quartz(type = "pdf", family = "Arial")` — cairo_pdf not available (no X11)

### Colour palette (in order)

| # | Hex | Assigned to |
|---|---|---|
| 1 | `#006384` | STITCH |
| 2 | `#9F0162` | BRENDA |
| 3 | `#FEAF16` | GEM:Human-GEM |
| 4 | `#BBCC33` | GEM:Recon3D |
| 5 | `#EA6572` | GEM:Human-GEM;GEM:Recon3D |
| 6 | `#009E73` | KEGG |
| 7 | `#99DDFF` | GEM:Human-GEM;KEGG |
| 8 | `#D03293` | GEM:Recon3D;KEGG |
| 9 | `#984EA3` | GEM:Human-GEM;GEM:Recon3D;KEGG |
| — | `#BEBEBE` | None / Unknown / NA |

Rules: use first N colours for N categories; grey (`#BEBEBE`) for unknown/NA; if only one colour needed, use grey.

---

## Figure 1 — `figure/fig1_DEM_PKN_groups.pdf`

- **Size:** 25.4 × 21.8 cm
- **Structure:** 2 stacked panels (KRAS top, EGFR bottom)
- **X-axis:** 4 groups — `up - Allos`, `down - Allos`, `up - GEM`, `down - GEM`
- **Y-axis:** Number of PKN edges
- **Fill:** Resource (9 categories, all palette colours)
- **Legend:** 2 rows, collected at bottom

**Suggested title:** *PKN connectivity of differentially expressed metabolites in KRAS- and EGFR-mutant lung cancer cell lines*

---

## Figure 2 — `figure/fig2_DEM_PKN_locations.pdf`

- **Size:** 21.2 × 11.5 cm
- **Structure:** 1 row, 4 panels — KRAS up | KRAS down | EGFR up | EGFR down
- **X-axis:** Subcellular location (cytoplasm, mitochondria, nucleus, ER, extracellular, Golgi, vesicle, peroxisome, lysosome)
- **Y-axis:** Number of PKN edges
- **Fill:** Resource (9 categories, all palette colours)
- **Legend:** 1 row, collected at bottom; key 0.25 cm, font size 7

**Suggested title:** *Subcellular compartment distribution of PKN interactions for differentially expressed metabolites*

---

## Open Questions / To-Do

- **GEM bidirectional edges:** Reversible GEM reactions appear as two rows in the PKN (one per direction: `reverse: False` enzyme→metabolite, `reverse: True` metabolite→enzyme). Currently counted as **2 separate edges**. Decision pending: keep both (directionally meaningful for COSMOS/MOON) or deduplicate by reaction ID.
- Figure 2 title not finalised.
