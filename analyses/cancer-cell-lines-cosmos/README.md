# Cancer cell lines: COSMOS PKN

Extension of the cancer cell lines use case (Yunfan Bai): builds the COSMOS
prior-knowledge network (PKN) with `omnipath_metabo` and connects the
differentially abundant metabolites (DEMs) from
`analyses/cancer-cell-lines/04_differential_analysis.Rmd` to it. Working
notes: `NOTES_2026-05-20.md`.

## Scripts

| Script | Reads | Writes |
|---|---|---|
| `01_cosmos_pkn.py` | DA xlsx (`data/derived/cancer-cell-lines/`) | `data/derived/cancer-cell-lines-cosmos/pkn_allosteric.csv`, `pkn_enzyme_metabolite.csv` |
| `02_connect_dem_pkn.R` | DA xlsx + PKN CSVs | nothing (prints tables) |
| `03_visualization.R` | DA xlsx + PKN CSVs | `output/` (gitignored) |

Figure 6 (panels E, F) reads the two PKN CSVs through
`case_study_cosmos_pkn()`. Its plots were ported from `03_visualization.R`
into `R/plots-stacked_bars.R`.

`01_cosmos_pkn.py` needs the `omnipath_metabo` package. It was run from an
`omnipath-metabo` checkout, e.g.:

```bash
cd ~/biotools/omnipath-metabo
uv run python <repo>/analyses/cancer-cell-lines-cosmos/01_cosmos_pkn.py
```

The first run takes about 40 min (PubChem name lookup).

## Open questions (for review)

- Are `02_connect_dem_pkn.R` and `03_visualization.R` still needed, given
  that Figure 6 now uses `R/plots-stacked_bars.R`?
- `03_visualization.R` opens `quartz()` devices, which exist only on macOS.
- `01_cosmos_pkn.py` filters DEMs by p < 0.05 and |logFC| > 0.5, while the
  notes and `02` use the top 10 up and top 10 down by t-statistic. Which
  selection is meant?
- Which `omnipath_metabo` version or commit was used? It should be pinned
  in the environment.
