# Azelate OmniPath FastAPI Queries, Figures, and Supplementary Tables

Azelate evidence retrieval from the OmniPath Metabo FastAPI (Daniele
Bottazzi), part of the cancer cell lines use case. It produces the data for
Figure 6 panels C and D and the azelate supplementary workbook.

## What Is Queried

The extraction queries `https://dev.omnipathdb.org/api` with exactly three seed handles:

| Query seed | Entity status | Entity PK | Notes |
|---|---|---|---|
| `Azelaic acid (PTF52203)` | `unresolved` | `834f4132-7b63-4502-d458-6dc442371796` | PTFI source-specific record. |
| `Azelaic acid; nonanedioic acid` | `unresolved` | `65c227ae-c250-c609-e519-395d821b1c62` | Recon3D source-specific record. |
| `Azelate` | `resolved` | `dae8e667-770c-5d34-a381-8fd25179bc58` | Resolved chemical entity with InChI key `BDJRBEYXGGNYIS-UHFFFAOYSA-N`. |

The two source-specific records are retained in the supplement tables as `unresolved`. Plots and figure-specific data sheets filter to the resolved `Azelate` entity only.

## Folder Contents

```text
analyses/azelate/                 code (this folder)
  azelaic_acid_query_tables.py      API queries -> tables, workbook, manifest
  azelaic_acid_query_tables.ipynb   notebook wrapper around the .py module
  azelaic_acid_plots.py             preview plots from the CSVs
  azelaic_acid_plots.ipynb          notebook wrapper around the .py module
  output/figures/                   preview plots (gitignored)

data/derived/azelate/             results (tracked)
  azelaic_acid_tables.xlsx
  manifest.json
  csvs/
```

Figure 6 reads `csvs/interaction_types_by_source.csv` (panel C) and
`csvs/cancer_assoc_by_sample_type.csv` (panel D) and draws its own panels in
R. The PNGs written by `azelaic_acid_plots.py` are only previews.

## Main Outputs

The main supplement workbook is:

```text
data/derived/azelate/azelaic_acid_tables.xlsx
```

Planned: the workbook is a paper deliverable and should become a `tables/`
artifact that `rebuild.sh` builds from the CSVs in `data/derived/azelate/`,
with a caption and provenance. Until then the query script writes it
directly.

Note: the tracked workbook and `manifest.json` date from the 2026-06-03
run, before the rename, so their sheet and table names are still
`fig3_panel_C_data` and `fig3_panel_D_data`. The CSVs have already been
renamed; the next run brings all of them in line.

Workbook sheets:
- `legend`: Sheet descriptions and row/column counts.
- `summary`: Collapsed entity-resolution and count summary for the three queried entities.
- `entity_annotations`: Seed entity annotations, including source memberships, external IDs, names, formulas, structures, and other identifiers.
- `relations`: Human-readable direct relations for each queried seed entity.
- `hmdb_classes`: HMDB/ChemOnt class relations.
- `food_occurrence`: FooDB/food occurrence relations.
- `metabolic_reactions`: Direct metabolic reaction or transport relations.
- `macdb_associations`: MACDB disease-association evidence in long format.
- `interaction_types_by_source`: Resolved-only data for the interaction-type/source stacked barplot.
- `cancer_assoc_by_sample_type`: Resolved-only cancer association counts by disease type and sample type.
- `fig_mindmap_data`: Resolved-only node/edge data for the mindmap.

The `csvs/` folder contains one CSV per workbook data sheet.

## Figures

Preview figures are written to `analyses/azelate/output/figures/`
(gitignored).

Main figure panels:
- `interaction_types_by_source.png`
- `cancer_assoc_by_sample_type.png`

Extra figures:
- `extra_evidence_atlas_mindmap.png`
- `extra_target_type_composition_by_source.png`

`figure_index.csv` and `figure_index.json` list generated figure paths and labels.

## Regenerate Tables And Workbook

From the repository root (paths are resolved relative to the scripts, so
any working directory works):

```bash
python3 -u analyses/azelate/azelaic_acid_query_tables.py
```

This regenerates, in `data/derived/azelate/`:
- `csvs/*.csv`
- `azelaic_acid_tables.xlsx`
- `manifest.json`

The script uses live FastAPI calls, so network access to `https://dev.omnipathdb.org/api` is required.

## Regenerate Figures

Run the figure script after regenerating the tables:

```bash
python3 -u analyses/azelate/azelaic_acid_plots.py
```

This reads the CSVs and writes preview figures into
`analyses/azelate/output/figures/`.

Plots use only rows marked `resolved`. Unresolved seed records remain in the supplement workbook/CSVs.

## Run Notebooks

From this folder:

```bash
jupyter nbconvert --execute --to notebook --inplace azelaic_acid_query_tables.ipynb --ExecutePreprocessor.timeout=1200
jupyter nbconvert --execute --to notebook --inplace azelaic_acid_plots.ipynb --ExecutePreprocessor.timeout=600
```

Run the extraction notebook before the visualization notebook.

## Current Run Summary

Current generated summary:

| Query seed | Status | Interactions | Evidence records |
|---|---:|---:|---:|
| `Azelaic acid (PTF52203)` | `unresolved` | 130 | 130 |
| `Azelaic acid; nonanedioic acid` | `unresolved` | 2 | 4 |
| `Azelate` | `resolved` | 51 | 92 |

