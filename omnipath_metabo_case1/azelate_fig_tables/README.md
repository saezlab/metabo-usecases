# Azelate OmniPath FastAPI Queries, Figures, and Supplementary Tables

This folder contains the azelate OmniPath metabo FastAPI workflow, including plots and tables generation code.

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
notebooks/
  azelaic_acid_query_tables.py
  azelaic_acid_query_tables.ipynb
  azelaic_acid_plots.py
  azelaic_acid_plots.ipynb

outputs/
  azelaic_acid_tables.xlsx
  manifest.json
  csvs/
  figures/
```

## Main Outputs

The main supplement workbook is:

```text
outputs/azelaic_acid_tables.xlsx
```

Workbook sheets:
- `legend`: Sheet descriptions and row/column counts.
- `summary`: Collapsed entity-resolution and count summary for the three queried entities.
- `entity_annotations`: Seed entity annotations, including source memberships, external IDs, names, formulas, structures, and other identifiers.
- `relations`: Human-readable direct relations for each queried seed entity.
- `hmdb_classes`: HMDB/ChemOnt class relations.
- `food_occurrence`: FooDB/food occurrence relations.
- `metabolic_reactions`: Direct metabolic reaction or transport relations.
- `macdb_associations`: MACDB disease-association evidence in long format.
- `fig3_panel_C_data`: Resolved-only data for the interaction-type/source stacked barplot.
- `fig3_panel_D_data`: Resolved-only cancer association counts by disease type and sample type.
- `fig_mindmap_data`: Resolved-only node/edge data for the mindmap.

The `csvs/` folder contains one CSV per workbook data sheet.

## Figures

Generated figures are in:

```text
outputs/figures/
```

Main figure panels:
- `fig3_panel_C_interaction_type_barplot.png`
- `fig3_panel_D_cancer_associations_sample_type.png`

Extra figures:
- `extra_evidence_atlas_mindmap.png`
- `extra_target_type_composition_by_source.png`

`figure_index.csv` and `figure_index.json` list generated figure paths and labels.

## Regenerate Tables And Workbook

From this folder:

```bash
python3 -u notebooks/azelaic_acid_query_tables.py
```

This regenerates:
- `outputs/csvs/*.csv`
- `outputs/azelaic_acid_tables.xlsx`
- `outputs/manifest.json`

The script uses live FastAPI calls, so network access to `https://dev.omnipathdb.org/api` is required.

## Regenerate Figures

Run the figure script after regenerating the tables:

```bash
python3 -u notebooks/azelaic_acid_plots.py
```

This reads the curated CSVs and writes figures into:

```text
outputs/figures/
```

Plots use only rows marked `resolved`. Unresolved seed records remain in the supplement workbook/CSVs.

## Run Notebooks

From this folder:

```bash
jupyter nbconvert --execute --to notebook --inplace notebooks/azelaic_acid_query_tables.ipynb --ExecutePreprocessor.timeout=1200
jupyter nbconvert --execute --to notebook --inplace notebooks/azelaic_acid_plots.ipynb --ExecutePreprocessor.timeout=600
```

Run the extraction notebook before the visualization notebook.

## Current Run Summary

Current generated summary:

| Query seed | Status | Interactions | Evidence records |
|---|---:|---:|---:|
| `Azelaic acid (PTF52203)` | `unresolved` | 130 | 130 |
| `Azelaic acid; nonanedioic acid` | `unresolved` | 2 | 4 |
| `Azelate` | `resolved` | 51 | 92 |

