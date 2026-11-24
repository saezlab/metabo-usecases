# metalinks-versions

MetaLinksDB 2.0 compared with MetaLinksDB 1.0: coverage of interactions,
metabolites and proteins, metabolite and protein classes, and relationship
types.

**Status (2026-11-24): builds** (about 3 min), with layout problems: the
x-axis labels of the metabolite-class panel overlap, and the legend of the
relationship-type panel is cut off at the right. MetaLinksDB 2.0 is read
from
`custom_views.metalinksdb_relations` on the `prod` Postgres deployment
(June 2026 data, opt-in). dev5 does not have that view. MetaLinksDB will become a
subset of the Parquet release (`docs/DB_ACCESS.md`); move the query there
when it lands.

## Inputs

| Input | What |
|---|---|
| `prod` Postgres, `custom_views.metalinksdb_relations` | MetaLinksDB 2.0 |
| OmnipathR (`metalinks_v1_snapshot()`) | MetaLinksDB 1.0 from the SQLite shipped with OmnipathR |
| `data/vendored/mpi-baselines/` | Optional normalized snapshots of other resources (CellPhoneDB, scConnect, STITCH). None exist; the build leaves those resources out and records that in the sidecar. |

## Outputs (`out/`, gitignored)

| File | Panel |
|---|---|
| `coverage` | Interactions, metabolites and proteins per resource |
| `metabolite_classes` | Coverage by metabolite class |
| `protein_classes` | Coverage by protein class |
| `metalinks_overview` | MetaLinksDB 2.0 by contributing resource |
| `relationship_types` | Transport, receptor and interaction relations |
| `metalinks-versions` | The composite |
| `supplementary/` | Per-panel CSVs and `legend.csv` |

Each figure as PDF and SVG, plus the captioned PDF, `caption.txt` and a
provenance sidecar.

## Run

```sh
./rebuild.sh metalinks-versions
```
