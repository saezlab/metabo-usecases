# cosmos-pkn

The COSMOS+ prior-knowledge network: which regulation types it encodes,
how it compares with the old COSMOS PKN, its coverage of subcellular
compartments, and the contribution of each resource.

**Status (2026-11-24):** reads vendored COSMOS+ exports of June 2026 and
MetaLinksDB 2.0 from dev5 (`custom_views.metalinksdb_relations`, which dev5
no longer has, so the MetaLinksDB panel fails). COSMOS and MetaLinksDB will
become subsets of the Parquet release (`docs/DB_ACCESS.md`); re-vendor or
read from there when they land.

## Panels

| Composite | Output (`out/`) | Shows | Source |
|---|---|---|---|
| A | — | Regulation-type schematic | `manual/regulation-types.png` |
| B | `panel_a` | Old COSMOS PKN vs COSMOS+ by interaction type | `cosmos_comparison_panel()` |
| C | `panel_b` | COSMOS+ interactions per compartment | `cosmos_compartment_panel()` |
| D | `panel_c_split` | Metabolites and proteins per resource | `cosmos_resource_contribution_panel()` |
| — | `panel_c` | Variant of D, stacked | |
| — | `panel_d` | MetaLinksDB 2.0 vs COSMOS+ (not in the composite) | `metalinks_cosmos_panel()` |

Plus `cosmos-pkn.{pdf,svg}` (composite, 180 × 200 mm), the captioned PDF,
`caption.txt` and a provenance sidecar. `composition.yaml` lists the
panels.

## Inputs

| Input | What |
|---|---|
| `data/vendored/cosmos/meta_network.RData` | Old COSMOS PKN (cosmosR 1.18.1) |
| `data/vendored/cosmos/cosmos_plus_{human,mouse}.csv` | COSMOS+ exports |
| dev5 Postgres | MetaLinksDB 2.0, and the snapshot id for provenance |
| `manual/regulation-types.{pdf,png}` | Schematic, see `manual/README.md` |

## Run

```sh
./rebuild.sh cosmos-pkn
```

## Known issues

- `cosmos_old_pkn()` installs cosmosR from GitHub at build time instead
  of reading the vendored `meta_network.RData` (see `RESTRUCTURING.md`).
- `cosmos_plus_data()` takes more than 40 minutes on the 29 MB CSV in
  the tests.
