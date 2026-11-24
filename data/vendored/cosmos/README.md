# Vendored COSMOS networks

Pinned inputs of `figures/cosmos-pkn/`. COSMOS will become a subset of the
OmniPath Parquet release (`docs/DB_ACCESS.md`); then COSMOS+ should be read
from there and these CSVs retired.

## meta_network.RData

The old COSMOS prior-knowledge network (metabolic reactions only), the
baseline the figure compares COSMOS+ against.

- Source: `data/meta_network.RData` of the R package cosmosR 1.18.1
  (<https://github.com/saezlab/cosmosR>)
- Retrieved: 2026-06-10
- MD5: `60c30310020aa4c47e150450d0c42b4d`

Note: `cosmos_old_pkn()` currently installs cosmosR from GitHub and uses
its `meta_network` instead of this file (see `RESTRUCTURING.md`, Known
issues).

## cosmos_plus_human.csv, cosmos_plus_mouse.csv

COSMOS+ networks for human (taxon 9606) and mouse (10090), exported with
the `cosmos-pkn` command of the former `omnipath-metabo` package:

```bash
cosmos-pkn export --all-columns --organism 9606  --no-stitch --output data/vendored/cosmos/cosmos_plus_human.csv
cosmos-pkn export --all-columns --organism 10090 --no-stitch --output data/vendored/cosmos/cosmos_plus_mouse.csv
```

Columns: `source`, `target`, `sign`, `interaction_type`, `resource`,
`source_type`, `target_type`, `locations`.

| File | SHA-256 | Rows |
|---|---|---|
| `cosmos_plus_human.csv` | `826ebe290f5d11d197625dc2f410192acda67fc4deb9f5fd95a0136f70fed4ff` | about 350,000 |
| `cosmos_plus_mouse.csv` | `239bc3ed179fef4a1d2c094bf0e114bc9b9a47bcad79543094923658872c1454` | about 253,000 |

History: first exported on 2026-06-11 from omnipath-metabo commit
`bc3ae3f7` (commit `8786c56` here); regenerated on 2026-06-17 after a fix
to the KEGG input, which added KEGG catalysis interactions (`ea73b7f`).
The omnipath-metabo commit of the regeneration was not recorded.
