# Vendored logo / stock-image sources

Per spec FR-005b + SC-010, every logo embedded in any figure is
committed to this repo so the pipeline is self-contained (no external
image fetches at build time).

| File | Source | Licence | Vendored on |
|---|---|---|---|
| `postgres.svg` | `~/omnipath_figures/stock-images/Postgresql_elephant.svg` (originally https://wiki.postgresql.org/wiki/Logo) | PostgreSQL trademark policy permits use in technical figures with attribution | 2026-05-28 |
| `python.svg` | `~/omnipath_figures/msb2021/python-logo-generic.svg` (originally https://www.python.org/community/logos/) | Python Software Foundation: free use in publications | 2026-05-28 |
| `r.svg` | `~/omnipath_figures/misc/Rlogo.svg` (originally https://www.r-project.org/logo/) | CC-BY-SA 4.0 (R Foundation) | 2026-05-28 |
| `rdkit.png` | `~/omnipath_figures/stock-images/rdkit-logo.png` (originally https://www.rdkit.org/) | RDKit project: BSD-licensed assets | 2026-05-28 |

## Still to vendor

| File | Source plan | Used by |
|---|---|---|
| `cytoscape.svg` | Source from https://cytoscape.org/press_kit.html (PNG or vector); convert to SVG if needed | Figure 1A clients band |
| `omnipath-resources.svg` (or substitute) | Project-internal — drawn for the Figure 1A "Download and processing" band | Figure 1A |
| `omnipath-build.svg` | Project-internal — drawn for the Figure 1A main column | Figure 1A |
| `omnipath-utils.svg`  | Project-internal — drawn for the Figure 1A utils column | Figure 1A |
| `omnipath-metabo.svg` | Project-internal — drawn for the Figure 1A metabo column | Figure 1A |
| `annnet.svg` | Existing repo (`saezlab/annnet`) marketing PNG converted to SVG | Figure 1A clients band |

Add an entry above for every new asset committed; never reach for the
file system at compile time (FR-005b).
