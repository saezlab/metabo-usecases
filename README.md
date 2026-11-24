# metabo-usecases

Code and data behind the figures and tables of the **OmniPath Metabo**
manuscript: the use-case analyses, and a pipeline that rebuilds every figure
and table from the OmniPath data and the analysis results.

The tag [`preprint_v0`](https://github.com/saezlab/metabo-usecases/tree/preprint_v0)
marks the state used for the first preprint. `RESTRUCTURING.md` records
what changed since, `AGENTS.md` the conventions for working here.

## How the repository works

```
data/raw/<name>/          external inputs, README with the source
      │
      ▼
analyses/<name>/          use-case analyses (R Markdown, Python, notebooks)
      │
      ▼
data/derived/<name>/      analysis results, written only by analyses/
      │                                  OmniPath data (docs/DB_ACCESS.md)
      ▼                                           │
figures/<name>/, tables/<name>/   build.R  ◄──────┘   data/vendored/
      │
      ▼
figures/<name>/out/, tables/<name>/out/   PDF/SVG/CSV + provenance (gitignored)
```

- **Analyses** do the substantive work and are run by their own entry
  script (see each folder's README). They write the files that figures and
  tables read.
- **Figures and tables** are built by one entry point, `./rebuild.sh`,
  which runs every `figures/*/build.R` and `tables/*/build.R` in the order
  of `manuscript.yaml`. Shared code (data loaders, plots, styles,
  provenance) is the R package `metabo.figures` in `R/`.
- **OmniPath data** comes from the Parquet release of OmniPath (default for
  new work) or, for now, the dev5 Postgres instance; see
  `docs/DB_ACCESS.md`.

## Folders

| Folder | Contents |
|---|---|
| `analyses/` | Use-case analyses, one folder each |
| `data/raw/`, `data/derived/`, `data/vendored/` | External inputs, analysis results, pinned third-party snapshots; each input folder has a README with its source |
| `figures/`, `tables/` | One folder per figure or table: `build.R`, `caption.tex`, `README.md`, `manual/` for hand-made assets |
| `manuscript.yaml` | Figure and table numbers; the only place they are defined |
| `R/`, `man/`, `tests/`, `DESCRIPTION`, `NAMESPACE` | R package `metabo.figures` and its tests |
| `inst/extdata/` | Package resources: palettes, logos, JSON schema, connection template |
| `python/` | Python environment (uv) and helpers |
| `tex/`, `lib/` | LaTeX caption style, bash logging helper |
| `docs/` | `DB_ACCESS.md` (backends), `CONFIGURATION.md` (pipeline connection), `CODING_STYLE.md` |
| `build/` | Pipeline output (gitignored): logs, build manifests, manuscript bundle |

### Analyses

| Folder | Content |
|---|---|
| `analyses/cancer-cell-lines/` | Shorthouse 2022 cancer cell line metabolomics: pre-processing, mutation status, contrasts, differential analysis, MetSigDB enrichment |
| `analyses/ramp-ambiguity/` | Exploration of ambiguous RaMP ID mappings (not runnable at the moment) |

### Figures and tables

| Folder | Content |
|---|---|
| `figures/architecture/` | Architecture diagram (manual asset) and database statistics |
| `figures/database-content/` | Database content overview |
| `figures/metalinks-versions/` | MetaLinksDB 2.0 vs 1.0 |
| `figures/cosmos-pkn/` | COSMOS+ prior-knowledge network vs the old COSMOS PKN |
| `figures/cancer-cell-lines/` | Cancer cell lines use case (not built yet) |
| `tables/id-resolving/` | Identifiers per resource and type |
| `tables/ramp-comparison/` | RaMP InChIKey conflicts |
| `tables/record-coverage/` | Resource × record-type coverage |

Each folder's README states what it shows, its inputs and its current
state.

## Running

The OmniPath data lives on `beauty` (loopback ports and local files), so
the pipeline runs there or through an SSH tunnel; see `docs/DB_ACCESS.md`.

```sh
. ~/.config/metabo-figures/env.sh          # on beauty: connection settings
nix-shell shell.nix                        # LaTeX, pdftk, libpq
./rebuild.sh --dry-run                     # list what would be built
./rebuild.sh metalinks-versions            # build one figure
./rebuild.sh --png --bundle --check        # build everything, with checks
tail -f build/logs/latest.log              # follow the pipeline log
```

Analyses run with their own entry script, e.g.
`Rscript analyses/cancer-cell-lines/run_analysis.R`.

## Environment

R dependencies of the pipeline are listed in `DESCRIPTION`; Python is
managed with uv (`python/pyproject.toml`, `python/uv.lock`); system tools
(LaTeX, libpq, pdftk) come from `shell.nix` on beauty. The analyses' R
dependencies are not declared yet. A pinned Docker image is planned
(`RESTRUCTURING.md`, Planned work).

## Conventions

See `AGENTS.md` and `docs/CODING_STYLE.md`. In short: every file in
`data/derived/` is written by a script in `analyses/`; no personal absolute
paths; names by content, numbers only in `manuscript.yaml`; every built
artifact gets a provenance sidecar; trunk-based development on `main`.

## License

BSD 3-Clause, see `LICENSE`.
