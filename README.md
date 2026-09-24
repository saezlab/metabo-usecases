# metabo-usecases

Code and data behind the figures and tables of the **OmniPath Metabo**
manuscript: the use-case analyses, and a pipeline that rebuilds every figure
and table from the OmniPath database and the analysis results.

> The tag [`preprint_v0`](https://github.com/saezlab/metabo-usecases/tree/preprint_v0) marks the state used for
> the first preprint. The `restructure` branch reorganises the repository
> for the next version (see `AGENTS.md` for the current conventions).

## How the repository works

```
data/raw/            external inputs (e.g. published supplementary data)
      │
      ▼
analyses/<name>/     use-case analyses (Rmd, Python, notebooks)
      │
      ▼
data/derived/<name>/ analysis results, written only by scripts in analyses/
      │                                        OmniPath Postgres (live)
      ▼                                                 │
figures/<id>/, tables/<id>/   build.R per figure/table ◄┘
      │
      ▼
figures/<id>/out/, tables/<id>/out/   PDF/SVG/CSV + provenance (gitignored)
```

- **Analyses** are run by hand, in the order given in each folder's README.
  They write the files that figures and tables read.
- **Figures and tables** are rebuilt by one entry point, `./rebuild.sh`,
  which runs every `figures/*/build.R` and `tables/*/build.R`. The shared
  code (data loaders, plots, styles, provenance) is the R package
  `metabo.figures` in `R/`.

## Folders

| Folder | Contents |
|---|---|
| `analyses/` | Use-case analyses, one folder each (see below). |
| `data/raw/` | External inputs, each folder with a README stating the source. |
| `data/derived/` | Results of the analyses; the input for figures and tables. |
| `figures/` | One folder per figure: `build.R`, `caption.tex`, `README.md`. |
| `tables/` | One folder per table, same structure. |
| `R/`, `man/`, `tests/`, `DESCRIPTION`, `NAMESPACE` | R package `metabo.figures` and its tests. |
| `inst/extdata/` | Vendored snapshots (old COSMOS PKN, COSMOS+, MPI baselines), manual assets, palettes, connection template. |
| `python/`, `tex/`, `lib/` | Python, LaTeX and bash helpers used by the pipeline. |
| `docs/` | `CONFIGURATION.md` (database connection config), `DB_ACCESS.md` (reaching the databases on beauty). |
| `logs/`, `manifests/` | Pipeline run logs and build manifests (contents gitignored). |

### Analyses

| Folder | Content | Used by |
|---|---|---|
| `analyses/cancer-cell-lines/` | Shorthouse 2022 cancer cell line metabolomics: preprocessing, Cellosaurus metadata, feature processing, differential analysis. | Figure 6 |
| `analyses/cancer-cell-lines-cosmos/` | Connects the differential metabolites to the COSMOS prior-knowledge network. | Figure 6 |
| `analyses/azelate/` | Azelate evidence retrieval from the OmniPath Metabo API; supplementary workbook. | Figure 6, supplement |
| `analyses/ramp-ambiguity/` | Exploration of ambiguous RaMP ID mappings (runs on dev2 only). | — |

### Figures and tables

| Folder | Content |
|---|---|
| `figures/fig01-architecture/` | Architecture diagram (manual asset) and database statistics |
| `figures/fig02-overview/` | Database content overview |
| `figures/fig04-metalinks-versions/` | MetaLinksDB v1 vs v2 and other metabolite–protein interaction resources |
| `figures/fig05-cosmos-pkn/` | Old COSMOS PKN vs COSMOS+ |
| `figures/fig06-lungcancer-usecase/` | Cancer cell lines use case |
| `tables/tab01-id-resolving/` | Identifier resolving across integrated resources |
| `tables/tab02-ramp-comparison/` | RaMP InChIKey conflicts |
| `tables/tab03-record-coverage/` | Resource × record-type coverage |

## Running

The OmniPath Postgres instances run on `beauty` and listen on localhost
only, so the pipeline runs there (or through an SSH tunnel). See
`docs/DB_ACCESS.md` and `docs/CONFIGURATION.md`.

```sh
./rebuild.sh --dry-run                     # list what would be built
./rebuild.sh fig06-lungcancer-usecase      # build one figure
./rebuild.sh --png --bundle --check        # build everything, with checks
tail -f logs/latest.log                    # follow the pipeline log
```

Analyses are run by hand; each folder's README gives the run order and
what the scripts read and write. Paths are resolved from the repository
root, so the working directory does not matter.

## Environment

Currently: R dependencies of the pipeline are listed in `DESCRIPTION`,
Python helpers are managed with `uv` (`python/uv.lock`), and system tools
(LaTeX, libpq, pdftk) come from `shell.nix` on beauty. The analyses'
dependencies are not yet declared. A pinned Docker image covering R,
Python and system tools is planned.

## Conventions

- Every file in `data/derived/` is written by a script in `analyses/`.
- No personal absolute paths; resolve paths from the repository root.
- Every built artifact gets a provenance sidecar; build manifests go to
  `manifests/`.
- R code follows the Saez lab R style guide
  (`human/guidelines/r-coding-style.md` in the saezverse repository).
- Trunk-based development on `main`.

## License

BSD 3-Clause, see `LICENSE`.
