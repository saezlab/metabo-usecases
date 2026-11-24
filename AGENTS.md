# AGENTS.md

Notes for coding agents and developers working in this repository: the
code and data behind the figures and tables of the OmniPath Metabo
manuscript.

Read first: `README.md` (overview), the README of the folder you work in,
and `docs/DB_ACCESS.md` if you need OmniPath data. `RESTRUCTURING.md`
records structural changes since the first preprint (tag `preprint_v0`)
and the known issues.

## Layout

| Path | Contents |
|---|---|
| `analyses/<name>/` | Use-case analyses (R Markdown, Python, notebooks) with an entry script. Reports and QC go to `analyses/<name>/output/` (gitignored). |
| `data/raw/<name>/` | External inputs, unchanged, with a README stating the source. Large or unlicensed inputs stay local; the README says how to get them. |
| `data/derived/<name>/` | Results written by the scripts in `analyses/<name>/`. |
| `data/vendored/<name>/` | Pinned third-party snapshots, with a README (source, version, fingerprint). |
| `figures/<name>/`, `tables/<name>/` | One folder per figure or table: `build.R`, `caption.tex`, `README.md`, `manual/` for hand-made assets. Outputs go to `out/` (gitignored). |
| `manuscript.yaml` | Figure and table numbers; `rebuild.R` orders builds and the bundle by it. `null` = not built here. |
| `R/`, `man/`, `tests/`, `DESCRIPTION`, `NAMESPACE` | R package `metabo.figures`: data loaders, plots, styles, provenance shared by the builds. |
| `python/` | The repository's Python environment (uv) and helpers. |
| `tex/`, `lib/` | Caption style for LaTeX; bash logging helper. |
| `inst/extdata/` | Package resources: palettes, logos, JSON schema, connection template. |
| `docs/` | `DB_ACCESS.md`, `CONFIGURATION.md`, `CODING_STYLE.md`. |
| `build/` | Pipeline output (gitignored): `logs/latest.log`, `manifests/`, `manuscript-bundle.pdf`. |

## How work flows

```
data/raw ─► analyses/<name> ─► data/derived/<name> ─► figures|tables/<name>/build.R ─► out/
                 ▲                                              ▲
         OmniPath data (docs/DB_ACCESS.md)          data/vendored, OmniPath data
```

- **Analyses do the substantive work**: queries, filtering, counting,
  statistics, and every decision, documented in notebooks with the plots
  and text that support it. They write plain results (TSV, RDS, Parquet)
  to `data/derived/<name>/`.
- **Figures and tables only present results**: `build.R` reads
  `data/derived/`, `data/vendored/` or OmniPath data, draws or formats,
  adds the caption, and writes a provenance sidecar. When a figure or
  table needs non-trivial counting or decisions, put that part into an
  analysis and let `build.R` read its output. The current database-backed
  figures and tables still query Postgres directly through `R/data-*.R`;
  move them to this split when they are updated.

### Adding or updating an analysis

1. Create `analyses/<name>/` with an entry script (e.g.
   `run_analysis.R`), a `README.md`, and `config/` if it has settings.
   Raw inputs go to `data/raw/<name>/` with a README stating the source.
2. Write results to `data/derived/<name>/`; record inputs and versions
   (e.g. a manifest JSON per step, as in `analyses/cancer-cell-lines/`).
3. Decision-relevant steps get a notebook that shows the evidence (plots,
   tables) and states the decision. Keep side audits separate from the
   steps that produce the final results.

### Adding or updating a figure or table

1. Folder `figures/<name>/` or `tables/<name>/`, named by content:
   `build.R`, `caption.tex`, `README.md`.
2. Shared code goes into the package (`R/data-*.R`, `R/plots-*.R`,
   `R/tables-*.R`) with roxygen documentation; export it in `NAMESPACE`
   and add its `man/` page.
3. Set its number in `manuscript.yaml`; nowhere else.
4. Build with `./rebuild.sh <name>` and check `out/` and the log.

### READMEs

Every folder in `analyses/`, `figures/`, `tables/`, `data/raw/` and
`data/vendored/` has a short README: what it is, a **Status** line with
the date (does it build or run, on which backend, known problems),
inputs, outputs, how to run. Update the status line when you change the
folder.

## OmniPath data

Details and commands: `docs/DB_ACCESS.md`. In short:

- **Default for new work: the Parquet release** (`2026.11.19`), read-only
  as one DuckDB database on beauty, or through the clients
  `omnipath.client` (R) / `omnipath_client` (Python) from the `monorepo`
  branch of saezlab/pypath, which is the ground truth for the API, the
  clients and the release layout.
- **Interim: dev5 Postgres** (`127.0.0.1:5405` on beauty) for the current
  figure builds and for what Parquet does not have yet (listed in
  `docs/DB_ACCESS.md`).
- Do not use dev3/dev4 (empty), the old `dev.omnipathdb.org` API, or
  public hostnames such as `dev3.omnipathdb.org:5403`.
- Backends live on beauty and listen on loopback only: work there or
  through an SSH tunnel.
- beauty's instances use lock files (`~/instances/lockctl status <inst>`).
  Reading needs no claim. Never change an instance (build, restart, write)
  without a claim, and as an agent ask your user before claiming. Never
  write into another team's data directory (e.g. metsig1).

## Rules

- Every file in `data/derived/` is written by a script in `analyses/`. Do
  not commit outputs that no script in the repository produces.
- Figures and tables read only from `data/derived/`, `data/vendored/` and
  OmniPath data. Paper deliverables come out of `figures/` and `tables/`.
- Name folders, ids, functions and files by content, never by figure or
  table number.
- No personal absolute paths. Analyses resolve paths from the repository
  root (`here::here()` in R, `Path(__file__)` in Python); the pipeline runs
  from the repository root.
- Generated files go to `out/`, `output/` or `build/`, not into tracked
  folders.
- Deleting files, analyses or figures needs the maintainer's approval;
  show what a file is before proposing to delete it. Record structural
  changes (moves, removals, renames) with their commits in
  `RESTRUCTURING.md`.
- Coding style: `docs/CODING_STYLE.md` (saezverse R and Python guidelines
  plus repository rules).

## Running

On beauty:

```bash
. ~/.config/metabo-figures/env.sh                       # connection settings
./rebuild.sh --dry-run                                  # what would be built
./rebuild.sh metalinks-versions                         # one figure or table
tail -f build/logs/latest.log
Rscript analyses/cancer-cell-lines/run_analysis.R --list
```

During development, load the package from the checkout instead of the
installed version:

```bash
Rscript -e 'pkgload::load_all(".", export_all = FALSE, quiet = TRUE); source("figures/metalinks-versions/build.R")'
```

`rebuild.R` loads the installed `metabo.figures`; the shared `omnipath`
account on beauty has an old version, so install a branch into a private
library (`R CMD INSTALL --library=.Rlib .`, then `R_LIBS=$PWD/.Rlib`).
Rendering R Markdown needs pandoc (`nix-shell -p pandoc` on beauty);
captions need `xelatex`.

Python: `uv run --project python ...`; see `python/README.md`.

## Checks

- R package tests: `Rscript -e 'devtools::test()'`. Several tests fail on
  `main` already (list in `RESTRUCTURING.md`, Known issues); do not add
  new failures.
- Lint: `lintr::lint_package()` for R, `uv run --project python ruff check
  python analyses` for Python.
- After changing `R/`: regenerate `man/` and `NAMESPACE` with
  `devtools::document()`, or edit them by hand in roxygen's format.

## Git

- Trunk-based development on `main`; short-lived branches for larger
  changes.
- Stage files explicitly and check `git status` before committing; local
  files (personal notes, `TASKS.md`, data you may not redistribute) must
  not end up in commits. Personal files can be listed in
  `.git/info/exclude`.
- Commit messages: a short imperative subject that says what changed,
  and a body that says why.
