# AGENTS.md

Notes for coding agents (and humans) working in this repository.

The repository is being restructured on the `restructure` branch; the
layout below reflects that branch. The tag `preprint_v0` marks the state
used for the first preprint.

## Layout

| Path | Contents |
|---|---|
| `analyses/<name>/` | Upstream analyses (Rmd, notebooks) that produce the inputs for figures. QC outputs go to `analyses/<name>/output/` (gitignored). |
| `data/raw/<name>/` | External inputs, with a README stating where they come from. |
| `data/derived/<name>/` | Files written by the scripts in `analyses/<name>/`; figures read from here. |
| `figures/`, `tables/` | One folder per figure/table: `build.R`, `caption.tex`, `README.md`. |
| `R/`, `man/`, `tests/`, `DESCRIPTION` | R package `metabo.figures`: shared data loaders, plots, styles, provenance. |
| `python/`, `tex/`, `lib/` | Python, LaTeX and bash helpers used by the pipeline. |
| `inst/extdata/` | Vendored snapshots and manually made assets. |
| `docs/` | `CONFIGURATION.md` (connection config), `DB_ACCESS.md` (database access on beauty). |
| `omnipath_metabo_case1/` | Legacy folder: COSMOS-PKN extension and azelate notebooks, still to be moved. |

## Rules

- Every file in `data/derived/` must be written by a script in `analyses/`.
  Do not commit outputs that no script in the repo produces.
- Use repo-root-relative paths (`here::here()` in R). No personal absolute
  paths.
- Trunk-based development on `main`.

## Running the pipeline and reaching the database

The OmniPath Postgres instances run on `beauty` and listen on loopback only,
so run there or open an SSH tunnel (see `docs/DB_ACCESS.md`). Do not connect
to public hostnames such as `dev3.omnipathdb.org:5403`.

On beauty, a wrapper sets up the connection environment by sourcing
`~/.config/metabo-figures/env.sh`:

```bash
~/.local/bin/with-metabo-figures-env ./rebuild.sh fig04-metalinks-versions
```

For development, load the package and source one build script:

```bash
~/.local/bin/with-metabo-figures-env Rscript -e \
  'pkgload::load_all(".", export_all = FALSE, quiet = TRUE); source("figures/fig04-metalinks-versions/build.R")'
```
