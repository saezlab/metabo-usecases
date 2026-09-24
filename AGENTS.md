# AGENTS.md

Notes for coding agents (and humans) working in this repository.

The repository is being restructured on the `restructure` branch; the
layout below reflects that branch. The tag `preprint_v0` marks the state
used for the first preprint; `RESTRUCTURING.md` maps old names and paths
to new ones.

## Layout

| Path | Contents |
|---|---|
| `analyses/<name>/` | Upstream analyses (Rmd, notebooks) that produce the inputs for figures. QC outputs go to `analyses/<name>/output/` (gitignored). |
| `data/raw/<name>/` | External inputs, with a README stating where they come from. |
| `data/derived/<name>/` | Files written by the scripts in `analyses/<name>/`; figures read from here. |
| `figures/`, `tables/` | One folder per figure/table: `build.R`, `caption.tex`, `README.md`. |
| `manuscript.yaml` | Figure/table numbering. Name folders, ids and functions by content, never by number. |
| `R/`, `man/`, `tests/`, `DESCRIPTION` | R package `metabo.figures`: shared data loaders, plots, styles, provenance. |
| `python/`, `tex/`, `lib/` | Python, LaTeX and bash helpers used by the pipeline. |
| `data/vendored/` | Pinned third-party snapshots (COSMOS PKNs, MPI baselines), each with a source note. |
| `inst/extdata/` | Package resources: palettes, logos, JSON schema, connection template. |
| `docs/` | `CONFIGURATION.md` (connection config), `DB_ACCESS.md` (database access on beauty). |

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
~/.local/bin/with-metabo-figures-env ./rebuild.sh metalinks-versions
```

For development, load the package and source one build script:

```bash
~/.local/bin/with-metabo-figures-env Rscript -e \
  'pkgload::load_all(".", export_all = FALSE, quiet = TRUE); source("figures/metalinks-versions/build.R")'
```
