# metabo-usecases

Reproducible pipeline that regenerates every figure and table for the
OmniPath Metabo manuscript — main text and supplementary — from a pinned
OmniPath Postgres snapshot.

The repository is structured as an R package (`metabo.figures`) with
adjacent figure/table workspace directories, Python helpers under
`python/`, TikZ + LaTeX sources under `tikz/` and `tex/`, and a single
top-level rebuild entrypoint.

## Quick start

```sh
ssh -p 2323 omnipath@omnipathdb.org   # development happens on beauty
cd ~/dev/metabo-usecases
./rebuild.sh --dry-run                 # print the work plan
./rebuild.sh --png --bundle --check    # full rebuild with checks
tail -f logs/latest.log                # follow the unified pipeline log
```

See the full quickstart and contracts in the saezverse spec:
- Spec: `saezverse/ai/specifications/omnipath-metabo-figures/specs/001-figures-pipeline/spec.md`
- Plan: `…/plan.md`
- Quickstart: `…/quickstart.md`
- Contracts: `…/contracts/`

## Deliverables

- **Figure 1** — workflow-architecture diagram (TikZ) + quantitative
  database-content panels
- **Figure 2** — webapp screenshots (manual asset)
- **Figure 3** — MetalinksDB v1 vs v2 comparison
- **Figure 4** — old vs new COSMOS PKN comparison + included schematic graphics
- **Figure 5** — lung-cancer use case (refactored from `omnipath_metabo_case1`)
- **Tables** — ID resolving (Methods), RaMP comparison (Methods),
  supplementary tables
- Optional **manuscript bundle**: a single PDF concatenating every artifact
  in submission order

## Conventions

- **R code style**: follows `human/guidelines/r-coding-style.md` in the
  saezverse repository (roxygen2 docstring order, < 80-char lines,
  hanging closing parens, pipes over intermediate variables, `snake_case`
  with no `get_` prefix, `logger::log_*` for output, `@importFrom` per
  function, single implicit return).
- **Trunk-based development on `main`**; the `UseCase2` branch is a
  backup for future case studies and is not pulled into `main`.
- **Logging**: every R, Python, bash, and captured xelatex line writes to
  one file pointed at by `METABO_FIGURES_LOG`; line format in
  `contracts/log-format.md`.
- **Reproducibility**: every artifact emits a provenance sidecar; per-build
  manifests and snapshot identifiers are written to `manifests/`.

## License

BSD 3-Clause — see `LICENSE`.
