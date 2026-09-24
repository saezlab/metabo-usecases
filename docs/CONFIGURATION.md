# Configuration

The pipeline reads its deployment + credentials configuration from
**either** `~/.config/metabo-figures/connection.yaml` **or** a
project-local `.env` (gitignored). Environment variables passed at
invocation time win over both.

## Deployment matrix (2026-06-09 cycle-001 + 002 handover)

The figure pipeline targets two OmniPath Postgres instances on `beauty`:

| Deployment | Port  | Role                | Used for |
|-----------:|:------|:--------------------|:---------|
| `dev3`     | 5403  | figures-default     | Default — gene-centric build with stored entity labels. |
| `dev4`     | 5404  | figures-structures  | FR-007a Structures facet, FR-007e, FR-007f, FR-015 — RDKit structural specificity + RaMP-conflict tables. |
| `dev5`     | 5405  | reserved            | Reserved for the ongoing integrated build — refused by `load_connection()`. |
| `prod`     | 5485  | opt-in              | Excluded from the active rotation; loadable only with `allow_optin = TRUE`. |
| `dev2`     | 5402  | opt-in              | Excluded from the active rotation; loadable only with `allow_optin = TRUE`. |

The per-panel routing is encoded in `R/utils-deployment.R::panel_deployment()`;
user-side `overrides:` block (see schema below) can shadow the built-in
registry when the integrated build promotes and a single deployment
collapses the split.

## connection.yaml schema

See also: `inst/extdata/connection.template.yaml`.

```yaml
default_deployment: dev3       # dev3 (default) | dev4 | dev5 (refused)
credentials_source: env        # env | file | pgpass
credentials_path: null         # required when credentials_source == "file"
overrides:                     # optional — user-side per-(panel, facet) routing
  database-content:
    structures: dev4
    panel_e: dev4
    ramp_conflict: dev4
  ramp-comparison:
    default: dev4
```

### `default_deployment`

The pinned OmniPath Postgres instance the pipeline targets when no
per-panel override applies. `dev3` is the cycle-002 gene-centric build
with stored entity labels — appropriate for ongoing development.
`dev5` is currently reserved for the integrated build and is refused by
`load_connection()`; `prod` / `dev2` are excluded from the active
rotation and require `allow_optin = TRUE` to load.

For backwards compatibility the legacy single-deployment field
`deployment:` is still read and promoted to `default_deployment`.

### `overrides`

Optional mapping from panel id → either `{ default: <name> }`
(panel-wide override) or `{ <facet>: <name>, … }` (per-facet
override). User overrides win over the built-in registry in
`R/utils-deployment.R::deployment_registry()`. Useful when:

- The integrated build promotes on `dev5` and all panels can collapse
  back to a single deployment (set `default_deployment: dev5`, remove
  the `overrides:` block).
- A specific panel needs to run against a one-off deployment for
  debugging (add a panel-level override).

### `credentials_source`

- `env` (default) — credentials read from `PGUSER`, `PGPASSWORD`,
  `PGHOST`, `PGPORT`, `PGDATABASE`. On beauty, `PGHOST=localhost` (no
  SSH tunneling needed).
- `file` — credentials read from the YAML file at `credentials_path`
  containing `user`, `password`, `host`, `port`, `database`.
- `pgpass` — credentials read from `~/.pgpass` per standard libpq
  semantics; only `host`, `port`, `database` come from `connection.yaml`.

## Required env (when `credentials_source: env`)

```sh
export PGHOST=localhost          # on beauty; otherwise the deployment hostname
export PGPORT_DEV3=5403          # gene-centric build (default)
export PGPORT_DEV4=5404          # protein-centric build with RDKit
export PGDATABASE=omnipath
export PGUSER=...
export PGPASSWORD=...
```

`PGPORT_DEV3` / `PGPORT_DEV4` / `PGPORT_DEV5` / `PGPORT_PROD` / `PGPORT_DEV2`
override the per-deployment port from the table above — useful when
the pipeline runs against tunnelled deployments. Plain `PGPORT` is
honoured only when no per-deployment override is set.

## Logging env

The pipeline writes a single unified log file per run, path communicated via:

```sh
METABO_FIGURES_LOG=<absolute-path-to-log-file>
```

Set automatically by `rebuild.R`. Scripts run outside `rebuild.R` (e.g.
during development) fall back to `build/logs/orphan-<PID>.log` with a stderr
warning. See `specs/001-figures-pipeline/contracts/log-format.md` for the
full contract.

## Verification

After configuring, verify the connections:

```sh
Rscript -e 'metabo.figures::check_deployment()'
```

The check resolves each deployment record listed in `c("dev3", "dev4")`
and logs a one-line summary per deployment. The full `build_manifest`
probe lives in `metabo.figures::build_manifest_for()` and lands during
the manifest refresh (commit 3 of the post-handover refresh, T026).
