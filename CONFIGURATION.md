# Configuration

The pipeline reads its deployment + credentials configuration from
**either** `~/.config/metabo-figures/connection.yaml` **or** a
project-local `.env` (gitignored). Environment variables passed at
invocation time win over both.

## connection.yaml schema

See also: `inst/extdata/connection.template.yaml`.

```yaml
deployment: dev3                # one of: prod, dev2, dev3, dev4, dev5
credentials_source: env         # one of: env, file, pgpass
credentials_path: null          # required when credentials_source == "file"
```

### `deployment`

The pinned OmniPath Postgres instance on the `beauty` workstation. See
`saezverse/human/dev-deployments.md` for the current deployment table.
`dev5` is currently empty and is refused by `load_connection()` until
populated.

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
export PGHOST=localhost   # on beauty; otherwise the deployment hostname
export PGPORT=5403        # 5402=dev2, 5403=dev3, 5404=dev4, 5405=dev5, 5485=prod
export PGDATABASE=omnipath
export PGUSER=...
export PGPASSWORD=...
```

## Logging env

The pipeline writes a single unified log file per run, path communicated via:

```sh
METABO_FIGURES_LOG=<absolute-path-to-log-file>
```

Set automatically by `rebuild.R`. Scripts run outside `rebuild.R` (e.g.
during development) fall back to `logs/orphan-<PID>.log` with a stderr
warning. See `specs/001-figures-pipeline/contracts/log-format.md` for the
full contract.

## Verification

After configuring, verify the connection:

```sh
Rscript -e 'metabo.figures::check_deployment()'
```

The check prints the deployment name, the three per-build snapshot
identifiers, and the `partial_build` flag for each. Implemented in
Phase 2A (T014) — until then this command stubs out.
