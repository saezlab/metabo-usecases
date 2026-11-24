# Pipeline connection configuration

How the `metabo.figures` pipeline (`rebuild.R`, `figures/*/build.R`,
`tables/*/build.R`) finds its Postgres database. Which backend to use for
what, and how to reach it, is in `docs/DB_ACCESS.md`.

Figures and tables that move to the Parquet release read it directly (see
`docs/DB_ACCESS.md`); this configuration then no longer applies to them.

## Deployments

Known deployments (`R/utils-config.R::known_deployments()`):

| Deployment | Port | Role | State (2026-11-24) |
|---|---|---|---|
| `dev5` | 5405 | `figures-default` | Full load (October 2026); the only deployment the pipeline uses |
| `dev3` | 5403 | `figures-legacy` | Empty |
| `dev4` | 5404 | `figures-legacy` | Empty |
| `prod` | 5485 | `opt-in` | Production database of June 2026; loads only with `allow_optin = TRUE` |
| `dev2` | 5402 | `opt-in` | No longer on beauty |

All ports listen on `127.0.0.1` on beauty only. Every panel routes to
`default_deployment`; the built-in per-panel registry
(`R/utils-deployment.R::deployment_registry()`) is empty. When the new
Postgres instance comes, add it to `known_deployments()` and point
`default_deployment` at it.

## connection.yaml

The pipeline reads the first of:

1. `connection.yaml` in the repository root (gitignored), for a local
   override;
2. `~/.config/metabo-figures/connection.yaml`.

Without either, `load_connection()` stops with an error. Start from
`inst/extdata/connection.template.yaml`:

```yaml
default_deployment: dev5
credentials_source: env      # env | file | pgpass
credentials_path: null       # required when credentials_source is file
overrides: {}
```

- `default_deployment`: the deployment every panel uses. The legacy field
  `deployment:` is read as `default_deployment`.
- `overrides`: optional routing per panel, either `{ default: <name> }` or
  `{ <facet>: <name>, ... }`, e.g.

  ```yaml
  overrides:
    ramp-comparison:
      default: dev5
  ```

  Overrides win over the built-in registry.
- `credentials_source`:
  - `env`: user, password and database from `PGUSER`, `PGPASSWORD`,
    `PGDATABASE`;
  - `file`: a YAML file at `credentials_path` with `user`, `password`,
    `host`, `port`, `database`;
  - `pgpass`: `~/.pgpass` (libpq rules).

## Environment

```sh
export PGHOST=localhost          # on beauty, or the local end of a tunnel
export PGDATABASE=omnipath
export PGUSER=omnipath
export PGPASSWORD=omnipath
export PGPORT_DEV5=5405          # per-deployment port overrides
```

`PGPORT_<DEPLOYMENT>` (e.g. `PGPORT_DEV5`) overrides that deployment's
port, which helps with tunnels on other local ports. Plain `PGPORT` is used
only when no per-deployment variable is set.

On beauty, `~/.config/metabo-figures/env.sh` sets these variables and
`~/.config/metabo-figures/connection.yaml` selects dev5 with
`credentials_source: file`:

```sh
. ~/.config/metabo-figures/env.sh
./rebuild.sh <target>
```

## Logging

`rebuild.R` sets `METABO_FIGURES_LOG` to the run's log file
(`build/logs/`). Scripts run outside `rebuild.R` write to
`build/logs/orphan-<PID>.log` and warn on stderr.

## Check

```sh
Rscript -e 'metabo.figures::check_deployment()'
```

resolves `dev5` and logs one summary line for it, without opening a
connection. Pass other names to check those, e.g.
`check_deployment(c("dev5", "dev3"))`; `prod` and `dev2` are refused.
