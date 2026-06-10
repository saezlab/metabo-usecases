# DB Access & Quick-Iteration Handover (figures pipeline)

Operational companion to **`CONFIGURATION.md`** (which covers the connection
schema + per-panel routing). This file is the *ground truth* for actually
opening a connection and iterating fast. If something here disagrees with code,
trust a live check (`docker exec … psql`).

## TL;DR

- The OmniPath Postgres instances run **on `beauty`**, each published on
  **`127.0.0.1:<port>` only** (not on a public hostname).
- Credentials are the same trivial dev creds everywhere: **user `omnipath`,
  password `omnipath`, database `omnipath`**. There is no secret to hunt for.
- Because the ports are loopback-only, you must either **run on beauty**
  (`PGHOST=localhost`) or **open an SSH tunnel** from your machine. Connecting to
  `dev3.omnipathdb.org:5403` directly will **not** work.

## Where the data is

| Deployment | Host port | Container | Build |
|-----------:|:----------|:----------|:------|
| **dev3** (default) | `5403` | `omnipath-present-dev3-postgres-1` | **gene-centric** full build, ~5.47M entities, stored `entity.label`, has `Gene:MI:0250` entity type |
| **dev4** | `5404` | `omnipath-present-dev4-postgres-1` | **protein-centric** full build, ~5.57M entities, RDKit structural-specificity + RaMP-conflict tables (no gene type) |
| dev5 | 5405 | (reserved) | the new **integrated** build (gene-centric + full chemical resolution + global protein anchoring + metabo) is being finalized in the *build* PG and will be promoted to a web slot later — **use dev3/dev4 for now** |

All three DBs are `omnipath`/`omnipath`/`omnipath`. Treat them as **read targets**
— don't run heavy DDL/`VACUUM FULL`/drops; they are shared.

## Connecting — pick one

### A. On beauty, via the container (no creds needed — simplest sanity check)
```sh
ssh -p2323 omnipath@omnipathdb.org            # off Uni network
# on Uni network instead: ssh omnipath@beauty
docker exec -it omnipath-present-dev3-postgres-1 psql -U omnipath -d omnipath
```

### B. On beauty, via the host port (for the R/Python pipeline)
`psql` is **not** on beauty's PATH — get it from nix when you need the CLI:
```sh
nix-shell -p postgresql --run \
  "PGPASSWORD=omnipath psql -h 127.0.0.1 -p 5403 -U omnipath -d omnipath"
```
Libraries connect straight to the loopback port:
```
postgresql://omnipath:omnipath@127.0.0.1:5403/omnipath      # dev3
postgresql://omnipath:omnipath@127.0.0.1:5404/omnipath      # dev4
```

### C. From your laptop, via an SSH tunnel (ports are loopback-only)
```sh
ssh -p2323 -N \
  -L 5403:127.0.0.1:5403 \
  -L 5404:127.0.0.1:5404 \
  omnipath@omnipathdb.org
# leave that running; in another shell, connect to the tunnelled local ports:
#   postgresql://omnipath:omnipath@127.0.0.1:5403/omnipath
```

## Wiring it into the pipeline (`credentials_source: env`)
```sh
export PGHOST=localhost          # on beauty, or the tunnelled localhost
export PGDATABASE=omnipath
export PGUSER=omnipath
export PGPASSWORD=omnipath
export PGPORT_DEV3=5403           # gene-centric (default)
export PGPORT_DEV4=5404           # protein-centric + RDKit
```
See `CONFIGURATION.md` for `connection.yaml`, per-panel `overrides:`, and the
`pgpass`/`file` credential sources.

## Quick-iteration flow (edit local → rsync → run on beauty → commit)

The pipeline must run **where it can reach `127.0.0.1:540x`** — i.e. on beauty
(no tunnel, lowest latency). Iterate without committing each attempt:

```sh
# 1. edit the figures code in your local clone
# 2. push it to your beauty checkout (fast, no commit):
rsync -az -e "ssh -p2323" --delete \
  --exclude .git --exclude '*.Rproj.user' --exclude 'logs/' --exclude '.venv' \
  ~/omnipath-metabo/usecases/   omnipath@omnipathdb.org:~/figures-work/usecases/
# 3. run it on beauty (DB is local there):
ssh -p2323 omnipath@omnipathdb.org \
  'cd ~/figures-work/usecases && PGHOST=localhost PGUSER=omnipath PGPASSWORD=omnipath \
   PGDATABASE=omnipath PGPORT_DEV3=5403 PGPORT_DEV4=5404 ./rebuild.sh <target>'
# 4. once it works, commit + push from your LOCAL clone (commits = checkpoints)
```

Notes:
- Create the beauty checkout once (`~/figures-work/` above is just a suggestion —
  pick a path and keep using it). `rsync` the working tree; don't `git pull` on
  beauty during iteration.
- Keep heavy downloads / DB writes off your laptop — run them on beauty.
- For a single throwaway query, recipe **A** (`docker exec … psql`) is fastest.

## Gotchas

- **Loopback-only ports.** `dev3.omnipathdb.org:5403` from outside fails — tunnel
  or run on beauty. (The hostnames in `R/utils-config.R` assume `PGHOST=localhost`
  on beauty.)
- **No `psql` on beauty's PATH** — `docker exec` into the container, or
  `nix-shell -p postgresql`.
- **dev3 vs dev4 are different builds**, not replicas: dev3 is gene-centric
  (genes + labels), dev4 is protein-centric (RDKit specificity + RaMP-conflict).
  Pick per the deployment matrix in `CONFIGURATION.md`.
- **SSH:** off the Uni network use `ssh -p2323 omnipath@omnipathdb.org`; on it,
  `ssh omnipath@beauty`.
