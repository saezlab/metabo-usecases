# Database access

How this repository reaches OmniPath data, and which backend to use for
what. Checked on beauty on 2026-11-24; backends change often, so verify
with the commands below before relying on a detail.

`docs/CONFIGURATION.md` covers the connection configuration of the
`metabo.figures` pipeline (Postgres).

## Backends

| Backend | Where | Release | Use for |
|---|---|---|---|
| **Parquet release** (default for new work) | beauty, `~/instances/metsig1/data-2026.11.19/` | `2026.11.19` | Analyses, and every figure or table that is updated |
| Public dev API | `https://omnipath-metabo-dev.schaul.click/api` | latest `2026.10.25.1` | Quick lookups from anywhere; client API mode; installing a local copy |
| **dev5 Postgres** (interim) | beauty, `127.0.0.1:5405` | dev5 load of October 2026 | The current figure/table builds and what Parquet does not have yet (see below) |
| New Postgres instance | — | — | Coming: loaded from a Parquet release by `packages/postgres` of the monorepo. Not available yet. |

The ground truth for releases, the API and the clients is the `monorepo`
branch of [saezlab/pypath](https://github.com/saezlab/pypath/tree/monorepo):
`packages/query` (R `omnipath.client`, Python `omnipath_client`, shared
SQL), `packages/api`, `packages/subsets` and `docs/versioning.md`.

Not usable: dev3 and dev4 (empty since autumn 2026), the old
`dev.omnipathdb.org` FastAPI, and public hostnames such as
`dev3.omnipathdb.org:5403`.

## Parquet release

A release pins one version of each resource. It is published as Parquet
files plus one DuckDB database that holds the whole release (the same
database the API serves).

```
~/instances/metsig1/data-2026.11.19/
  releases/2026.11.19.json              the pins: resource -> version
  resources/<resource>/<version>/       Parquet per resource (entity, relation, ...)
  serving/2026.11.19.duckdb             the release as one DuckDB database (18 GB)
  subsets/2026.11.19/metsigdb/          MetSigDB memberships
  subsets/2026.11.19/exclusion_metabolites/
```

No API serves `2026.11.19` yet; read the files on beauty. metsig1 is an
agent sandbox of the build team: **read only**, never write into it (see
[Instance locks](#instance-locks)). When the release moves, update the path
here and in the analyses' configuration.

### Tables of the serving database

| Table | One row per |
|---|---|
| `entity` | merged entity: `entity_key`, `label`, `entity_type`, `taxon`, `resources` (list), counts, `summary` (struct) |
| `entity_identifier` | identifier of an entity (`ns`, `id`, e.g. `inchikey`, `hmdb`, `chebi`, `uniprot`) |
| `entity_annotation` | annotation of an entity (`term`, `value`, quantities) |
| `entity_evidence` | source record stating an entity, with `resource` |
| `relation` | merged relation: `subject_n`, `object_n`, `predicate`, `category`, `interaction_class`, `sign`, `resources` |
| `relation_annotation`, `relation_evidence` | annotations and source records of relations |
| `relation_endpoint` | a relation under each of its ends (fast lookup by entity) |
| `ontology_edge` | `child_key` → `parent_key` |
| `metsigdb_membership` | metabolite in a set (`resource`, `set_label`, `hmdb`, `chebi`, `inchikey`, ...) |
| `exclusion_metabolite` | metabolite on the exclusion list, with its class |
| `resources`, `database_info`, `subset_info` | versions and build information |

Entities and relations are merged across resources; filter with
`list_contains(resources, 'hmdb')`. There are no per-resource tables.
`packages/query/python/README.md` of the monorepo describes the tables in
full.

### Reading it

On beauty, open the database read-only. Several sessions can read at once.

```r
con <- DBI::dbConnect(
    duckdb::duckdb(),
    "/home/omnipath/instances/metsig1/data-2026.11.19/serving/2026.11.19.duckdb",
    read_only = TRUE
)
DBI::dbGetQuery(con, "SELECT release FROM database_info")
DBI::dbGetQuery(
    con,
    "SELECT set_label, hmdb FROM metsigdb_membership WHERE resource = ? LIMIT 5",
    list("Reactome")
)
DBI::dbDisconnect(con)
```

```python
import duckdb

con = duckdb.connect(
    "/home/omnipath/instances/metsig1/data-2026.11.19/serving/2026.11.19.duckdb",
    read_only=True,
)
con.sql("SELECT release FROM database_info").fetchall()
```

To read rows by number fast, put the numbers into the SQL as a literal list
in a query of their own (`WHERE entity_n IN (1, 2, 3)`); DuckDB uses its
index only then (see `packages/query/README.md`).

### Through the clients

The R package `omnipath.client` and the Python package `omnipath_client`
(monorepo, `packages/query/`) answer the explorer's questions (search,
relations, entity details, MetSigDB) with the same SQL as the API:

- **API mode** asks the public dev API (release `2026.10.25.1` at the time
  of writing); nothing is downloaded.
- **Local mode** downloads a release into a data directory and writes
  `omnipath.duckdb` there (all resources: about 18 GB and 20 minutes; or a
  selection by topic or name).

```r
# install from a checkout of the monorepo branch
install.packages("packages/query/r", repos = NULL, type = "source")

library(omnipath.client)
api <- op_connect()                                    # API mode
op_metsigdb(api, resource = "Reactome", entity = "HMDB0000538")
op_init(c("hmdb", "chebi"))                            # local mode
con <- op_connect(op_data_directory())
op_sql(con, "SELECT count(*) FROM entity")
```

Neither client is installed on beauty yet. They are not on CRAN or PyPI:
install them from the monorepo checkout, in a separate environment for
Python (`pip install omnipath-client` installs an unrelated older package).

### What Parquet does not have yet

Compared with the tables the pipeline queries on dev5:

| Missing | Needed by | Status |
|---|---|---|
| Chemical classes (metabolite, drug, food, lipid), metabolic domains | architecture statistics, database-content, record-coverage | Postgres derivation (`classify/chemical_class.yaml`); the `summary.entityFacetHints` fields in Parquet are empty |
| Facet bitmaps | database-content, record-coverage, specificity | Postgres speed-up; replace with SQL on `relation` / `entity_annotation` |
| RDKit structural specificity, RaMP InChIKey conflicts | architecture statistics, database-content, ramp-comparison | Missing on dev5 too |
| MetaLinksDB | metalinks-versions, cosmos-pkn | Will become a Parquet subset; port verified against Postgres (`packages/subsets/docs/reports/metalinksdb-from-parquet-20261026.md`). Missing on dev5 (`custom_views.metalinksdb_relations`). |
| COSMOS PKN | cosmos-pkn | Will become a Parquet subset; prototype only |
| Resolved vs unresolved records | id-resolving | Resolution happens before publication, so every entity is resolved; needs a new definition |

## dev5 Postgres (interim)

The figure and table builds of `metabo.figures` connect to dev5
(`default_deployment: dev5`). Credentials are the shared dev credentials
`omnipath` / `omnipath`, database `omnipath`.

On beauty, set the connection environment and run a build:

```sh
. ~/.config/metabo-figures/env.sh
./rebuild.sh metalinks-versions
```

`env.sh` sets `PGHOST=localhost`, the user, password and the per-deployment
ports (`PGPORT_DEV5=5405`, ...). `~/.config/metabo-figures/connection.yaml`
selects dev5 and reads the credentials from `pg-credentials.yaml` next to
it. (The wrapper `with-metabo-figures-env` named in older notes no longer
exists.)

`psql` is not on beauty's PATH. For a quick query:

```sh
docker exec -it omnipath-present-dev5-postgres-1 psql -U omnipath -d omnipath
# or
nix-shell -p postgresql --run \
  "PGPASSWORD=omnipath psql -h 127.0.0.1 -p 5405 -U omnipath -d omnipath"
```

From a laptop, tunnel the loopback port:

```sh
ssh -N -L 5405:127.0.0.1:5405 beauty
# then postgresql://omnipath:omnipath@127.0.0.1:5405/omnipath
```

dev5 is shared and held read-only for collaborators: no DDL, `VACUUM FULL`
or drops. Do not run two builds against it in parallel (Postgres shared
memory ran out when this was tried).

## beauty

```sh
ssh -p 2323 omnipath@omnipathdb.org    # off the university network
ssh omnipath@beauty                    # on it
```

Everyone logs in as the shared `omnipath` user. Keep heavy downloads and
queries on beauty; the Parquet release and the Postgres ports are only
there.

Iterating on the figures: edit locally, copy the working tree to a checkout
on beauty with `rsync` (exclude `.git`, `build/`, `.venv`), run there,
commit from the local clone.

## Instance locks

The instances under `~/instances/` on beauty use a lock file per instance
(`~/instances/AGENTS.md`). Check before using one:

```sh
~/instances/lockctl status metsig1
~/instances/lockctl status dev5
```

Reading data does not need a claim. Anything that changes an instance
(builds, restarts, writes to its database or data directory) needs a claim
by the person doing it, and an agent must ask its user first. Never use
`--force` on someone else's claim.
