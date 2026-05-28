# Test fixtures

Synthetic, version-pinned inputs used by the local + CI test suite when a
real `beauty` Postgres deployment is not reachable.

## `omnipath-mini.sql`

A minimal subset of the OmniPath public schema (4 resources, ~50 rows
total) sufficient to exercise:

- `R/data/postgres.R::pg_query` (returns rows + records the result hash).
- `R/provenance/manifest.R::build_manifest_for` (introspects
  `data_source` and produces a stable Snapshot Identifier).
- The Figure 1 quantitative panels' query shape (so the rebuild driver
  can exercise the end-to-end path during CI without beauty access).

Load into a local Postgres for `R CMD check`:

```sh
createdb omnipath_test
psql omnipath_test -f tests/fixtures/omnipath-mini.sql
PGUSER=$USER PGPASSWORD= PGHOST=localhost PGPORT=5432 \
    PGDATABASE=omnipath_test \
    Rscript -e 'devtools::test()'
```

Used by `.github/workflows/rebuild.yml` to drive a smoke rebuild against
the fixture instead of beauty.

## Updating the fixture

When the upstream OmniPath schema changes, regenerate the fixture from
beauty (only the table shapes; do NOT carry real data):

```sh
ssh -p 2323 omnipath@omnipathdb.org \
    docker exec omnipath-present-dev3-postgres-1 \
    pg_dump --schema-only -U "$POSTGRES_USER" -d "$POSTGRES_DB" \
    > tests/fixtures/omnipath-schema.sql
```

Then hand-curate the synthetic rows in `omnipath-mini.sql` to cover the
panels you need.
