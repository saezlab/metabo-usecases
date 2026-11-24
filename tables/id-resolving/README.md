# id-resolving

Identifiers contributed by each resource, split by identifier type: the ten
most frequent types as columns, the rest as `Misc`, and a `Total` column and
row.

**Status (2026-11-24):** queries dev5 Postgres (`identifier_evidence`,
`entity_evidence_identifier`, about 75 s). To be redone: in the Parquet
release every entity is resolved before publication, so the table needs a
new definition (`docs/DB_ACCESS.md`). When it is updated, the counting goes
into `analyses/` and this folder only formats the result (see `AGENTS.md`).

## Inputs

dev5 Postgres through `R/tables-id_resolving.R` (`tbl_id_resolving_*()`).

## Outputs (`out/`, gitignored)

| File | What |
|---|---|
| `id-resolving.pdf`, `.csv` | The table (gt, LaTeX) and its data |
| `id-resolving-with-caption.pdf`, `caption.txt` | With typeset caption; plain-text caption |
| `id-resolving.pdf.provenance.json` | Provenance sidecar |

## Run

```sh
./rebuild.sh id-resolving
```

Needs `xelatex` on the PATH.
