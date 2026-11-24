# record-coverage

Which resource contributes which record type (proteins, metabolites,
lipids, drugs, diseases, transport, signalling, structures, ...): a
checkmark where the resource has at least N records of the type (N = 1,
recorded in the sidecar).

**Status (2026-11-24):** queries the facet bitmaps on dev5 Postgres. The
chemical-class rows (metabolites, lipids, drugs, food) and the structures
row need Postgres derivations that the Parquet release does not have
(`docs/DB_ACCESS.md`). When it is updated, the counting goes into
`analyses/` and this folder only formats the result (see `AGENTS.md`).

## Rows

Defined in `R/tables-record_coverage.R::record_coverage_row_definitions()`:
each record type maps to an entity facet (`entity_type`, `chemical_class`,
`ontology_id`, `structural_specificity`) or a relation predicate. Types
without a distinct facet in the database (xenobiotics, GEM metabolism,
allosteric regulation, ligand-receptor, drug-target, TF-target, subcellular
localisation, literature) are left empty rather than derived from other
tables.

## Outputs (`out/`, gitignored)

| File | What |
|---|---|
| `record-coverage.pdf`, `.csv` | The table (hand-built LaTeX, rotated column labels) and its data |
| `record-coverage-with-caption.pdf`, `caption.txt` | With typeset caption; plain-text caption |
| `record-coverage.pdf.provenance.json` | Provenance sidecar |

## Run

```sh
./rebuild.sh record-coverage
```

Needs `xelatex` on the PATH.
