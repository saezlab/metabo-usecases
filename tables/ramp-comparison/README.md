# ramp-comparison

RaMP identifiers that map to several distinct InChIKeys, grouped by the
reason of the conflict (stereo, specificity, tautomer, similar, unrelated):
RaMP ids and conflict pairs per reason, their share, and two examples each.

**Status (2026-11-24): fails.** It reads `metabo_ramp_inchikey_conflict`, a
table computed with the RDKit Postgres cartridge; neither dev5 nor the
Parquet release has it (`docs/DB_ACCESS.md`). The caption still describes the old dev4 build. When it is
updated, the conflict classification goes into `analyses/` (e.g. RDKit in
Python on the RaMP identifiers of the Parquet release) and this folder only
formats the result (see `AGENTS.md`).

## Inputs

`metabo_ramp_inchikey_conflict` through `R/tables-ramp_comparison.R`.

## Outputs (`out/`, gitignored)

| File | What |
|---|---|
| `ramp-comparison.pdf`, `.csv` | The table (gt, LaTeX) and its data |
| `ramp-comparison-with-caption.pdf`, `caption.txt` | With typeset caption; plain-text caption |
| `ramp-comparison.pdf.provenance.json` | Provenance sidecar |

## Run

```sh
./rebuild.sh ramp-comparison
```

Needs `xelatex` on the PATH.
