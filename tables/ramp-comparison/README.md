# ramp-comparison — RaMP InChIKey conflicts (structure-based identifier comparison)

**Mode**: pipeline (FR-015, Methods table)

**Deployment**: `dev4` (default for `ramp-comparison` — see
`R/utils-deployment.R::deployment_registry`). The
`metabo_ramp_inchikey_conflict` table is dev4-only.

## Shape

| Aspect | Choice |
|---|---|
| Rows | Conflict reason (`stereo`, `specificity`, `tautomer`, `similar`, `unrelated`); rows sorted by conflict-pair count descending. |
| Columns | Distinct RaMP ids · Conflict pairs · Share of RaMP ids (%) · Example RaMP ids. |
| Cells | Aggregated counts plus two example RaMP ids per reason (`array_to_string` of the sorted `LIMIT 2` subquery). |
| Grand-summary row | Column sums for `n_ramp_ids` and `n_conflict_pairs`. |
| Renderer | `gt` via the LaTeX backend (`gt::as_latex()` → xelatex). |

## Caption draft (pre-FR-040 source)

> RaMP identifiers that map to multiple distinct InChIKeys, classified
> by reason. Classification is computed in-database by the cycle-001
> RDKit Postgres cartridge and persisted in
> `metabo_ramp_inchikey_conflict` on `dev4`. Each row of the source
> table is one conflicting pair `(ramp_id, inchikey_a, inchikey_b,
> conflict_reason)`; this Methods table aggregates per reason and
> records distinct RaMP-id count, conflict-pair count, share of the
> total conflicting-RaMP population, and two representative example
> RaMP identifiers. The grand-summary row at the bottom carries the
> column sums.

The canonical caption consumed by the build lives in this folder as
`caption.tex` (FR-040, FR-040a).

## Files

- `build.R` — orchestration: resolves dev4, runs
  `tbl_ramp_comparison_summary()`, renders via
  `tbl_ramp_comparison_gt()`, saves PDF + CSV via
  `tables_save_pdf_csv()`, composes the caption-and-table PDF via
  `tables_compose_caption()`, writes the provenance sidecar via
  `write_sidecar()`.
- `caption.tex` — FR-040a caption source.
- `out/` — generated artifacts:
  - `ramp-comparison.pdf` — typeset table (xelatex via gt LaTeX).
  - `ramp-comparison.csv` — wide-format CSV (FR-026).
  - `ramp-comparison-with-caption.pdf` — table + caption (FR-041).
  - `caption.txt` — plain-text caption (FR-041a).
  - `ramp-comparison.pdf.provenance.json` — sidecar (FR-031).

## Performance

Cheap — the per-reason aggregation runs in tens of ms against the
~24 K-row `metabo_ramp_inchikey_conflict` table on dev4; no derived
table needed.
