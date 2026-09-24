# id-resolving — Identifier resolving across integrated resources

**Mode**: pipeline (FR-014, Methods table)

**Deployment**: `dev3` (default for `id-resolving` — see
`R/utils-deployment.R::deployment_registry`).

## Shape

| Aspect | Choice |
|---|---|
| Rows | Resources (sorted by total identifier count, descending). |
| Columns | Top-10 identifier types by global frequency, then `Misc` long-tail collapse, then a `Total` column. |
| Cells | Distinct identifier count for the (resource, id_type) pair. |
| Renderer | `gt` via the LaTeX backend (`gt::as_latex()` → xelatex). |

## Caption draft (pre-FR-040 source)

> Identifiers contributed by each integrated OmniPath resource, split
> by identifier type. Counts are distinct identifiers per (resource,
> identifier type) pair on `dev3`. The ten most-populated identifier
> types are shown as named columns; the long tail is collapsed into a
> single `Misc` column. The `Total` column is the row sum across all
> identifier types (including `Misc`). The grand-summary row at the
> bottom of the body is the column sum across resources.

The canonical caption consumed by the build lives in this folder as
`caption.tex` (FR-040, FR-040a).

## Files

- `build.R` — orchestration: resolves the deployment, runs
  `tbl_id_resolving_counts()`, pivots to wide via
  `tbl_id_resolving_wide()`, renders via `tbl_id_resolving_gt()`,
  saves PDF + CSV via `tables_save_pdf_csv()`, composes the
  caption-and-table PDF via `tables_compose_caption()`, writes the
  provenance sidecar via `write_sidecar()`.
- `caption.tex` — FR-040a caption source (main title + one panel
  `(a)` sub-caption since the table is a single panel).
- `out/` — generated artifacts:
  - `id-resolving.pdf` — typeset table (xelatex via gt LaTeX).
  - `id-resolving.csv` — wide-format CSV (FR-026).
  - `id-resolving-with-caption.pdf` — table + caption (FR-041).
  - `caption.txt` — plain-text caption (FR-041a).
  - `id-resolving.pdf.provenance.json` — sidecar (FR-031).

## Performance

The underlying join is the same shape FR-007a Identifiers uses (~76 s
on dev3 against the cycle-001 build). A DB-side proposal at
`saezverse/human/plans/omnipath-improvements-2026-06-identifier-source-count.md`
proposes a derived `identifier_source_count` table that would cut
this to a single sequential scan; T070 will adopt it once it lands.
