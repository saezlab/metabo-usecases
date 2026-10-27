# record-coverage — Resource × record-type coverage (FR-015a)

**Mode**: pipeline (FR-015a)

**Deployment**: `dev3` (default for `record-coverage`) + `dev4`
(Structures row only — `structural_specificity` facet bitmap is
dev4-only). See `R/utils-deployment.R::deployment_registry`.

## Shape

| Aspect | Choice |
|---|---|
| Rows | 19 record types per FR-015a (Proteins, Metabolites, Lipids, Drugs, Food compounds, Xenobiotics, Diseases, Phenotypes, Molecular classes, Metabolism (GEMs), Allosteric regulation, Transport, Ligand-receptor, Drug-target, TF-target, Signaling, Subcellular localization, Structures, Literature references). |
| Columns | Resources, rendered at 60° rotation to fit page width. |
| Cells | `\checkmark` when the per-cell count ≥ N (default N = 1); blank otherwise. The threshold N is carried into the sidecar's `parameters` block. |
| Renderer | Hand-built LaTeX `tabular` with `\rotatebox{60}{...}` column labels — gt's LaTeX backend escapes backslashes inside column labels, so a custom render is the simpler path. |

## Row mapping (R/tables-record_coverage.R::record_coverage_row_definitions)

| Record type | Kind | Facet | Values | Deployment |
|---|---|---|---|---|
| Proteins | entity | `entity_type` | `Gene:MI:0250`, `Protein:MI:0326` | dev3 |
| Metabolites | entity | `chemical_class` | `metabolite` | dev3 |
| Lipids | entity | `chemical_class` | `lipid` | dev3 |
| Drugs | entity | `chemical_class` | `drug` | dev3 |
| Food compounds | entity | `chemical_class` | `food` | dev3 |
| Xenobiotics | blank | — | (deferred — `chemical_class.xenobiotic` is empty) | — |
| Diseases | entity | `ontology_id` | `mondo` | dev3 |
| Phenotypes | entity | `ontology_id` | `hpo` | dev3 |
| Molecular classes | entity | `entity_type` | `Cv Term:OM:0012` | dev3 |
| Metabolism (GEMs) | blank | — | (GEM not yet exposed as a distinct facet) | — |
| Allosteric regulation | blank | — | (predicate vocab collapses to `Other` — FR-007d) | — |
| Transport | relation | `predicate` | `transports` | dev3 |
| Ligand-receptor | blank | — | (predicate vocab collapses to `Other`) | — |
| Drug-target | blank | — | (predicate vocab collapses to `Other`) | — |
| TF-target | blank | — | (predicate vocab collapses to `Other`) | — |
| Signaling | relation | `predicate` | `controls`, `positively_regulates`, `negatively_regulates` | dev3 |
| Subcellular localization | blank | — | (subcellular component not yet a distinct facet) | — |
| Structures | entity | `structural_specificity` | five non-`no_structure` levels | dev4 |
| Literature references | blank | — | (`facet_evidence_bitmap` not yet derived) | — |

Blank rows render as a row of empties — the spec rule is to NOT
synthesise finer classes from secondary tables (FR-007d "Other"
honesty constraint). The caption documents the build state.

## Caption draft (pre-FR-040 source)

> Coverage of 19 manuscript record types across the integrated
> resources. Cell content is a checkmark when the underlying record
> count is ≥ N (default 1, recorded in the artifact sidecar). Rows
> are mapped to combinations of cycle-001 `facet_entity_bitmap`
> (`entity_type`, `chemical_class`, `ontology_id`,
> `structural_specificity`) and `facet_relation_bitmap`
> (`predicate`); the Structures row is the only one that routes to
> `dev4`. Rows whose underlying vocabulary has not yet landed in the
> cycle-001 build (Xenobiotics, Allosteric regulation,
> Ligand-receptor, Drug-target, TF-target, Subcellular localization,
> Literature references, Metabolism (GEMs)) render as empty rows
> rather than be synthesised from secondary tables (FR-007d "Other"
> honesty constraint).

## Files

- `build.R` — orchestration: resolves dev3 + dev4, runs
  `record_coverage_long()` (one query per non-blank row), pivots to
  wide via `record_coverage_wide()`, builds the LaTeX body via
  `record_coverage_latex()`, saves PDF + CSV via
  `tables_save_latex_pdf_csv()`, composes the caption-and-table PDF
  via `tables_compose_caption()`, writes the provenance sidecar via
  `write_sidecar()`.
- `caption.tex` — FR-040a caption source.
- `out/` — generated artifacts:
  - `record-coverage.pdf` — typeset table (xelatex).
  - `record-coverage.csv` — wide-format CSV (FR-026).
  - `record-coverage-with-caption.pdf` — table + caption (FR-041).
  - `caption.txt` — plain-text caption (FR-041a).
  - `record-coverage.pdf.provenance.json` — sidecar (FR-031,
    `parameters.threshold` records N).
