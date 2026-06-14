# Figure 4: MetaLinksDB resource comparison

This directory builds the FR-010 family of panels for Figure 4 (post-
2026-06-14 six-figure renumbering — was Figure 3).

Inputs:
- MetaLinksDB v2.0 from `custom_views.metalinksdb_relations` on the
  `dev5` integrated build.
- MetaLinksDB v1.0 from the SQLite shipped by `OmnipathR`.
- Vendored normalized baseline snapshots under
  `inst/extdata/fig04-mpi-baselines/` for external MPI resources
  (CellPhoneDB, scConnect, STITCH).

Outputs:
- `out/coverage.{pdf,svg}` (Panel A — FR-010a)
- `out/metabolite_classes.{pdf,svg}` (Panel B — FR-010b)
- `out/protein_classes.{pdf,svg}` (Panel C — FR-010c)
- `out/metalinks_overview.{pdf,svg}` (Panel D — FR-010d, MetaLinksDB
  2.0 grouped-bar overview)
- `out/relationship_types.{pdf,svg}` (Panel E — FR-010e, relationship
  types restricted to transport / receptor / interaction)
- `out/fig04-metalinks-versions.{pdf,svg}` — five-panel composite
- `out/fig04-metalinks-versions.pdf.provenance.json` — provenance
  sidecar carrying the dev5 build_id (FR-032) and per-resource
  evidence-dimension availability (FR-010g).
- `out/supplementary/4A..4E/*.csv` — per-panel raw-interaction CSV
  exports (FR-010h).
- `out/supplementary/legend.csv` — central legend for every
  supplementary CSV (FR-010i).
