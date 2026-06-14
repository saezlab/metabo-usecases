# Figure 3: MetaLinksDB resource comparison

This directory builds the first implementation slice for Figure 3.

Inputs:
- MetaLinksDB v2.0 from `custom_views.metalinksdb_relations` on the `dev4` deployment.
- MetaLinksDB v1.0 from the SQLite shipped by `OmnipathR`.
- Vendored normalized baseline snapshots under `inst/extdata/fig03-mpi-baselines/` for external MPI resources.

Outputs:
- `out/coverage.{pdf,svg}`
- `out/metabolite_classes.{pdf,svg}`
- `out/protein_classes.{pdf,svg}`
- `out/evidence.{pdf,svg}`
- `out/source_relationship.{pdf,svg}`
- `out/fig03-metalinks-versions.{pdf,svg}`
- `out/fig03-metalinks-versions.pdf.provenance.json`
