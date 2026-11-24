# ramp-ambiguity

Exploratory analysis (Jonathan Schaul, June 2026) of ambiguous RaMP ID
mappings in the OmniPath chemical resolver: classifies ambiguous RaMP IDs by
InChIKey blocks and annotates the cases with distinct connectivity with
ClassyFire/ChemOnt and ChEBI metadata.

**Status (2026-11-24):** kept as is; may change soon. Not part of the figure
pipeline, and it cannot run at the moment: it read the resolver Parquet
files of the `dev2` instance (`/home/omnipath/instances/dev2/data/chemicals/`),
which no longer exists, and needs `pypath` with `pypath.inputs_v2.rampdb` in
the `omnipath-build` environment. The RaMP identifiers and InChIKeys are in
the Parquet release (`entity_identifier`, namespace `ramp`;
`docs/DB_ACCESS.md`), which is the place to read them from when it is
updated.

## Files

| File | What |
|---|---|
| `ramp_ambiguity_investigation.ipynb` | The analysis |
| `output/` | TSVs written by the notebook (gitignored) |
