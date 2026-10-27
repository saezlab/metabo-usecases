# RaMP ambiguity investigation

Exploratory analysis (Jonathan Schaul, June 2026) of ambiguous RaMP ID
mappings in the OmniPath chemical resolver: classifies ambiguous RaMP IDs by
InChIKey blocks and annotates the distinct-connectivity cases with
ClassyFire/ChemOnt and ChEBI metadata.

## Requirements

Not part of the figure pipeline, and not runnable from this repo alone:

- runs on `dev2`, in the `omnipath-build` virtual environment
  (needs `duckdb`, `pandas` and `pypath` with `pypath.inputs_v2.rampdb`);
- reads the resolver parquet files under
  `/home/omnipath/instances/dev2/data/chemicals/`.

Output TSVs go to `output/` next to the notebook (gitignored).
