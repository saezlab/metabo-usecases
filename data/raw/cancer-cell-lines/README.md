Shorthouse 2022 publication: https://pubmed.ncbi.nlm.nih.gov/36321551/

Data: 

msb202211006-sup-0002-datasetev1.xlsx	- injections/samples
msb202211006-sup-0003-datasetev2.xlsx	- ions/features
msb202211006-sup-0004-datasetev3.xlsx	- intensity values

Cellosaurus local dependency:

- `analyses/cancer-cell-lines/02_add_metadata.Rmd` expects a local file at `data/raw/cancer-cell-lines/cellosaurus.txt`.
- Download source: `https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt`
- After download, remove the first 57 lines before using the file in the workflow.
- The prepared file is intentionally kept out of Git because it is a large static local dependency.

Example preparation command:

```bash
curl -L https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt | tail -n +58 > data/raw/cancer-cell-lines/cellosaurus.txt
```

Running `analyses/cancer-cell-lines/02_add_metadata.Rmd` also creates `data/derived/cancer-cell-lines/cellosaurus_entry_lookup.RData` as a cache of the parsed Cellosaurus entries (tracked in git).
