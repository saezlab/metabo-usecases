Raw inputs of the cancer cell lines use case. Files are kept as downloaded;
do not edit them.

## Shorthouse et al. 2022

Shorthouse D. et al. (2022) Heterogeneity of the cancer cell line metabolic
landscape. Mol Syst Biol 18:e11006. <https://pubmed.ncbi.nlm.nih.gov/36321551/>
Supplementary files, downloaded 2025-08-05. The study reanalyses the
metabolomics data of Cherkaoui et al. 2022.

| File | Content |
|---|---|
| `msb202211006-sup-0001-appendix.pdf` | Appendix (supplementary methods and figures) |
| `msb202211006-sup-0002-datasetev1.xlsx` | Dataset EV1: data for sample injections |
| `msb202211006-sup-0003-datasetev2.xlsx` | Dataset EV2: metabolite annotations for each measured peak |
| `msb202211006-sup-0004-datasetev3.xlsx` | Dataset EV3: measured peak levels for each injection |
| `msb202211006-sup-0005-datasetev4.xlsx` | Dataset EV4: t-statistics for each metabolite / cancer driver gene mutation regression model |
| `msb202211006-sup-0006-datasetev5.xlsx` | Dataset EV5: t-statistics for specific KRAS mutations against metabolites |
| `msb202211006-sup-0007-datasetev6.xlsx` | Dataset EV6: log10 p-values (hypergeometric test), PROGENy pathway scores vs SMPDB pathway activity |
| `msb202211006-sup-0008-datasetev7.xlsx` | Dataset EV7: log10 p-values, DoRothEA A+B regulon activity vs SMPDB pathway activity |
| `msb202211006-sup-0009-datasetev8.xlsx` | Dataset EV8: log10 p-values, IC50 values vs SMPDB pathway activity |
| `msb202211006-sup-0010-datasetev9.xlsx` | Dataset EV9: Pearson correlation of pathway sensitivities between drug pairs |

EV1–EV3 are used by `analyses/cancer-cell-lines/`. EV4–EV9 and the appendix
were added on 2026-10-28 as reference material; unused files may be removed
later.

## Ortmayr et al. 2019 (`ortmayr_2019/`)

Ortmayr K., Dubuis S., Zampieri M. (2019) Metabolic profiling of cancer cells
reveals genome-wide crosstalk between transcriptional regulators and
metabolism. Nat Commun 10:1841. <https://doi.org/10.1038/s41467-019-09695-9>
Supplementary files (downloaded 2025-08-05). Not used by any analysis yet and
not in git (41 MB); kept locally until it is used. To set it up, download the
files below from the article's supplementary information into
`data/raw/cancer-cell-lines/ortmayr_2019/`.

| File | Content |
|---|---|
| `41467_2019_9695_MOESM2_ESM.xlsx` | Supplementary Data 1: metabolomics of 54 adherent NCI-60 cell lines. Sheets: cell line metadata (tissue, doubling time, cell size), relative metabolite abundances, their standard errors, regression p-values, tissue-specific metabolites, glucose uptake and lactate secretion rates, enzyme–metabolite distances, mouse tissue abundances |
| `41467_2019_9695_MOESM8_ESM.docx` | Titles and descriptions of Supplementary Data 1–5 (Data 2–5, e.g. TR activity profiles, were not downloaded) |

## Cellosaurus (local dependency)

- `analyses/cancer-cell-lines/02_add_metadata.Rmd` expects a local file at `data/raw/cancer-cell-lines/cellosaurus.txt`.
- Download source: `https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt`
- After download, remove the first 57 lines before using the file in the workflow.
- The prepared file is intentionally kept out of Git because it is a large static local dependency.

Example preparation command:

```bash
curl -L https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt | tail -n +58 > data/raw/cancer-cell-lines/cellosaurus.txt
```

Running `analyses/cancer-cell-lines/02_add_metadata.Rmd` also creates `data/derived/cancer-cell-lines/cellosaurus_entry_lookup.RData` as a cache of the parsed Cellosaurus entries (tracked in git).
