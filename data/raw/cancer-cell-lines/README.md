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

EV1–EV3 are the input of the pre-processing. EV4 (genes tested by the
authors) is used by the analysis and audit 04, EV5 by audit 04. EV6–EV9 and
the appendix were added on 2026-10-28 as reference material and are not used
by any script.

## Shorthouse et al. GitHub data (`shorthouse_github/`, local only)

`Cellline_mappings_shorthouse.csv`, `Metabolite_levels_shorthouse.csv`,
`Mutations_in_celllines.csv` from `Data/` of
<https://github.com/shorthouse-mrc/CellLine_Metabolomics> at commit
`cce3b15776c28f6a45efe74b3566e9e4efd56d3a` (downloaded 2026-10-28; see
`SOURCE.txt` there). The repository has no licence, so the files are not
committed. Used by audits 02, 03 and 04.

## Cherkaoui et al. 2022 primary data (`cherkaoui_2022/`)

Cherkaoui S., Durot S., Bradley J., Critchlow S., Dubuis S., Masiero M. M.,
Wegmann R., Snijder B., Othman A., Bendtsen C., Zamboni N. (2022) A functional
analysis of 180 cancer cell lines reveals conserved intrinsic metabolic
programs. Mol Syst Biol 18:e11033. <https://doi.org/10.15252/msb.202211033>

Data record: Zamboni N., Cherkaoui S. (2021) *Epithelial-mesenchymal
transition status is the main driver of intrinsic metabolism in cancer cell
lines*. ETH Zurich Research Collection.
<https://doi.org/10.3929/ethz-b-000511784>, licence CC BY 4.0. Cited as the
data source in the paper's Data Availability section (local copy:
`literature/Cherkaoui_2022_PMC9627673.xml`).

`Metabolomics_180CCL.xlsx` (MD5 `32d8509f90c192367730d6097fbf509c`) is taken
unchanged from `PrimaryAnalysis.zip` of that record, downloaded 2026-11-27 from
<https://www.research-collection.ethz.ch/bitstreams/c3737fa0-7e05-4749-b4a2-af4df64ed054/download>.
The other files of the archive (`metabolomics_180CCL.mzTab`,
`OmicsIntegration_180CCL.csv`, `PathwayScore_180CCL.csv`) are not kept.

| Sheet | Content |
|---|---|
| `injections` | 1153 injections: sample name, perturbation (cell line), confluency, injection sequence, **culture batch (1–7)**, cell line, tissue |
| `ions` | 1809 ions: m/z, quality, top annotation |
| `annotation` | candidate metabolites per ion (formula, adduct, score, HMDB and KEGG IDs) |
| `intensities` | ion × injection intensities of the primary analysis |

The Shorthouse et al. data are a reanalysis of these measurements with more
repeat cultures (EV1: sub-batches such as `Batch1A`). Where both list the same
sample, the batch labels agree; the majority batch per cell line agrees for
173 of 176 shared lines (checked 2026-11-27).

## Ortmayr et al. 2019 (`ortmayr_2019/`)

Ortmayr K., Dubuis S., Zampieri M. (2019) Metabolic profiling of cancer cells
reveals genome-wide crosstalk between transcriptional regulators and
metabolism. Nat Commun 10:1841. <https://doi.org/10.1038/s41467-019-09695-9>
Supplementary files (downloaded 2025-08-05). Used by audits 03 and 06. Not
in git (41 MB): download the files below from the article's supplementary
information into `data/raw/cancer-cell-lines/ortmayr_2019/`.

| File | Content |
|---|---|
| `41467_2019_9695_MOESM2_ESM.xlsx` | Supplementary Data 1: metabolomics of 54 adherent NCI-60 cell lines. Sheets: cell line metadata (tissue, doubling time, cell size), relative metabolite abundances, their standard errors, regression p-values, tissue-specific metabolites, glucose uptake and lactate secretion rates, enzyme–metabolite distances, mouse tissue abundances |
| `41467_2019_9695_MOESM8_ESM.docx` | Titles and descriptions of Supplementary Data 1–5 (Data 2–5, e.g. TR activity profiles, were not downloaded) |

## DepMap 24Q4 Public (`depmap_24Q4/`, local dependency)

Source of mutation and copy-number status (decision in
`analyses/cancer-cell-lines/audit/04_mutation_sources.Rmd`). DepMap 24Q4 Public,
Figshare <https://doi.org/10.25452/figshare.plus.27993248.v1>, licence CC BY 4.0.
Not in git because of its size (about 1.2 GB). Download the files below from
the Figshare record into `data/raw/cancer-cell-lines/depmap_24Q4/` and check the
MD5 sums given there:

| File | Content |
|---|---|
| `README.txt` | Description of all files of the release |
| `Model.csv` | Cell line (model) annotation, including RRID |
| `OmicsSomaticMutations.csv` | Somatic point mutations and indels (Mutect2) with annotation |
| `OmicsSomaticMutationsMatrixHotspot.csv` | Model × gene: hotspot mutations |
| `OmicsSomaticMutationsMatrixDamaging.csv` | Model × gene: likely loss-of-function mutations |
| `OmicsAbsoluteCNGene.csv` | Gene-level absolute copy number (PureCN) |
| `OmicsExpressionProteinCodingGenesTPMLogp1.csv` | Model × gene: RNA expression of protein-coding genes, log2(TPM + 1); for the NRF2 activity audit. Downloaded 2026-11-27 (Figshare file 51065489, MD5 `71794802b750ce77c422dad0720a40af`) |

## Cellosaurus (local dependency)

- `analyses/cancer-cell-lines/preprocessing/02_cell_lines_mutations.Rmd` expects a local file at `data/raw/cancer-cell-lines/cellosaurus.txt` (release 53.0 is used, see `config.yml`).
- Download source: `https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt`
- After download, remove the first 57 lines before using the file in the workflow.
- The prepared file is intentionally kept out of Git because it is a large static local dependency.

Example preparation command:

```bash
curl -L https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt | tail -n +58 > data/raw/cancer-cell-lines/cellosaurus.txt
```
