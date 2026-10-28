# Pre-processing changes (Shorthouse 2022)

Changes of the refactored pre-processing (`preprocessing/`, output
`data/derived/cancer-cell-lines/shorthouse_v1/`) relative to the previous
pipeline (`01_preprocess.Rmd` … `03_feature_processing.Rmd`, which stays
in place for Figure 6 until the differential analysis is replaced). Bug ids
refer to the review of 2026-10-28.

Verification: `preprocessing/verify_against_previous.R` compares both
pipelines; every difference below is attributed there, none is unexplained.

## Structure

| Previous | New | Why |
|---|---|---|
| `01_preprocess.Rmd` | `preprocessing/01_import_shorthouse.Rmd` | Keeps the injection level, joins by key, logs filters. |
| `02_add_metadata.Rmd` | `preprocessing/02_cell_lines_mutations.Rmd` | Accessions from EV1, variations from the local Cellosaurus file, long mutation table. |
| `03_feature_processing.Rmd` | `preprocessing/03_feature_ids.Rmd` | Targeted fixes only (see B6, B16); ClassyFire classes moved here from 02. |
| RData objects | `01_*.rds`, `02_*.tsv`, `03_features.rds` + `0x_manifest.json` | Versioned folder; manifests record input checksums, git commit and package versions. |
| paths in each Rmd | `config/config.yml`, `R/use_case.R` | One place for paths and parameters. |
| knit each Rmd by hand | `run_analysis.R` | Renders all or selected steps in fresh R sessions. |

## Changes to the data

| Id | Change | Effect (verified) |
|---|---|---|
| B14 | MCF7 removed, as in Shorthouse et al. 2022 (QC line, 64 injections, biased the data) | 1298 → 1234 injections; matches the authors' final data set exactly |
| B3 | CORL88A/B, MDAMB231/MDAMB231B, NCIH441/NCIH441B merged (same Cellosaurus ID in EV1; counted as one line by the authors) | 178 → 175 lines; merged values are the injection-weighted means of the previous per-name values |
| B3 | NCIH460DNP53 (engineered, dominant-negative TP53) keeps the parental accession but receives no sequence variations | 4 parental variations no longer assigned |
| B2 | NCIH441B was mapped to NCI-H322 (CVCL_1556); it is NCI-H441 (CVCL_1561, KRAS p.G12V) | resolved through the EV1 Cellosaurus name; merged with NCIH441 |
| B15 | 977 was mapped to SNU-977 (CVCL_5109, rectal); EV1 annotates it as 97-7 (CVCL_8625, bladder) | variations of SNU-977 no longer assigned |
| B1 | Cellosaurus "None_reported" (reported wild type) counted as mutated | 14 line–gene pairs, among them 13 TP53 |
| B1 | Gene deletions counted as mutated | 27 line–gene pairs (e.g. CDKN2A) now status "deletion" |
| B1 | Gene fusions entered `mutated_genes` as tokens like "ABL1 + HGNC" | 18 tokens; fusions now status "fusion" with partner gene |
| B8 | Variations from the local Cellosaurus flat file (release 53.0) instead of the live API | no differences to the API results of 2026-10-27; results tied to one release |
| B6 | Equivalent identifiers joined by ion index instead of the non-unique name (port of thesis fix METRIC 65f78bd) | ions 797 and 815 no longer exchange identifiers |
| B10 | Injections, ions and lines matched by key (`dsIdx`, `V<ionIdx>`) instead of by position | no value changes (positions happened to match) |
| B11 | `"\bNA\b"` (a backspace) replaced by `"\\bNA\\b"` | no change in tissue labels |
| B16 | addEquivalent UpSet plot saved the traverse plot | QC plot only |

Unchanged on purpose (decided later, see `TASKS.md`):

- Intensities stay on the published, normalised scale; negative values are
  kept. The transformation is chosen in the data inventory (Phase 1b).
- The three ions with duplicate annotation are still dropped
  (`drop_duplicate_annotation_ions`); the authors kept all 1099 ions.
- Tissue rules are those of the previous pipeline; the authors' 11 tissue
  labels are compared in Phase 1b.
- Equivalent identifiers with 50 or more expansions are still discarded (B12).
