# Analysis plan

Plan for the cancer cell lines use case from the analysis-ready data set to
the metabolomics differential analysis of the chosen mutation contrast(s).
Decisions already taken are in `decision_log.md`, open points in
`open_questions.md`, data caveats in `caveats.md`. Approved 2026-10-28.

## Overview

| Step | File | Phase | Ends with |
|---|---|---|---|
| p01–p03 | `preprocessing/01_import_shorthouse.Rmd`, `02_cell_lines_mutations.Rmd`, `03_feature_ids.Rmd` | done | — |
| p04 | `preprocessing/04_depmap_status.Rmd` | 2→3 | mutation and copy-number tables |
| p05 | `preprocessing/05_analysis_dataset.Rmd` | 2→3 | one versioned analysis-ready data set |
| a01–a06 | `audit/*.Rmd` | 1b, done | decisions D7–D13 |
| 02 | `analysis/02_mutation_landscape.Rmd` | 3 | description of mutations in the metabolomics lines |
| 03 | `analysis/03_contrast_selection.Rmd` | 3 | ranked candidate contrasts → **Checkpoint 2** |
| 04 | `analysis/04_metabolomics_qc.Rmd` | 4 | QC for the chosen design |
| 05 | `analysis/05_metabolomics_differential.Rmd` | 4 | full result tables → **Checkpoint 3** |

Analysis steps live in `analysis/` (the previous `01`–`04` still occupy the
top level). All steps run through `run_analysis.R`; outputs go to
`data/derived/cancer-cell-lines/shorthouse_v1/` (data handed on) and
`analyses/cancer-cell-lines/output/<step>/` (tables, figures as PDF + PNG,
rendered reports). Thresholds and parameters go into `config/config.yml`.
Every Rmd has purpose, inputs, methods text, commented code, previews of the
data, interpretation and `sessionInfo()`.

## p04 – DepMap mutation and copy-number status

- **Inputs:** `02_cell_lines.tsv`; DepMap 24Q4 `Model.csv`, `OmicsSomaticMutations.csv`, `OmicsAbsoluteCNGene.csv` (local, `config.yml: depmap`).
- **Method:** match lines to DepMap models by RRID (report unmatched lines;
  169 of 174 matched in audit/04). For each line and gene, record whether it
  carries any protein-changing variant, a variant flagged by DepMap as hotspot,
  likely loss of function, Hess driver, or high-impact oncogene/tumour
  suppressor variant, and the protein changes. From copy number, record
  homozygous deletion (absolute copy number 0) and copy number per gene.
  The engineered line NCIH460DNP53 gets no status.
- **Outputs:** `04_mutation_status.tsv` (long: line, gene, flags, protein changes),
  `04_copy_number.tsv` (long: line, gene, copy number, deleted), manifest.
- **No decisions here;** which flag defines "mutated" is chosen per gene in step 03.

## p05 – Analysis-ready data set

- **Inputs:** `01_lines.rds`, `02_cell_lines.tsv`, `03_features.rds`, `04_*.tsv`,
  batch per injection (authors' GitHub file, local).
- **Method:** apply the analysis exclusions (`analysis.exclude_lines`: 184A1)
  and remove lines without tissue or without DepMap match from the modelling
  set (kept with a flag); log2 after half-minimum imputation (D7) with a flag per
  ion for imputed values; majority batch per line (D8); join feature
  annotation, the ratio of between- to within-line variance per ion (reported
  alongside results, D13), and the per-line mutation/copy-number tables.
- **Output:** `shorthouse_analysis_v1.rds` (matrix lines × ions, line table, ion
  table, mutation and copy-number tables) + manifest. All analysis steps read
  only this file.

## 02 – Mutation landscape (Phase 3)

Restricted to lines in the analysis data set; repeated for the intersection
with other modalities once these exist.

1. **Frequencies** of mutations and homozygous deletions overall and per
   tissue (oncoprint-style overview), for the genes tested by Shorthouse et
   al. (EV4 list, 658 genes) plus MTAP. Both definitions (any protein-changing,
   likely driver) shown side by side.
2. **Tissue specificity:** Fisher's exact test per gene × tissue (tissue vs
   all other lines), Benjamini–Hochberg (BH) adjusted, for genes mutated in at
   least `min_mutant_total` lines.
3. **Co-occurrence and mutual exclusivity** between candidate genes, accounting
   for tissue with the Cochran–Mantel–Haenszel test stratified by tissue (2 × 2
   tables per tissue), BH adjusted, with the common odds ratio as effect size.
   Chosen because co-occurrence pooled across tissues mostly reflects shared
   tissue frequencies; a tissue-adjusted logistic model is the alternative if
   strata are too sparse.

## 03 – Contrast selection (Phase 3) → Checkpoint 2

4. **Effective n per candidate gene:** number of mutant and wild-type lines in
   tissues where both groups have at least `min_per_group_per_tissue` lines,
   and the list of contributing tissues; minimum sizes in `config.yml`.
5. **Pathway-level aggregation** (e.g. RTK/RAS) only if you approve a curated,
   cited gene-set source; otherwise skipped.
6. **Positive controls:** KEAP1/NFE2L2 (glutathione ions measured
   unambiguously), STK11; MTAP deletion has genotypes in DepMap but its
   readout (MTA) is not measured; IDH has hardly any mutant lines and an
   ambiguous readout. Report availability and n.
7. **Ranking** of candidate contrasts by effective n, number of contributing
   tissues, co-mutation structure (covariates chosen by likely metabolic
   relevance, e.g. STK11/KEAP1 for KRAS in lung), balance across batches,
   biological relevance and positive-control value. For each candidate: the
   design formula `~ batch + tissue + covariates + mutation`, its rank, and
   which tissues inform the mutation coefficient.
8. **Mutation definition per gene:** oncogenes by DepMap's hotspot / high-impact
   oncogene flags, tumour suppressors by any protein-changing variant, likely
   loss of function or homozygous deletion, using DepMap's own oncogene and
   tumour-suppressor annotation (decision point, see open questions).

The contrast is chosen on genotype and design grounds **before** looking at
metabolomics results for it. An exploratory screen of metabolic separation
across many candidates is not planned; if done, it is labelled as exploratory
and its effect on the statistics of the selected contrast is recorded.

**Checkpoint 2:** ranked candidates, recommendation and design formula; you
choose.

## 04 – Metabolomics QC for the chosen design (Phase 4)

PCA of the modelling set coloured by tissue, batch, mutation status and
covariates; variance of leading components explained by each factor (as in
audit/02); check that lineage or batch does not align with the mutation groups;
overview of group sizes per tissue and batch.

## 05 – Differential analysis (Phase 4) → Checkpoint 3

- **Model:** limma on log2 intensities (D7) with the agreed design
  `~ batch + tissue + covariates + mutation`; moderated t-statistics, empirical
  Bayes variance moderation; BH over all ions per contrast. Design rank and
  estimability are checked; the tissues and number of lines informing the
  mutation coefficient are reported. MetaProViz differential analysis is used
  instead if it supports this design (to be checked in the installed version).
- **Missing values:** none beyond the imputed values <= 0 (D7); ions with
  imputed values are flagged.
- **Sensitivity analyses:** leave-one-tissue-out; covariate-adjusted vs
  exclusion-based design (the previous approach); without ions that contain
  imputed values; log2 on ions with only positive values.
- **Positive control:** where applicable (e.g. KEAP1/NFE2L2 and glutathione
  ions), checked before interpreting the main contrast.
- **Outputs:** full result table for all ions (logFC, t, p, adjusted p,
  imputation flag, replicate signal ratio, annotation), volcano and MA plots,
  sensitivity comparisons. Full signed statistics are kept for later
  footprint and network methods.

**Checkpoint 3:** results and caveats.

## Later phases (not started without approval)

Same contrast and design in other modalities, restricted to shared lines;
footprint methods (decoupleR + CollecTRI for TF activities, kinase–substrate
prior knowledge for kinase activities, OCEAN for enzyme activities) on full
statistic vectors; COSMOS with signed upstream activities and downstream
metabolites, optionally extended with MetalinksDB.

## Housekeeping

- Comments above the code in `preprocessing/03_feature_ids.Rmd`.
- Previous `01`–`04` and Figure 6 stay until step 05 replaces the differential
  analysis; then Figure 6 is switched to the new results (separate decision).
- Reproducibility: versions recorded in manifests now; Docker image as agreed
  for the repository.
