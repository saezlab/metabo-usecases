# Checkpoint 3: results and caveats

State of 2026-10-30, branch `cancer-cell-lines-v2` at `4712ced`. Shorthouse et
al. (2022) metabolomics (flow-injection MS) of cancer cell lines, DepMap 24Q4
genotypes, OmniPath Metabo (dev5, build `dac1e0f66ee8`) for prior knowledge.
Decisions D1–D36 in `decision_log.md`, caveats 1–24 in `caveats.md`.

## Data and model

- **167 lines** in models (of 175), **1096 ions**. Values: authors' normalised
  intensities (Lowess for injection order and confluency, mass correction,
  batch/plate-associated metabolites removed), half-minimum imputation of
  values <= 0, log2 (D7).
- Genotypes from DepMap 24Q4 matched by RRID; status by gene role (D19), copy
  number for MTAP and CDKN2A (D25, D27); lines without DepMap mutation calls
  excluded (D24).
- limma with **culture batch in every model**; tissue as blocking factor in
  pan-tissue contrasts; **lung histology** as blocking factor in lung contrasts
  (D28); co-mutations as covariates. BH over all ions per contrast.

## What stands

### 1. NRF2 activation in lung lines: the one strong mutation contrast

`analysis/05`; 24 NRF2-altered (19 KEAP1, 5 NFE2L2) vs 61 lung lines.

- **The positive control passes at the ion level.** All five pre-specified
  readouts move as expected: glutathione log2FC +0.84 (adjusted p 0.005, rank
  14 of 1096), oxidised glutathione up, glutamate down (−0.65, adjusted p
  0.064). The glutathione ions as a set are clearly higher (fry p 0.002).
  γ-glutamylcysteine, cysteinylglycine and glycine are also higher (not
  pre-specified).
- 17 ions at adjusted p < 0.05.
- The NRF2 grouping is supported by literature and by OmniPath Metabo
  (20 resources link KEAP1 and NFE2L2, 6 with directed inhibition; shared
  pathways). Mutual exclusivity in DepMap is *not* shown (audit/07).
- The KEAP1 screen in lung shows the same leading ions (FAD, 6-hydroxy-FAD,
  orotidylic acid).

### 2. Lineage contrasts: the clearest metabolic differences

`analysis/07` (exploratory, "robust" = adjusted p < 0.05 without imputed values).

| Contrast | Lines | Robust ions |
|---|---|---|
| Neuroendocrine vs non-small cell lung | 14 / 69 | 55 (cystathionine +2.9, S-adenosylhomocysteine, purine nucleotides, phospholipids) |
| Colon vs rest | 6 / 161 | 19 |
| Kidney vs rest | 10 / 157 | 13 |
| Skin vs rest | 8 / 159 | 10 |
| Adenocarcinoma vs squamous lung | 39 / 15 | 0 |
| Lung, blood, bladder, breast, ovary, pancreas vs rest | | 0–3 |

### 3. Main genotype contrasts: null at the ion level

`analysis/06`, pre-specified (D22).

| Contrast | Lines | Ions at adjusted p < 0.05 |
|---|---|---|
| KRAS, tissue-blocked | 30 / 66 | 0 (smallest adjusted p 0.51) |
| MTAP deletion (9p21), tissue-blocked | 32 / 89 | 0 (0.97) |
| STK11 in lung | 22 / 63 | 0 (0.12; top ions phospholipids) |

- Neither histology blocking (audit/10) nor the previous exclusion-based
  design (audit/09) changes this. KRAS vs EGFR is not testable (2–4 EGFR lines).
- KRAS subsets (KL, KP, KC; Skoulidis et al. 2015): no significant ion.

### 4. Exploratory screen of all genes

`analysis/06`, 21 genes, `~ batch + tissue + TP53 + gene`.

- Hits: SMARCA4 (12), CDKN2A (4, led by γ-glutamylcysteine), KEAP1 (3),
  STK11 (2). Co-mutations are not separated in this model.
- **Agreement with Shorthouse et al. (Dataset EV4):** Spearman 0.68–0.86 for
  the well-defined drivers (SMARCA4, NRAS, STK11, SMAD4, PTEN, RB1, KRAS,
  TP53, KEAP1). The direction of the published associations is reproduced.
- Screens within tissues with ≤5 altered lines are driven by single
  hypermutated lines (NCIH2342, NCIH1155) and are not interpretable at the
  gene level.

### 5. Enrichment on MetSigDB: rankings, not significance

`analysis/08`.

- Ions are linked to MetSigDB sets through OmniPath Metabo from the original
  identifiers (D30–D35). 699 of 1096 ions reach a KEGG, Reactome or
  WikiPathways set, against 529 with KEGG IDs alone. 519 ions are unambiguous.
- camera (correlation estimated, non-directional and directional): **no set
  at FDR < 0.05 in any contrast** (smallest 0.087, neuroendocrine vs NSCLC).
- **The glutathione sets rank near the top for NRF2** (WikiPathways WP100 4th
  of 139, Reactome 14th of 638), but are not significant after correction.
  camera and ULM rank sets alike (Spearman 0.87–0.98).
- The ions are truly correlated (audit/13): there is no technical factor, but
  there is one strong, probably biological axis of variation between lines.
  Sets dilute single-ion signals: KEGG glutathione metabolism links 24 ions,
  of which 6 change.

## Candidate final contrasts

1. **NRF2 activation in lung lines** (KEAP1/NFE2L2), the positive control and
   the only mutation contrast with a clear, interpretable metabolic signal.
2. **Neuroendocrine vs non-small cell lung lines**, the strongest lineage
   contrast and biologically distinct.
3. **Colon, kidney and skin vs the other lines.** Lineage identity; small
   groups (6–10 lines).

Optional, exploratory only: SMARCA4 (screen hits, not yet adjusted for its
co-mutations in lung).

## Open caveats

- **Mutation definitions:** "any variant" for tumour suppressors counts
  passengers in hypermutated lines (caveat 21); the oncogene rule misses EGFR
  exon 19 deletions and counts truncating oncogene variants (caveat 22, kept
  by D29).
- **Genotype source:** DepMap stocks, not the measured stocks (STR-verified
  identity only; audit/08).
- **Ion annotation:** 577 ions have several candidate metabolites; set
  results carry this ambiguity (share of unambiguous ions reported per set).
- **Prior knowledge:** OmniPath Metabo accessed by SQL on dev5, interim
  until MetSigDB is published as Parquet (caveat 23); eight database issues
  to report (caveat 24).
- **Normalisation:** only the authors' normalised data are available
  (caveat 2).

## Decisions needed

1. **Final contrasts:** confirm 1–3 above (and whether SMARCA4 is kept as
   exploratory).
2. **Next step:** crystallise the workflow around the final contrasts, or
   first extend to other modalities or footprint methods (later phases of the
   plan).
3. **Figure 6** of the manuscript: switch to the new results (separate
   decision).
4. **Report** the database and client issues to the OmniPath Metabo
   developers.
