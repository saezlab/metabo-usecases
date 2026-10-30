# Cancer cell lines use case: what was planned, what we found, and why we diverged

Written 2026-10-30 at the end of the metabolomics phase on branch
`cancer-cell-lines-v2` (state: commit `ef92bad` plus this document). Meant to
be read on its own after a few weeks away. Details live in the documents
listed at the end; decision numbers (D…) refer to `decision_log.md`, caveat
numbers to `caveats.md`.

## In short

- The **metabolomics of 167 cancer cell lines** (Shorthouse et al. 2022) is
  cleaned, re-processed and verified, with DepMap genotypes and a reproducible
  pipeline (`run_analysis.R`).
- **One mutation contrast works: NRF2 activation in lung lines** (KEAP1 or
  NFE2L2 altered). It recovers the known rise in glutathione and serves as the
  positive control.
- **The planned main mutation contrasts (KRAS, MTAP deletion, STK11) show
  nothing at the ion level.** The clearest metabolic differences are between
  **lineages**: neuroendocrine vs non-small cell lung lines, and colon, kidney
  and skin vs the rest.
- **Classical pathway enrichment gives rankings, not significant sets.**
  Metabolite sets mix substrates and products, ions are truly correlated, and
  ions often have several candidate metabolites.
- The way forward is **reaction-aware analysis** (signed enzyme footprints)
  and **other omics** (transcriptomics, proteomics). That is a larger,
  exploratory direction and moves to a new branch (`handoff_multiomics.md`).

## Where we started

**The first preprint** (tag `preprint_v0`, Figure 6 in `manuscript.yaml`)
showed KRAS and EGFR volcano plots for lung lines (old
`04_differential_analysis.Rmd`). Lines were filtered by co-mutations (e.g.
"TP53 mutant, EGFR wild type, CDKN2A wild type, STK11 wild type"). The
mutation status came from Cellosaurus, and the story relied on a single
metabolite of interest (azelate, panels C–D). The old pipeline had bugs: a
positional feature-annotation join, Cellosaurus parsing that counted wild
type and deletions as mutations, and cwd-relative paths. These were fixed on
the `restructure` branch, which has been merged into `main`.

**The plan for this phase** (planning prompt, 2026-10-27, local only):

1. clean and verify the Shorthouse metabolomics;
2. choose one mutation-based contrast **on genotype and design grounds,
   before looking at metabolomics results**, with tissue as blocking factor;
3. run it in metabolomics and later in the other omics;
4. integrate with prior-knowledge network methods (COSMOS, footprints).

Checkpoints after each phase.

## What was done

| Step | File | Question | Key result |
|---|---|---|---|
| p01–p05 | `preprocessing/01`–`05` | Import, cell line identity, ion IDs, DepMap status, analysis data set | 175 lines, 1096 ions; 167 lines in models |
| p06 | `preprocessing/06_metsigdb_sets.Rmd` | Link ions to MetSigDB sets via OmniPath Metabo | 699 ions in pathway sets (529 via KEGG IDs alone) |
| m02 | `analysis/02_mutation_landscape.Rmd` | Which genes are altered, where, together? | Tissue specificity matches known biology; only CDKN2A–MTAP co-occur |
| m03 | `analysis/03_contrast_selection.Rmd` | Which contrasts are feasible? | KRAS, MTAP, STK11, NRF2 ranked; chosen at Checkpoint 2 (D18–D22) |
| m04 | `analysis/04_metabolomics_qc.Rmd` | Are the designs estimable, confounded? | All estimable; batch is the largest factor |
| m05 | `analysis/05_positive_control_nrf2.Rmd` | Does the pipeline find a known effect? | Yes: glutathione up, glutamate down |
| m06 | `analysis/06_metabolomics_differential.Rmd` | Main contrasts + screen of 21 genes | Main contrasts null; screen agrees in direction with the paper |
| m07 | `analysis/07_exploratory_contrasts.Rmd` | Wider exploration (histology, lineage, subsets) | Lineage contrasts clear; small screens line-driven |
| m08 | `analysis/08_enrichment.Rmd` | Which MetSigDB sets change? | No set at FDR < 0.05; glutathione sets near the top |
| a01–a13 | `audit/` | Decision-relevant side analyses | See below |

The audits answer one question each:

- **01** scale and transform;
- **02** batch and tissue structure;
- **03** tissue labels;
- **04** mutation sources (Cellosaurus vs DepMap);
- **05** metabolite annotation;
- **06** Ortmayr data set (descriptive);
- **07** evidence for gene groupings (NRF2, 9p21);
- **08** genotype reliability;
- **09** the preprint's lung contrasts redone;
- **10** lung histology;
- **11** oncogene variant rule;
- **12** metabolite identifiers;
- **13** correlation between ions.

## Findings

1. **NRF2 activation in lung lines** (24 vs 61):
   - 17 ions at adjusted p < 0.05;
   - glutathione log2FC +0.84 (adjusted p 0.005);
   - glutathione ions as a set p 0.002;
   - glutamate down;
   - γ-glutamylcysteine, cysteinylglycine and glycine also up.

   The grouping of KEAP1 and NFE2L2 rests on literature and OmniPath Metabo
   interactions; mutual exclusivity in DepMap is not shown (audit/07).
2. **Lineage contrasts** (robust ions):
   - neuroendocrine vs non-small cell lung 55 (cystathionine, S-adenosylhomocysteine, purine nucleotides, phospholipids);
   - colon 19, kidney 13, skin 10 vs rest;
   - adenocarcinoma vs squamous 0.
3. **Main mutation contrasts null**:
   - KRAS 30/66;
   - MTAP deletion 32/89;
   - STK11 in lung 22/63.

   This holds also with histology blocking and with the preprint's
   exclusion-based design.
4. **Screen** of 21 genes:
   - hits for SMARCA4, CDKN2A, KEAP1 and STK11;
   - the direction agrees with Shorthouse et al. (Spearman 0.68–0.86 for drivers).
5. **Enrichment**: camera finds no set at FDR < 0.05. The glutathione sets
   rank high for NRF2 (WikiPathways 4th of 139). camera and ULM rank sets
   alike.
6. **Structure of the data**: after batch and tissue, one strong, probably
   biological axis remains (phosphatidylglycerols, citrate, nucleotides; 20 %
   of residual variance). Ions are correlated beyond chance; there is no
   technical normalisation gap (audit/13).

## Where and why we diverged from the plan

| Planned | What we did instead | Why | Decisions |
|---|---|---|---|
| One mutation contrast at the centre | NRF2 as positive control; lineage contrasts as main candidates | The pre-specified mutation contrasts are null; lineage gives the clear signal | D18, D22; Checkpoint 3 |
| Exclusion of lines with co-drivers (preprint) | Co-mutations as covariates, all lines kept | Filtering left 1–7 mutant lines and produced artefacts from imputed values | audit/09 |
| `~ tissue + covariates + mutation` | Plus culture batch everywhere; lung histology as blocking factor | Batch dominates the data; RB1, CDKN2A and STK11 depend on histology | D8, D28 |
| Cellosaurus mutation status | DepMap 24Q4 | Cellosaurus misses most co-mutations and all MTAP deletions | audit/04 |
| ID checks and translation with MetaProViz/RaMP | Linking through OmniPath Metabo entity resolution | The IDs of an ion are alternative candidates; the compatibility check removed valid ones | D30–D35, audit/12 |
| Classical pathway enrichment | camera non-directional and directional; ULM only as ranking | Substrates and products move in opposite directions; ions are correlated | D36, audit/13 |
| Mutation contrasts only | Broader exploration (histology, lineage, KRAS subsets, screens within tissues) | Asked for, to find what the data support | `analysis/07` |
| Footprints later, from transcriptomics | Reaction-aware metabolite analysis proposed as the next step | Fits metabolites better than pathway sets | `handoff_multiomics.md` |

## Open caveats

These are the ones that matter most later; all are listed in `caveats.md`.

- **21** "Any variant" counts passengers in hypermutated lines for tumour suppressors.
- **22** The oncogene rule misses EGFR exon 19 deletions and counts truncating variants.
- Genotypes are those of the DepMap stocks, not the measured stocks (audit/08).
- 577 of 1096 ions have several candidate metabolites; set results carry this ambiguity.
- **23** OmniPath Metabo is accessed by SQL on dev5, as an interim solution.
- **24** Eight database issues are to be reported to the OmniPath Metabo developers.
- **2** Only the authors' normalised data exist for Shorthouse.

## Possible next steps

On this branch (metabolomics only):

1. Confirm the final contrasts (Checkpoint 3).
2. Crystallise the pipeline around them: keep only what supports them,
   archive the rest.
3. Decide whether Figure 6 switches to the new results.
4. Report database and client issues (caveat 24, `op_metabo_id_handling.md`).

On a new branch (`handoff_multiomics.md`):

- signed enzyme footprints from OmniPath Metabo reactions;
- product/substrate ratios;
- explaining the strong biological axis;
- transcriptomics and proteomics for the same contrasts;
- later COSMOS.

## Where to find things

- `docs/checkpoint_3.md`: results summary with all numbers
- `docs/decision_log.md`: D1–D36
- `docs/caveats.md`: caveats 1–24
- `docs/analysis_plan.md`: approved plan (Phase 2)
- `docs/preprocessing_changes.md`: bug fixes against the thesis pipeline
- `docs/op_metabo_id_handling.md`: ID handling via OmniPath Metabo, SQL, client test
- `docs/handoff_multiomics.md`: handoff for the next branch
- `config/config.yml`: all parameters and contrast definitions
- Reports: `output/reports/*.html` (rendered on beauty, gitignored)
- `TASKS.md`: running task list (local only)
