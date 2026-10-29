# Decision log

Every design decision of the cancer cell lines use case, with date, rationale
and where the evidence is. Newest last. Bug fixes are listed in
`preprocessing_changes.md`.

| # | Date | Decision | Rationale | Evidence |
|---|---|---|---|---|
| D1 | 2026-10-28 | Shorthouse et al. 2022 is the main metabolomics data set; Ortmayr et al. 2019 is only described | Largest data set, focus of the thesis work | — |
| D2 | 2026-10-28 | Refactored pre-processing in `preprocessing/`, entry point `run_analysis.R`, config `config/config.yml`, docs in `docs/`; previous `01`–`04` kept for Figure 6 until replaced | Approved layout (Checkpoint 1a) | `README.md` |
| D3 | 2026-10-28 | MCF7 removed | QC line whose 64 injections biased the data; removed by the authors as well | Shorthouse et al. Methods; `preprocessing_changes.md` (B14) |
| D4 | 2026-10-28 | Additional cultures merged (CORL88A/B, MDAMB231/B, NCIH441/B) | Same Cellosaurus ID in EV1; the authors count them as one line | `preprocessing_changes.md` (B3) |
| D5 | 2026-10-28 | NCIH460DNP53 keeps the parental accession but no parental variations | Engineered derivative (dominant-negative TP53) | `preprocessing_changes.md` (B3) |
| D6 | 2026-10-28 | Cell line accessions from EV1 (RRID, then Cellosaurus name), two manual | Uses the authors' own annotation; avoids manual mapping errors (B2, B15) | `preprocessing/02` |
| D7 | 2026-10-28 | Values <= 0 set to half the metabolite's minimum (MetaProViz `mvi_imputation()`), then log2; all ions kept; imputed ions flagged in results | Standard, citable (Wei et al. 2018); affects 0.15 % of values; example comparisons unchanged | `audit/01` |
| D8 | 2026-10-28 | Culture batch is one covariate in every differential model; lines in several batches get their majority batch | Batch explains more variance than tissue and is partly confounded with it; adjusting is the standard approach | `audit/02` |
| D9 | 2026-10-28 | Tissue: rule-based labels renamed to the 11 tissues of Shorthouse et al., 8 lines overridden with documented reasons | Same labels as the authors, independent of their unlicensed repository | `audit/03`, `config.yml` |
| D10 | 2026-10-28 | Mutation and copy-number status from DepMap 24Q4 (matched by RRID); Cellosaurus for identity and as cross-check | Cellosaurus misses most co-mutations and all MTAP deletions; the authors' CCLE calls contain artefacts | `audit/04` |
| D11 | 2026-10-28 | Batch and tissue labels of the authors' GitHub repository used locally and cited, not redistributed | Needed for D8; no licence | Checkpoint 1b, `caveats.md` |
| D12 | 2026-10-28 | 184A1 excluded from all analyses | Non-malignant line | Checkpoint 1b, `config.yml` |
| D13 | 2026-10-28 | No filtering of ions by replicate reproducibility | limma's variance moderation handles noisy ions; avoids an extra data-driven choice | `audit/01` section 5 |
| D14 | 2026-10-28 | Other modalities chosen before the integration phase | Metabolomics first | Checkpoint 1b |
| D15 | 2026-10-28 | Analysis plan approved (`analysis_plan.md`) | Checkpoint Phase 2 | `analysis_plan.md` |
| D16 | 2026-10-29 | Lines without DepMap match or without tissue stay in the analysis data set with a flag but are not used in mutation models | Keeps the data complete; models need both annotations (Q4) | `open_questions.md` |
| D17 | 2026-10-29 | The three ions with duplicate annotation stay dropped (1096 ions) | As in the previous pipeline (Q8) | `preprocessing/01` |
| D18 | 2026-10-29 | Phase 4 starts with the positive control: NRF2 activation (KEAP1 or NFE2L2 altered) vs not in lung, `~ batch + KRAS + STK11 + TP53 + NRF2`; the main contrast is chosen after it | Validates the pipeline on an expected effect (glutathione ions measured unambiguously) before the main result | `analysis/03`, Checkpoint 2 |
| D19 | 2026-10-29 | "Mutated" by gene role: oncogenes (KRAS, NRAS, BRAF, PIK3CA, EGFR, NFE2L2) by DepMap likely-driver variant; tumour suppressors (TP53, STK11, KEAP1, SMARCA4, CDKN2A, PTEN, RB1) by any protein-changing variant or homozygous deletion | Likely-driver flag misses many KEAP1 variants; any-variant definition inflates passenger-prone genes (Q1) | `analysis/03`, `config.yml` |
| D20 | 2026-10-29 | Minimum group sizes: >= 2 mutant and >= 2 wild-type lines per contributing tissue, >= 8 effective mutant lines (Q2) | Defaults confirmed | `config.yml` |
| D21 | 2026-10-29 | KEAP1 and NFE2L2 combined as "NRF2 activation" for the positive control (Q3) | Both activate NRF2 (KEAP1 loss, NFE2L2 gain) | `config.yml` |
| D22 | 2026-10-29 | Pre-specified contrasts: NRF2 in lung (positive control), KRAS (tissue-blocked), MTAP deletion (tissue-blocked), STK11 in lung; plus an exploratory screen of all genes meeting D20 with one standard model `~ batch + tissue + TP53 + gene`, BH within each contrast, labelled hypothesis-generating and compared with Dataset EV4 | Many contrasts are of interest; fixing the main contrasts before any result keeps them confirmatory, the screen is declared as exploratory (Q11) | Checkpoint 2 |

| D23 | 2026-10-29 | The NRF2 group (D21) and the 9p21 reading of the MTAP contrast are kept; D19 is kept despite passenger variants in hypermutated lines, to be checked once for the final contrasts | NRF2: verified references, 20 OmniPath Metabo resources (6 directed KEAP1 ⊣ NFE2L2, none stimulatory), shared KEGG/Reactome/WikiPathways pathways; DepMap mutual exclusivity not shown (driven by hypermutated models). 9p21: 97 % of MTAP-deleted DepMap models are CDKN2A-deleted, verified references | `audit/07`, caveat 21 |
| D24 | 2026-10-29 | Lines whose DepMap model has no mutation calls are not used in mutation models (NCIH1755) | The model ACH-000282 has no rows in the 24Q4 mutation file, so the line read as wild type for every gene although Cellosaurus and the CCLE calls report TP53, KEAP1 and STK11 variants | `audit/08`, `preprocessing/05` |
| D25 | 2026-10-29 | Genes defined by deletion alone (MTAP) are evaluated only on lines with DepMap copy number; for tumour suppressors, lines without copy number count as not deleted | 30 lines in models lack copy number; MTAP status is unknown for them, while homozygous deletions of KEAP1, STK11 and TP53 occur in only 1–2 of 137 lines with copy number | `audit/08`, `config.yml` |
| D26 | 2026-10-29 | Differential analysis with limma directly (moderated t, BH over all ions per contrast); positive control passes if the pre-specified glutathione ions are higher as a set (`limma::fry`, one-sided p < 0.05) | MetaProViz `dma()` 4.99.0 fits `~0 + condition` only, without covariates; criterion fixed before the fit | `analysis/05`, `analysis_plan.md` |
