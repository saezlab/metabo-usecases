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
