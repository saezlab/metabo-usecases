# Open questions

Questions that change the design, each with the recommended default. Answered
questions move to `decision_log.md`. State of 2026-10-28.

| # | Question | When needed | Recommended default |
|---|---|---|---|
| Q1 | Which DepMap flag defines "mutated" per gene: any protein-changing variant, or likely driver (hotspot / likely LoF / Hess driver / high-impact)? | Phase 3 (step 03) | Oncogenes: hotspot or high-impact oncogene flag; tumour suppressors: any protein-changing variant, likely LoF or homozygous deletion. Shown for both definitions before choosing. |
| Q2 | Minimum group sizes: mutant and wild-type lines per contributing tissue, and in total | Phase 3 | At least 2 mutant and 2 wild-type lines per contributing tissue, at least 8 mutant lines in total; reported for alternatives. |
| Q3 | Pathway-level aggregation of mutations (e.g. RTK/RAS) to increase n? | Phase 3 | Not by default; only with a curated, cited gene-set source you approve. |
| Q5 | Use Ortmayr et al. 2019 as a replication set for the chosen contrast? | after Checkpoint 3 | Decide once the main results exist. |
| Q6 | Install OCEAN (and check its identifier coverage) | before the integration phase | Install then. |
| Q7 | Which sources for transcriptomics, proteomics and phosphoproteomics? | before the integration phase | DepMap 24Q4 RNA expression is available; others to be identified. |
| Q9 | When to replace Figure 6 with the new results? | after Checkpoint 3 | Separate decision with the co-authors. |
| Q10 | Ask the authors for a licence for their GitHub data? | before publication | Optional; the data are used locally and cited (D11). |
