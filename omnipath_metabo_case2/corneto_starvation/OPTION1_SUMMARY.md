# Option 1: metabolic-task curation on a small curated network

Status summary as of 2026-11-16, for picking this work up later.
Branch `corneto` (orphan, single commit) contains only this pipeline:
scripts 19-29, `result/option1_*`, inputs (`result/gene_scores`,
`result/metabolite_scores`), vendored `src/corneto_mcmt`, uv env files.
Option 2 / COSMOS-MOON / raw data live on branch `UseCase2`.

## What Option 1 is

The colleague's notebook (`saezlab/corneto-metabo-mct`,
`notebooks/model-example.ipynb`) defines two routes. Option 1 = metabolic
tasks on a SMALL hand-curated network (never the full Mouse-GEM — that is
Option 2's defining feature, per the notebook's own text). We reproduced
the tutorial, then expanded the network in three rounds with real
biochemistry, then gated it with our own proteomics per genotype and per
timepoint, and validated predictions against our own tissue metabolomics.

Dataset: Morita et al. starvation time course, WT vs ob/ob mouse liver,
8 timepoints (0,2,4,6,8,12,16,24h) -> 16 states.

## Network evolution

| Round | Script | Reactions | Tasks | Added |
|---|---|---|---|---|
| tutorial | 19, 20 | 58 | 6 / 4 | literal reproduction; starvation tasks exposed 3 structural gaps |
| 1 | 21 | 69 | 12 | ALT fixed to real transamination (+GLUD1), lumped urea synthesis, glycogen store, glycerol route (AQP9/GK/GPD1/TPI), glutamine route (GLS2), AKG_TO_OAA |
| 2 | 24 | 95 | 15 | full 5-enzyme urea cycle (CPS1/OTC/SLC25A15/ASS1/ASL/ARG1, Asp as 2nd N donor, fumarate->OAA), BCAA oxidation (lumped per Leu/Ile/Val), methionine cycle + transsulfuration (MAT1A..CDO1->taurine, GSH) |
| 3 | 27 | 117 | 20 | SAM methylation economy: real GNMT (+SARDH, glycine cleavage), creatine (GATM/GAMT), polyamines (AMD1/ODC1/SRM/SMS + ORN_SYN), betaine remethylation (CHDH/BHMT/DMGDH), NAD salvage (NAMPT/NMNAT1) + NNMT |

Naming: glycerol species renamed `GLY_c` -> `GLYR_c` in round 2 so
`GLY_c` = glycine (reaction ids EX_GLY_IN/GLY_KIN kept).

Recurring design rule learned by hitting infeasibilities: waste disposal
is an obligate co-product. Tasks consuming amino acids need EX_UREA_OUT
(+EX_ASP_IN) open; tasks consuming SAM need a homocysteine exit
(transsulfuration or remethylation) open. Closing them deadlocks the task.

## Pipeline (per expansion round; v3 = scripts 27/28/29)

1. **Base network, ungated** — joint solve across all tasks:
   `MultiSampleFBA(lambda_reg=1.0, beta_reg=0.0)`,
   `build_many(graph, reaction_bounds=task_bounds)`,
   explicit `problem.add_objective(problem.expr.edge_has_flux_OR.sum())`
   (needed in pinned CORNETO 1.0.0rc8 — notebook's own workaround),
   stage 2: fix union size, minimize total |flux|.
   Result: all 20 tasks jointly feasible, 100/117, <1s.
   Purpose: any later infeasibility is attributable to gene expression,
   not network structure.
2. **Gene expression gating (proteomics)** — single solves, no coupling:
   GPR human symbols `.capitalize()` -> mouse;
   `evaluate_gpr_expression(default_value=1)`; "low" = measured AND
   <= 0.10 quantile of that genotype's own network-genes x 8-states
   distribution (genes absent from the ~3,452-gene panel stay available
   — "not detected" != "not expressed"); unsupported GPR -> bounds (0,0)
   merged into each task's reaction_bounds;
   `MultiSampleFBA(lambda_reg=0.0).build()`, no objective, feasibility
   only. Per genotype (pooled medians): 2x20=40 solves — WT 12/20,
   ob_ob 15/20. Per state (threshold fixed per genotype, call on each
   state's own value): 16x20=320 solves -> feasibility matrix/heatmap.
3. **Per-state network reconstruction** — for each of the 16 states,
   take ITS feasible task subset + ITS closures, rerun the same
   two-stage joint solve -> 16 state-specific smallest joint networks
   (42-94 of 117) + per-task parsimonious flux.
   Joint solves are always across TASKS, never across states/genotypes.
4. **Validation** — check whether tissue-metabolome trajectories
   (independent modality) show the same tendency as predicted task
   feasibility.

Visualization: ALWAYS `plot_task_flux_gallery` (union panel + per-task
flux-labeled panels) — the notebook's own method; `plot_selected_network`
alone was explicitly rejected by the user. Figures saved under
`result/.../figures/` (not scratchpad). SVGs extracted from the gallery
HTML (`<?xml` .. `</svg>`), rendered with
`rsvg-convert --background-color=white`. Use neato layout (dot hangs).

## Findings (v3, 117 reactions)

Trajectories = 0h->24h, 1 = feasible:

```
                             WT         ob_ob
Glycogenolysis               11111000   11111111   Pygl lost at 12h in WT
Gluconeogenesis, glycerol    11111111   00000100   Aqp9
Glucose->lipid / fed state   00000000   11111100   Acly/Acaca/Fasn/Me1 low in WT
Transsulfuration->taurine    00000011   00000000   Mat1a recovers 16h (WT); Cdo1 never (ob_ob)
Creatine synthesis           00000011   00000000
Polyamine synthesis          00000011   00000000   Srm in ob_ob
Remethylation via betaine    00000011   11111111   mirror image of WT
Lactate gluconeo/ketogenesis/FA-ox: all 1, both genotypes (robust core)
```

- WT = ordered program: glycogen-supported 0-8h -> switch at 12h
  (glycogen depletion) -> methylation tasks all ON at 16-24h, every one
  traced to the single gene Mat1a (low early, recovers 12h+).
- ob_ob = no transition: permanent blocks (glycerol/Aqp9,
  transsulfuration/Cdo1, creatine, polyamine/Srm), fed-state lipogenesis
  retained to 12h, oscillating task counts (16,15,12,15,14,16,10,10);
  only the remethylation arm stays open (opposite arm to WT).
- Per-state networks: WT grows 57 (8h) -> 94 (24h) as methylation
  modules activate; ob_ob shrinks to 62 at 16-24h.

Metabolome cross-check (tissue, `result/metabolite_scores/metabolite/`):
- supports: Hypotaurine WT 5-8x ob_ob from 4h; Cystathionine depressed
  in ob_ob; Ornithine/Citrulline/Urea all lower in ob_ob;
  Guanidinoacetate WT 2-3x; liver Creatine peaks exactly at 16h in WT;
  WT GSH(red) falls 16-24h while ob_ob maintains.
- mixed: Putrescine higher in ob_ob (accumulation above blocked Srm —
  supports) but MTA also higher in ob_ob (contradicts). Unresolved.
- contested vs literature: glycerol/AQP9 (some papers report AQP9
  ELEVATED in ob/ob; our internal control — lactate/alanine routes
  sharing the same downstream machinery work fine in ob_ob — keeps the
  finding credible, but it stays an open question for domain experts).

## Caveats (state these when reporting)

- Per-state calls are single measurements, bottom-decile hard threshold,
  no replicates/significance testing. Mat1a/Cdo1/Srm are single-gene
  calls — check raw proteomics trajectories before relying on them.
- Feasibility is binary; flux values are parsimonious solutions, their
  magnitudes carry no quantitative meaning.
- 117-reaction toy model, not Mouse-GEM; absent pathways can't fail.
- ob_ob "oscillation" could be noise; orthogonal check would be OCEAN
  enzyme-activity estimation (see project plan
  `saezverse/plans/omnipath-metabo-case-study-2.md`).

## Where things are

- `scripts/19..29` — numbered in pipeline order; each docstring explains
  what changed and why (including bugs hit and fixed).
- `result/option1_sam_methylation_expansion/` — current base network
  (build_summary.txt, task_flux.csv, task_flux_gallery.html).
- `result/option1_genotype_task_feasibility_v3/` — per-genotype gating
  (gene_detection.csv shows which gene closed which reaction).
- `result/option1_temporal_task_feasibility_v3/` — 16-state gating:
  feasibility_status_matrix.csv, feasibility_heatmap.png,
  figures/networks/<state>.png (16 union networks),
  figures/networks/<state>/<task>.png (205 flux-labeled panels),
  task_flux_gallery_<state>.html.
- Older v1 (69-rxn) / v2 (95-rxn) result dirs kept for provenance.
- Solver: GUROBI (academic license). Env: `uv run`, PYTHONPATH=src.

## Not done / possible next steps

- Purine salvage/degradation module (AMP->...->Urate, all intermediates
  measured) — flagged, not built.
- Re-check Mat1a early-low-in-WT against literature (counterintuitive;
  not yet verified).
- Replicate-aware thresholding / significance for the gene calls.
- OCEAN enzyme-activity estimation as orthogonal validation.
