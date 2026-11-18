# Option 1: metabolic-task curation on a small curated network

For picking this work up later. Branch `corneto` contains only this
pipeline: `scripts/01-03`, their `result/` directories, the
`result/gene_scores` input they read, vendored `src/corneto_mcmt`, and the
uv environment files. Option 2 / COSMOS-MOON / raw data live on branch
`UseCase2`.

## What Option 1 is

The colleague's notebook (`saezlab/corneto-metabo-mct`,
`notebooks/model-example.ipynb`) defines two routes. Option 1 = metabolic
tasks on a small hand-curated network (never the full Mouse-GEM — that is
Option 2's defining feature). A metabolic task asks whether the network can
carry out a function: it opens the task's substrate inputs, requires a
product output, and asks for a steady-state flux within the reaction bounds.

Dataset: Morita et al. starvation time course, WT vs ob/ob mouse liver,
8 timepoints (0, 2, 4, 6, 8, 12, 16, 24h) -> 16 states.

## The network (117 reactions, 20 tasks)

Central carbon metabolism (glycolysis/gluconeogenesis, PPP, TCA, oxidative
phosphorylation), lipogenesis, beta-oxidation and ketogenesis, glycogen
storage, the full urea cycle, BCAA catabolism, the methionine cycle with
transsulfuration and glutathione, and the SAM methylation economy (GNMT,
creatine, polyamines, betaine remethylation, NAD salvage and NNMT).

Design rule, learned from infeasibilities: waste disposal is an obligate
co-product. Tasks consuming an amino acid need EX_UREA_OUT (+EX_ASP_IN)
open; tasks consuming SAM need a homocysteine exit (transsulfuration or
remethylation) open. Closing them deadlocks the task.

## Pipeline

1. **`01_metabolic_task_network`** — build the network and find the
   smallest reaction set supporting all 20 tasks jointly:
   `MultiSampleFBA(lambda_reg=1.0)`, `build_many`, minimize the reaction
   union, then fix its size and minimize total |flux|. All 20 tasks are
   jointly feasible (100/117 reactions). This ungated baseline means any
   later infeasibility is attributable to gene expression, not structure.
2. **`02_temporal_task_feasibility`** — close a reaction when its enzyme is
   low-expressed (measured AND <= the 0.10 quantile of that genotype's
   network genes x 8 states; genes absent from the panel stay available).
   The threshold is fixed per genotype, but the call uses each state's own
   value, so task feasibility is tracked across the 16 (genotype x
   timepoint) states (feasibility matrix + heatmap). A per-genotype overview
   is derived from the real states — tasks feasible at every timepoint,
   never feasible, and feasible at the 24h endpoint — rather than from a
   pooled timepoint median, which in a time course corresponds to no real
   physiological state.

3. **`03_state_network_reconstruction`** — for each state, take its feasible
   tasks (from 02) and its reaction closures, and rerun the two-stage solve
   of step 1 on that subset. This gives a state-specific minimal network
   (WT 57-94, ob/ob 62-84 of 117) that shows how the active metabolism
   rewires over time: WT expands as fasting modules switch on, ob/ob stays
   smaller. Infeasible tasks are excluded (a joint solve including one would
   itself be infeasible). Renders a union-network figure per state and a
   per-task flux panel, with `plot_task_flux_gallery`.

Feasibility solves (step 2) are single-task (no coupling); the joint
two-stage solve (steps 1 and 3) couples across the tasks in question.

## Findings

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

- WT = ordered program: glycogen-supported 0-8h -> switch at 12h (glycogen
  depletion) -> methylation tasks all on at 16-24h, every one traced to the
  single gene Mat1a (low early, recovers 12h+).
- ob_ob = no transition: permanent blocks (glycerol/Aqp9,
  transsulfuration/Cdo1, creatine, polyamine/Srm), fed-state lipogenesis
  retained to 12h; only the remethylation arm stays open (opposite arm
  to WT).

## Caveats (state these when reporting)

- Per-state calls are single measurements, bottom-decile hard threshold,
  no replicates or significance testing. Mat1a/Cdo1/Srm are single-gene
  calls — check raw proteomics trajectories before relying on them.
- Feasibility is binary; flux values are parsimonious solutions, their
  magnitudes carry no quantitative meaning.
- 117-reaction toy model, not Mouse-GEM; absent pathways cannot fail.

## Where things are

- `scripts/01-03` — pipeline order; each docstring explains what it does.
- `result/metabolic_task_network/` — the ungated base network
  (build_summary.txt, task_flux.csv, selected_reactions.csv,
  task_flux_gallery.html).
- `result/temporal_task_feasibility/` — 16-state gating
  (feasibility_status_matrix.csv, gene_detection_by_state.csv showing which
  gene closed which reaction, feasibility_heatmap.png, and a per-genotype
  overview in build_summary.txt).
- `result/state_network_reconstruction/` — per-state minimal networks
  (selected_reactions_by_state.csv, reconstruction_summary.csv, and
  figures/<state>.{svg,png} + figures/<state>/<task>.{svg,png}).
- `result/gene_scores/` — proteomics gene scores, the input read by 02.
- Solver: GUROBI (academic license). Figures: rsvg-convert + graphviz neato.
  Run with `uv run`.

## Possible next steps

- Purine salvage/degradation module (all intermediates measured).
- Re-check Mat1a early-low-in-WT against literature (counterintuitive).
- Replicate-aware thresholding / significance for the gene calls.
- OCEAN enzyme-activity estimation as orthogonal validation.
