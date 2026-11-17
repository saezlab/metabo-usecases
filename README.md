# Option 1: metabolic-task network inference for liver starvation

CORNETO-based metabolic-task analysis of a WT vs ob/ob mouse-liver
starvation time course (Morita et al., 8 timepoints: 0, 2, 4, 6, 8, 12,
16, 24h). This is the "Option 1" route from the
[corneto-metabo-mct](https://github.com/saezlab/corneto-metabo-mct)
notebook: metabolic tasks on a small hand-curated network, rather than the
full Mouse-GEM.

A metabolic task asks whether the network can carry out a function — it
opens the task's substrate inputs, requires a product output, and asks for
a steady-state flux within the reaction bounds. We build a 117-reaction
liver network supporting 20 fed- and fasted-state tasks, then close
reactions whose enzymes are low-expressed in the measured proteomics to see
which tasks each genotype can still perform, and how that changes over the
time course.

See [`OPTION1_SUMMARY.md`](OPTION1_SUMMARY.md) for the methods and findings.

## Scripts

Run in order; each writes to the matching directory under `result/`.

| Script | Does |
|---|---|
| `scripts/01_metabolic_task_network.py` | Build the 117-reaction network and the 20 tasks; find the smallest reaction set supporting all tasks jointly. |
| `scripts/02_temporal_task_feasibility.py` | Gate the network by each genotype's proteomics at every timepoint; track task feasibility across the 16 states, with a per-genotype overview (feasible throughout / never / at 24h). |

Script 02 reads `result/gene_scores/` (the proteomics gene scores) and
imports the network from script 01.

## Running

Requires Python ≥ 3.13, [uv](https://docs.astral.sh/uv/), and a GUROBI
license (academic licenses are free). The CORNETO analysis package is
vendored under `src/corneto_mcmt`.

```bash
uv run python scripts/01_metabolic_task_network.py
uv run python scripts/02_temporal_task_feasibility.py
```

## Layout

```
scripts/   01-03, the pipeline
src/       corneto_mcmt, the CORNETO analysis package
result/    outputs, one directory per script + gene_scores input
```
