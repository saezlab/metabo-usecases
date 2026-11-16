"""Script 23, rerun on Script 24's expanded 95-reaction network (full urea
cycle + BCAA catabolism + methionine cycle, 15 tasks instead of 12).
Mechanism is IDENTICAL to Script 23 (per-state low-expression calls against
a fixed per-genotype threshold, preserving the time trajectory) -- see
Script 23's own docstring for the full reasoning; not repeated here.

Outputs (corneto_starvation/result/option1_temporal_task_feasibility_v2/):
  gene_detection_by_state.csv, task_feasibility_by_state.csv,
  feasibility_status_matrix.csv, feasibility_heatmap.png, build_summary.txt
  -- same set as Script 23.
"""

from __future__ import annotations

import sys
from pathlib import Path

PKG_SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(PKG_SRC))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import importlib  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from corneto.methods import MultiSampleFBA  # noqa: E402

script22 = importlib.import_module("22_option1_genotype_task_feasibility")  # noqa: E402
script23 = importlib.import_module("23_option1_temporal_task_feasibility")  # noqa: E402
script24 = importlib.import_module("24_option1_amino_nitrogen_expansion")  # noqa: E402

PKG_ROOT = Path(__file__).resolve().parents[1]
OUTDIR = PKG_ROOT / "result" / "option1_temporal_task_feasibility_v2"

GENOTYPES = script23.GENOTYPES
STATE_SUFFIXES = script23.STATE_SUFFIXES
TIMEPOINT_HOURS = script23.TIMEPOINT_HOURS
STATUS_GOOD = script23.STATUS_GOOD
STATUS_CRITICAL = script23.STATUS_CRITICAL


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    graph = script24.build_full_metabolic_network()
    task_bounds_all = script24.build_task_bounds()
    task_names = list(task_bounds_all)
    print(f"Network: {graph.num_edges} reactions, {len(task_names)} tasks, "
          f"{len(GENOTYPES) * len(STATE_SUFFIXES)} states")

    network_genes = script22.network_mouse_genes(graph)
    wide_scores = pd.read_csv(script22.GENE_SCORES_PATH, index_col=0)
    low_by_state = script23.low_expression_genes_per_state(wide_scores, network_genes)

    gene_rows = []
    feasibility_rows = []
    status_matrix = pd.DataFrame(
        index=task_names,
        columns=[f"{g}_{t}" for g in GENOTYPES for t in STATE_SUFFIXES],
        dtype=bool,
    )

    for genotype in GENOTYPES:
        for timepoint in STATE_SUFFIXES:
            state = f"{genotype}_{timepoint}"
            availability = script22.reaction_availability(graph, low_by_state[state])
            for reaction_id, info in availability.items():
                gene_rows.append({"state": state, "genotype": genotype, "timepoint": timepoint,
                                   "reaction": reaction_id, **info})
            genotype_closure = {
                rid: (0.0, 0.0) for rid, info in availability.items() if not info["available"]
            }

            for task_name, bounds in task_bounds_all.items():
                merged_bounds = {**bounds, **genotype_closure}
                method = MultiSampleFBA(lambda_reg=0.0, beta_reg=0.0)
                problem = method.build(graph, reaction_bounds=merged_bounds)
                result = problem.solve(solver="GUROBI", verbosity=0)
                status = str(result.status).lower()
                feasible = status == "optimal"
                status_matrix.loc[task_name, state] = feasible
                feasibility_rows.append({
                    "genotype": genotype, "timepoint": timepoint, "state": state,
                    "task": task_name, "status": status,
                })
            n_feasible = int(status_matrix.loc[:, state].sum())
            print(f"  {state}: {n_feasible}/{len(task_names)} tasks feasible")

    pd.DataFrame(gene_rows).to_csv(OUTDIR / "gene_detection_by_state.csv", index=False)
    feasibility_df = pd.DataFrame(feasibility_rows)
    feasibility_df.to_csv(OUTDIR / "task_feasibility_by_state.csv", index=False)
    status_matrix.to_csv(OUTDIR / "feasibility_status_matrix.csv")

    # --- Heatmap: task x timepoint, one panel per genotype, status colors ---
    fig, axes = plt.subplots(1, 2, figsize=(13, 6.5), sharey=True)
    for ax, genotype in zip(axes, GENOTYPES):
        cols = [f"{genotype}_{t}" for t in STATE_SUFFIXES]
        sub = status_matrix[cols].astype(float)
        ax.imshow(sub, aspect="auto", cmap=plt.matplotlib.colors.ListedColormap([STATUS_CRITICAL, STATUS_GOOD]),
                  vmin=0, vmax=1)
        ax.set_xticks(range(len(STATE_SUFFIXES)))
        ax.set_xticklabels([f"{h}h" for h in TIMEPOINT_HOURS], fontsize=9)
        ax.set_yticks(range(len(task_names)))
        ax.set_yticklabels(task_names, fontsize=8)
        ax.set_title(genotype, fontsize=11, fontweight="bold")
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        for i in range(len(task_names)):
            for j in range(len(STATE_SUFFIXES)):
                feasible = sub.iloc[i, j] > 0.5
                ax.text(j, i, "✓" if feasible else "✗", ha="center", va="center",
                         fontsize=8, color="white", fontweight="bold")
        for j in range(len(STATE_SUFFIXES) + 1):
            ax.axvline(j - 0.5, color="white", linewidth=1.5)
        for i in range(len(task_names) + 1):
            ax.axhline(i - 0.5, color="white", linewidth=1.5)

    legend_elements = [
        Rectangle((0, 0), 1, 1, facecolor=STATUS_GOOD, label="feasible"),
        Rectangle((0, 0), 1, 1, facecolor=STATUS_CRITICAL, label="infeasible"),
    ]
    fig.legend(handles=legend_elements, loc="lower center", ncol=2, frameon=False, fontsize=9,
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Option 1 (expanded, 95-reaction network) task feasibility over the starvation time course", fontsize=12)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(OUTDIR / "feasibility_heatmap.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUTDIR / 'feasibility_heatmap.png'}")

    # --- Summary: which tasks change status at some point ---
    summary_lines = [f"Network: {graph.num_edges} reactions, {len(task_names)} tasks\n"]
    for task_name in task_names:
        row = status_matrix.loc[task_name]
        wt_vals = row[[f"WT_{t}" for t in STATE_SUFFIXES]]
        ob_vals = row[[f"ob_ob_{t}" for t in STATE_SUFFIXES]]
        if wt_vals.nunique() > 1 or ob_vals.nunique() > 1 or not np.array_equal(wt_vals.to_numpy(), ob_vals.to_numpy()):
            wt_traj = "".join("1" if v else "0" for v in wt_vals)
            ob_traj = "".join("1" if v else "0" for v in ob_vals)
            summary_lines.append(f"{task_name}:\n  WT    {wt_traj} (0h->24h)\n  ob_ob {ob_traj} (0h->24h)")
    summary = "\n".join(summary_lines) + "\n"
    (OUTDIR / "build_summary.txt").write_text(summary)
    print("\n" + summary)


if __name__ == "__main__":
    main()
