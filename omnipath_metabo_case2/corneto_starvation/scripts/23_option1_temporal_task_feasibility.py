"""Option 1 task feasibility, per (genotype, timepoint) STATE -- extends
Script 22 (which pooled all 8 timepoints into one static per-genotype
snapshot) to preserve the temporal trajectory. Pooling hid exactly the
kind of signal this starvation time course is about: a gene can be normal
at T0 (fed) and only drop low by T24 (late fasting), and WT vs ob/ob might
diverge in WHEN that happens, not just WHETHER it eventually does.

Mechanism (same two fixes from Script 22, now applied per state instead of
pooled across the genotype's 8 states):
  1. Low-expression threshold is still computed ONCE per genotype (its own
     full 8-state x network-genes distribution, bottom decile -- same
     LOW_EXPRESSION_QUANTILE=0.10 reasoning as Script 22: a fixed reference
     line, not a per-timepoint-local one, so crossing it is comparable
     across time instead of each timepoint grading itself on its own curve).
  2. That FIXED threshold is then applied to each INDIVIDUAL state's own
     value (not that state's median across other states -- there's only one
     state now) to call a gene low AT THAT TIMEPOINT specifically.
  3. All 12 of Script 21's tasks solved (single-task feasibility only, not
     the joint union -- a per-timepoint minimum-network union isn't
     obviously meaningful here, this is about tracking WHEN each task
     individually becomes feasible/infeasible) for each of the 16 states.

Outputs (corneto_starvation/result/option1_temporal_task_feasibility/):
  gene_detection_by_state.csv   per (genotype, timepoint, reaction): low genes, available?
  task_feasibility_by_state.csv per (genotype, timepoint, task): feasible?
  feasibility_heatmap.png       task x timepoint status heatmap, WT and ob_ob panels
  build_summary.txt             counts, which tasks ever flip status, when
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

script21 = importlib.import_module("21_option1_expanded_toy_model")  # noqa: E402
script22 = importlib.import_module("22_option1_genotype_task_feasibility")  # noqa: E402

PKG_ROOT = Path(__file__).resolve().parents[1]
OUTDIR = PKG_ROOT / "result" / "option1_temporal_task_feasibility"

GENOTYPES = ("WT", "ob_ob")
STATE_SUFFIXES = ("0h", "2h", "4h", "6h", "8h", "12h", "16h", "24h")
TIMEPOINT_HOURS = (0, 2, 4, 6, 8, 12, 16, 24)

# dataviz skill reference palette: status colors (reserved, not reused for
# categorical series) -- good/critical, always paired with a label/icon.
STATUS_GOOD = "#0ca30c"
STATUS_CRITICAL = "#d03b3b"


def low_expression_genes_per_state(
    wide_scores: pd.DataFrame, network_genes: list[str],
) -> dict[str, set[str]]:
    """Same fixed-per-genotype threshold as Script 22
    (LOW_EXPRESSION_QUANTILE, bottom decile of that genotype's own 8-state
    x network-genes distribution), but the per-gene CALL now uses each
    state's own individual value, not the cross-state median -- preserves
    which specific timepoints cross the line."""
    present_genes = [g for g in network_genes if g in wide_scores.index]
    low_by_state: dict[str, set[str]] = {}
    for genotype in GENOTYPES:
        state_cols = [f"{genotype}_{t}" for t in STATE_SUFFIXES]
        state_expression = wide_scores.loc[present_genes, state_cols]
        threshold = np.nanquantile(
            state_expression.to_numpy(dtype=float), script22.LOW_EXPRESSION_QUANTILE
        )
        for state in state_cols:
            values = state_expression[state]
            low_by_state[state] = set(values.index[values <= threshold])
    return low_by_state


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    graph = script21.build_expanded_core_metabolic_network()
    task_bounds_all = script21.build_task_bounds()
    task_names = list(task_bounds_all)
    print(f"Network: {graph.num_edges} reactions, {len(task_names)} tasks, "
          f"{len(GENOTYPES) * len(STATE_SUFFIXES)} states")

    network_genes = script22.network_mouse_genes(graph)
    wide_scores = pd.read_csv(script22.GENE_SCORES_PATH, index_col=0)
    low_by_state = low_expression_genes_per_state(wide_scores, network_genes)

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
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)
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
    fig.suptitle("Option 1 task feasibility over the starvation time course", fontsize=12)
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
