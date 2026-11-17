"""Test task feasibility at each (genotype, timepoint) state, to follow how
metabolic capacity changes over the starvation time course.

Close a reaction when its enzyme is low-expressed, then solve each task per
state. A gene can be normal at 0h (fed) and only drop low by 24h (late
fasting), and WT and ob/ob may differ in when that happens, so feasibility is
resolved per timepoint rather than pooled. The low-expression threshold is set
once per genotype (bottom decile of its 8-state distribution, a fixed
reference line), but the call uses each state's own value. The 20 tasks are
solved at each of the 16 states.

The per-genotype overview in build_summary.txt is read off the real states
(tasks feasible at every timepoint, never feasible, and feasible at the 24h
endpoint), not a pooled timepoint median -- in a time course the median
corresponds to no real physiological state.

Outputs (result/temporal_task_feasibility/):
  gene_detection_by_state.csv, task_feasibility_by_state.csv,
  feasibility_status_matrix.csv, feasibility_heatmap.png, build_summary.txt
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
from corneto.methods.metabolism._utils import evaluate_gpr_expression, get_genes_from_gpr  # noqa: E402

network = importlib.import_module("01_metabolic_task_network")  # noqa: E402

PKG_ROOT = Path(__file__).resolve().parents[1]
OUTDIR = PKG_ROOT / "result" / "temporal_task_feasibility"
GENE_SCORES_PATH = PKG_ROOT / "result" / "gene_scores" / "all_states_gene_scores.csv"

GENOTYPES = ("WT", "ob_ob")
STATE_SUFFIXES = ("0h", "2h", "4h", "6h", "8h", "12h", "16h", "24h")
TIMEPOINT_HOURS = (0, 2, 4, 6, 8, 12, 16, 24)
LOW_EXPRESSION_QUANTILE = 0.10

# dataviz skill reference palette: status colors (reserved, not reused for
# categorical series) -- good/critical, always paired with a label/icon.
STATUS_GOOD = "#0ca30c"
STATUS_CRITICAL = "#d03b3b"


def network_mouse_genes(graph) -> list[str]:
    """All unique mouse-converted gene symbols referenced in the network's GPRs
    (human all-caps symbols .capitalize()-d to mouse Title case)."""
    genes: set[str] = set()
    for i in range(graph.num_edges):
        gpr = graph.get_attr_edge(i).get("GPR") or ""
        if gpr.strip():
            genes.update(g.capitalize() for g in get_genes_from_gpr(gpr))
    return sorted(genes)


def reaction_availability(graph, low: set[str]) -> dict[str, dict]:
    """Per-reaction GPR evaluation against a low-gene set. default_value=1
    (available) for any gene not in `low`, including genes never measured.
    Reactions with no GPR are always available."""
    results = {}
    for i in range(graph.num_edges):
        attrs = graph.get_attr_edge(i)
        reaction_id = attrs["id"]
        gpr = attrs.get("GPR") or ""
        if not gpr.strip():
            results[reaction_id] = {"gpr": "", "available": True, "low_genes": ""}
            continue
        genes = set(get_genes_from_gpr(gpr))
        mouse_genes = {g: g.capitalize() for g in genes}
        symbol_values = {mouse_genes[g]: int(mouse_genes[g] not in low) for g in genes}
        mouse_gpr = gpr
        for human, mouse in mouse_genes.items():
            mouse_gpr = mouse_gpr.replace(human, mouse)
        available = bool(evaluate_gpr_expression([mouse_gpr], symbol_values, default_value=1)[0])
        low_here = [mouse_genes[g] for g in genes if mouse_genes[g] in low]
        results[reaction_id] = {
            "gpr": gpr, "available": available, "low_genes": ";".join(sorted(low_here)),
        }
    return results


def low_expression_genes_per_state(
    wide_scores: pd.DataFrame, network_genes: list[str],
) -> dict[str, set[str]]:
    """Fixed per-genotype threshold (bottom decile of that genotype's own
    8-state x network-genes distribution), but the per-gene CALL uses each
    state's own individual value -- preserves which timepoints cross the line."""
    present_genes = [g for g in network_genes if g in wide_scores.index]
    low_by_state: dict[str, set[str]] = {}
    for genotype in GENOTYPES:
        state_cols = [f"{genotype}_{t}" for t in STATE_SUFFIXES]
        state_expression = wide_scores.loc[present_genes, state_cols]
        threshold = np.nanquantile(
            state_expression.to_numpy(dtype=float), LOW_EXPRESSION_QUANTILE
        )
        for state in state_cols:
            values = state_expression[state]
            low_by_state[state] = set(values.index[values <= threshold])
    return low_by_state


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    graph = network.build_full_metabolic_network()
    task_bounds_all = network.build_task_bounds()
    task_names = list(task_bounds_all)
    print(f"Network: {graph.num_edges} reactions, {len(task_names)} tasks, "
          f"{len(GENOTYPES) * len(STATE_SUFFIXES)} states")

    network_genes = network_mouse_genes(graph)
    wide_scores = pd.read_csv(GENE_SCORES_PATH, index_col=0)
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
            availability = reaction_availability(graph, low_by_state[state])
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
    fig.suptitle("Task feasibility over the starvation time course", fontsize=12)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(OUTDIR / "feasibility_heatmap.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUTDIR / 'feasibility_heatmap.png'}")

    # --- Per-genotype overview from real states (no pooling): tasks feasible
    #     at every timepoint, never feasible, and feasible at the 24h endpoint.
    summary_lines = [f"Network: {graph.num_edges} reactions, {len(task_names)} tasks\n"]
    for genotype in GENOTYPES:
        cols = [f"{genotype}_{t}" for t in STATE_SUFFIXES]
        sub = status_matrix[cols]
        throughout = [t for t in task_names if sub.loc[t].all()]
        never = [t for t in task_names if not sub.loc[t].any()]
        at_24h = [t for t in task_names if status_matrix.loc[t, f"{genotype}_24h"]]
        summary_lines.append(
            f"{genotype}: {len(throughout)}/{len(task_names)} feasible throughout, "
            f"{len(never)} never feasible, {len(at_24h)}/{len(task_names)} feasible at 24h"
        )
        summary_lines.append(f"  never feasible: {', '.join(never) or '(none)'}")
    summary_lines.append("")

    # --- Tasks that change status over time or differ between genotypes ---
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
