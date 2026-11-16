"""Script 22, rerun on Script 24's expanded 95-reaction network (full urea
cycle + BCAA catabolism + methionine cycle, 15 tasks instead of 12).
Mechanism is IDENTICAL to Script 22 -- same genotype-pooled low-expression
gating, same LOW_EXPRESSION_QUANTILE=0.10, same default-available-if-
unmeasured convention -- only the network and task set changed. See Script
22's own docstring for the full methodological reasoning; not repeated
here.

Outputs (corneto_starvation/result/option1_genotype_task_feasibility_v2/):
  gene_detection.csv, task_feasibility.csv, selected_reactions_<genotype>.csv,
  build_summary.txt -- same set as Script 22.
"""

from __future__ import annotations

import sys
from pathlib import Path

PKG_SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(PKG_SRC))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import importlib  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from corneto.methods import MultiSampleFBA  # noqa: E402

script22 = importlib.import_module("22_option1_genotype_task_feasibility")  # noqa: E402
script24 = importlib.import_module("24_option1_amino_nitrogen_expansion")  # noqa: E402

PKG_ROOT = Path(__file__).resolve().parents[1]
OUTDIR = PKG_ROOT / "result" / "option1_genotype_task_feasibility_v2"

GENOTYPES = script22.GENOTYPES
LOW_EXPRESSION_QUANTILE = script22.LOW_EXPRESSION_QUANTILE


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "figures").mkdir(exist_ok=True)

    graph = script24.build_full_metabolic_network()
    task_bounds_all = script24.build_task_bounds()
    reaction_ids = list(graph.get_attr_from_edges("id"))
    print(f"Network: {graph.num_edges} reactions, {len(task_bounds_all)} tasks")

    network_genes = script22.network_mouse_genes(graph)
    print(f"Network GPR genes: {len(network_genes)}")

    wide_scores = pd.read_csv(script22.GENE_SCORES_PATH, index_col=0)
    low_by_genotype = script22.low_expression_genes_per_genotype(wide_scores, network_genes)
    for genotype, genes in low_by_genotype.items():
        print(f"{genotype}: {len(genes)}/{len(network_genes)} network genes called low-expression")

    gene_detection_rows = []
    feasibility_rows = []

    for genotype in GENOTYPES:
        availability = script22.reaction_availability(graph, low_by_genotype[genotype])
        for reaction_id, info in availability.items():
            gene_detection_rows.append({"genotype": genotype, "reaction": reaction_id, **info})

        genotype_closure = {
            rid: (0.0, 0.0) for rid, info in availability.items() if not info["available"]
        }
        print(f"\n{genotype}: {len(genotype_closure)}/{graph.num_edges} reactions closed "
              f"(GPR genes measured and low)")

        per_task_status = {}
        for task_name, bounds in task_bounds_all.items():
            merged_bounds = {**bounds, **genotype_closure}
            method = MultiSampleFBA(lambda_reg=0.0, beta_reg=0.0)
            problem = method.build(graph, reaction_bounds=merged_bounds)
            result = problem.solve(solver="GUROBI", verbosity=0)
            status = str(result.status).lower()
            per_task_status[task_name] = status
            blocking = [
                rid for rid, (lo, hi) in bounds.items()
                if lo > 0 and rid in genotype_closure
            ]
            feasibility_rows.append({
                "genotype": genotype, "task": task_name, "status": status,
                "blocked_by_missing_genes": ";".join(blocking),
            })
            print(f"  {task_name}: {status}" + (f"  (blocked: {blocking})" if blocking and status != 'optimal' else ""))

        if all(s == "optimal" for s in per_task_status.values()):
            task_method = MultiSampleFBA(lambda_reg=1.0, beta_reg=0.0)
            merged_task_bounds = {
                task: {**bounds, **genotype_closure} for task, bounds in task_bounds_all.items()
            }
            task_problem = task_method.build_many(graph, reaction_bounds=merged_task_bounds)
            task_problem.add_objective(
                task_problem.expr.edge_has_flux_OR.sum(), name="explicit_minimum_network_size"
            )
            first_solve = task_problem.solve(solver="GUROBI", verbosity=0)
            if str(first_solve.status).lower() == "optimal":
                selected = np.asarray(task_problem.expr.edge_has_flux_OR.value).ravel() > 0.5
                selected_df = pd.DataFrame({
                    "reaction": np.asarray(reaction_ids)[selected],
                }).set_index("reaction")
                selected_df.to_csv(OUTDIR / f"selected_reactions_{genotype}.csv")
                print(f"  Joint minimum network for {genotype}: {int(selected.sum())}/{graph.num_edges} reactions")
        else:
            infeasible = [t for t, s in per_task_status.items() if s != "optimal"]
            print(f"  {genotype}: NOT all tasks jointly feasible (infeasible: {infeasible}); "
                  f"skipping joint union solve")

    gene_detection_df = pd.DataFrame(gene_detection_rows)
    gene_detection_df.to_csv(OUTDIR / "gene_detection.csv", index=False)
    feasibility_df = pd.DataFrame(feasibility_rows)
    feasibility_df.to_csv(OUTDIR / "task_feasibility.csv", index=False)

    summary_lines = [f"Network: {graph.num_edges} reactions, {len(task_bounds_all)} tasks\n"]
    for genotype in GENOTYPES:
        sub = feasibility_df[feasibility_df["genotype"] == genotype]
        n_feasible = (sub["status"] == "optimal").sum()
        summary_lines.append(f"{genotype}: {n_feasible}/{len(sub)} tasks feasible")
        infeasible_tasks = sub[sub["status"] != "optimal"]
        for _, row in infeasible_tasks.iterrows():
            summary_lines.append(f"  INFEASIBLE: {row['task']} (blocked by: {row['blocked_by_missing_genes']})")
    summary = "\n".join(summary_lines) + "\n"
    (OUTDIR / "build_summary.txt").write_text(summary)
    print("\n" + summary)


if __name__ == "__main__":
    main()
