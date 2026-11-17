"""Test which metabolic tasks each genotype can carry out, given its own
measured proteomics.

We take the task network from 01_metabolic_task_network and close a reaction
when its enzyme is low-expressed: a gene counts as low only if it is measured
AND falls at or below the bottom decile of that genotype's expression over
the network's GPR genes and 8 timepoints. Genes absent from the proteomics
panel stay available (not detected is not the same as not expressed). A
reaction whose GPR then evaluates to unsupported is fixed to zero flux for
that genotype; a task that needs it becomes infeasible. The 20 tasks are
solved per genotype.

Outputs (result/genotype_task_feasibility/):
  gene_detection.csv, task_feasibility.csv, selected_reactions_<genotype>.csv,
  build_summary.txt
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
from corneto.methods.metabolism._utils import evaluate_gpr_expression, get_genes_from_gpr  # noqa: E402

network = importlib.import_module("01_metabolic_task_network")  # noqa: E402

PKG_ROOT = Path(__file__).resolve().parents[1]
OUTDIR = PKG_ROOT / "result" / "genotype_task_feasibility"
GENE_SCORES_PATH = PKG_ROOT / "result" / "gene_scores" / "all_states_gene_scores.csv"

GENOTYPES = ("WT", "ob_ob")
STATE_SUFFIXES = ("0h", "2h", "4h", "6h", "8h", "12h", "16h", "24h")
LOW_EXPRESSION_QUANTILE = 0.10


def low_expression_genes_per_genotype(
    wide_scores: pd.DataFrame, network_genes: list[str],
) -> dict[str, set[str]]:
    """Genes MEASURED (present in the proteomics panel) and at or below the
    bottom-decile threshold of that genotype's own 8-state x network-genes
    distribution -- genuine low-expression evidence, not panel absence."""
    present_genes = [g for g in network_genes if g in wide_scores.index]
    low: dict[str, set[str]] = {}
    for genotype in GENOTYPES:
        state_cols = [f"{genotype}_{t}" for t in STATE_SUFFIXES]
        state_expression = wide_scores.loc[present_genes, state_cols]
        threshold = np.nanquantile(state_expression.to_numpy(dtype=float), LOW_EXPRESSION_QUANTILE)
        gene_median = state_expression.median(axis=1, skipna=True)
        low[genotype] = set(gene_median.index[gene_median <= threshold])
    return low


def network_mouse_genes(graph) -> list[str]:
    """All unique mouse-converted gene symbols referenced in the network's GPRs.
    GPRs use human all-caps symbols (tutorial convention); mouse orthologs are
    the same symbol in Title case for the vast majority, so `.capitalize()`."""
    genes: set[str] = set()
    for i in range(graph.num_edges):
        gpr = graph.get_attr_edge(i).get("GPR") or ""
        if gpr.strip():
            genes.update(g.capitalize() for g in get_genes_from_gpr(gpr))
    return sorted(genes)


def reaction_availability(graph, low: set[str]) -> dict[str, dict]:
    """Per-reaction GPR evaluation against a low-gene set. default_value=1
    (available) for any gene not in `low`, including genes never measured --
    absence from the panel must NOT be read as evidence of absence. Reactions
    with no GPR (pure exchange/maintenance) are always available."""
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


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "figures").mkdir(exist_ok=True)

    graph = network.build_full_metabolic_network()
    task_bounds_all = network.build_task_bounds()
    reaction_ids = list(graph.get_attr_from_edges("id"))
    print(f"Network: {graph.num_edges} reactions, {len(task_bounds_all)} tasks")

    network_genes = network_mouse_genes(graph)
    print(f"Network GPR genes: {len(network_genes)}")

    wide_scores = pd.read_csv(GENE_SCORES_PATH, index_col=0)
    low_by_genotype = low_expression_genes_per_genotype(wide_scores, network_genes)
    for genotype, genes in low_by_genotype.items():
        print(f"{genotype}: {len(genes)}/{len(network_genes)} network genes called low-expression")

    gene_detection_rows = []
    feasibility_rows = []

    for genotype in GENOTYPES:
        availability = reaction_availability(graph, low_by_genotype[genotype])
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
