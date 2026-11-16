"""Option 1 task feasibility, per genotype -- still the same small
hand-curated network (Script 21's 69-reaction expansion, NOT the full
Mouse-GEM; confirmed Option 1 is defined as the small-network method, see
Script 21's docstring), but now gated by each genotype's OWN measured
proteomics.

A first attempt gated reactions on simple gene DETECTION (any finite
measured value for that gene, any state) and found every single task
infeasible for BOTH genotypes identically, with exactly the same 3,452
genes "detected" in each. That's not a real biological signal -- this
proteomics panel covers 3,452 of the mouse genome's >20,000 genes (~17%,
typical of shotgun MS depth), and the reactions it zeroed out include
basic housekeeping enzymes (Gpi, Eno1, Pc) that are certainly expressed in
real liver but just aren't in THIS experiment's coverage. "Not detected"
was being conflated with "not expressed" -- wrong default. Fixed here:
genes absent from the panel entirely default to AVAILABLE (unknown, not
penalized); only genes that ARE measured and fall in the LOW-expression
range (same threshold convention as Scripts 03/16) count as evidence the
enzyme is genuinely down, and close the reaction.

Mechanism:
  1. For each genotype, pool that genotype's own 8 states' log2 abundance
     per gene (median, skipna) and apply a strict bottom-decile threshold
     (LOW_EXPRESSION_QUANTILE = 0.10 -- see that constant's comment for why
     Scripts 03/16's own (0.50, 0.76) convention is too aggressive for a
     hard-knockout use case) to call each MEASURED gene "low" or not. Genes absent
     from the 3,452-gene panel are NOT called low -- default assumed
     available.
  2. Each reaction's GPR string in Script 21's network uses HUMAN gene
     symbol convention (all-caps, e.g. "PCK1" -- the tutorial's own
     choice). Mouse orthologs are the same symbol in Title case for the
     vast majority of genes (`"PCK1".capitalize() == "Pck1"`), so GPR
     symbols are converted with `.capitalize()` before matching against
     the (already mouse-symbol) proteomics index.
  3. `corneto.methods.metabolism._utils.evaluate_gpr_expression` evaluates
     each reaction's AND/OR GPR logic against that genotype's low-gene set
     (0 i.e. "absent/low" only for genes measured-and-low; default 1 i.e.
     "available" for every other gene, including ones never measured at
     all). Reactions with no GPR (pure exchange/maintenance reactions with
     no named enzyme in this toy model) are always available.
  4. Reactions whose GPR evaluates to unsupported are forced to (0, 0) for
     that genotype, THEN merged with each task's own bounds with the
     genotype closure taking priority: `{**task_bounds[task],
     **genotype_closure}`. If a task specifically requires a reaction that
     genotype closure has just zeroed out, that task becomes infeasible
     for that genotype -- the real, informative outcome this script is
     built to surface (not a bug to work around).
  5. All 12 tasks solved per genotype, individually (not jointly across
     genotypes -- there is no biological reason to force WT and ob/ob onto
     the same reaction union, unlike Script 21's joint union across
     scenarios of the SAME organism).

Outputs (corneto_starvation/result/option1_genotype_task_feasibility/):
  gene_detection.csv       per reaction, per genotype: GPR string, genes
                            used, detected?, reaction available?
  task_feasibility.csv     per task, per genotype: feasible? (and which
                            specific genotype-closed reactions blocked it,
                            if infeasible)
  selected_reactions_<genotype>.csv   per-genotype minimum joint network
                                       (only for genotypes where all tasks
                                       are jointly feasible)
  build_summary.txt        counts, solver status per genotype
"""

from __future__ import annotations

import sys
from pathlib import Path

PKG_SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(PKG_SRC))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from corneto.backend import VarType  # noqa: E402
from corneto.methods import MultiSampleFBA  # noqa: E402
from corneto.methods.metabolism._utils import evaluate_gpr_expression, get_genes_from_gpr  # noqa: E402

import importlib  # noqa: E402
script21 = importlib.import_module("21_option1_expanded_toy_model")  # noqa: E402

PKG_ROOT = Path(__file__).resolve().parents[1]
GENE_SCORES_PATH = PKG_ROOT / "result" / "gene_scores" / "all_states_gene_scores.csv"
OUTDIR = PKG_ROOT / "result" / "option1_genotype_task_feasibility"

GENOTYPES = ("WT", "ob_ob")
STATE_SUFFIXES = ("0h", "2h", "4h", "6h", "8h", "12h", "16h", "24h")


# Script 03/16's own (0.50, 0.76) convention calls roughly HALF of all
# genes "low" by construction (0.50 is the median) -- fine for soft iMAT
# scoring spread across a large, redundant network where most reactions
# tolerate a "low" call without becoming globally infeasible. It is far
# too aggressive for a HARD knockout in a small, tightly-coupled ~69-
# reaction network: any single required reaction has a ~50% chance of
# being "low" by construction, and most tasks need several enzymes in an
# unbranched sequence, so at least one link getting cut became almost
# certain -- a first attempt at this script found EVERY task infeasible
# for BOTH genotypes using 0.50, which is a quantile-construction artifact,
# not a real genotype signal. Tightened to the bottom decile (0.10): only
# genuinely extreme low outliers count as "the enzyme probably isn't
# functioning", not "below the middle of the pack".
LOW_EXPRESSION_QUANTILE = 0.10


def low_expression_genes_per_genotype(
    wide_scores: pd.DataFrame, network_genes: list[str],
) -> dict[str, set[str]]:
    """Genes MEASURED (present in the proteomics panel) and falling at or
    below a low-quantile threshold -- genuine evidence of low expression,
    not mere absence from the panel's coverage. Threshold is computed the
    same way Scripts 03/16 do (global quantile over the flattened genes x
    states matrix) but with a stricter cutoff (see LOW_EXPRESSION_QUANTILE)
    and scoped to just this network's own GPR genes (not the full
    3,452-gene panel) -- "low relative to what's typical among genes this
    pathway actually uses", matching Script 03/16's own scoping to
    network-relevant genes rather than the whole proteome."""
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
    """All unique mouse-converted gene symbols referenced in the network's GPRs."""
    genes: set[str] = set()
    for i in range(graph.num_edges):
        gpr = graph.get_attr_edge(i).get("GPR") or ""
        if gpr.strip():
            genes.update(g.capitalize() for g in get_genes_from_gpr(gpr))
    return sorted(genes)


def reaction_availability(graph, low: set[str]) -> dict[str, dict]:
    """Per-reaction: GPR string, human->mouse-converted genes, availability.
    default_value=1 (available) for any gene not in `low` -- includes genes
    never measured at all, which should NOT be treated as evidence of
    absence (see module docstring)."""
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
        # evaluate_gpr_expression expects values keyed by the SAME symbols
        # appearing in the GPR string; substitute the mouse-converted GPR text.
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

    graph = script21.build_expanded_core_metabolic_network()
    task_bounds_all = script21.build_task_bounds()
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
