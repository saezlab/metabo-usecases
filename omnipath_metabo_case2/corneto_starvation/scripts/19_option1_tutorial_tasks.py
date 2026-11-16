"""Faithful reproduction of the colleague's notebook Option 1
(notebooks/model-example.ipynb, cells 13-18: "Example option 1: metabolic
tasks on a small graph"), using the tutorial's own 58-reaction network and
the tutorial's own 6 tasks -- no custom curation yet. This is step one of
two: confirm the literal example reproduces correctly before deciding
whether the generic tasks need replacing with starvation-specific ones for
our own dataset.

Network: Script 06's `build_core_metabolic_network()` (corneto_mcmt.
flux_benchmark) -- same ~58-reaction hand-curated hepatic core model used
throughout the notebook. Note this function also returns its own two
built-in MetabolicTasks (fed_glycolytic_lipogenic / fasted_gluconeogenic_
ketogenic) -- those are NOT used here; only the graph is taken, and the
notebook's own 6 ad hoc tasks (defined directly in its Option 1 cell) are
used instead, exactly as written there.

Method (cells 15-16): MultiSampleFBA(lambda_reg=1.0, beta_reg=0.0),
build_many() across the 6 tasks as 6 samples (a genuinely joint multi-
sample union -- unlike Option 2, which is single-scenario). lambda_reg's
structural-sparsity term is supposed to penalize the union of reactions
ever selected across all 6 tasks, but per the notebook's own comment, the
pinned CORNETO version (1.0.0rc8, same version this project uses) needs
that count added EXPLICITLY as a second objective
(`edge_has_flux_OR.sum()`) to actually have the solver apply it -- the
automatic lambda_reg term alone doesn't do so reliably in this version.

Two-stage solve, exactly as the notebook does it:
  1. Minimize the reaction union subject to all 6 tasks being feasible.
  2. Fix the resulting minimum network size, then minimize total absolute
     flux (a second continuous variable, |flow| via lb/ub constraints) for
     a parsimonious, readable flux solution at that fixed network size --
     this doesn't change which reactions are selected, only picks among
     possibly-many flux solutions that use exactly that reaction set.

This is a genuinely different computational profile from Option 2: only
58 reactions (vs Mouse-GEM's 12,987), so the same joint-multi-sample
mechanic that failed catastrophically at full-genome scale (first Option 2
attempt, 14 scenarios x 12,987 reactions) should be trivial here.

Outputs (corneto_starvation/result/option1_tutorial_tasks/):
  task_flux.csv            flux value per reaction per task (58 x 6)
  selected_reactions.csv   the reactions in the minimum joint-feasible union
  task_flux_gallery.html   colleague's plot_task_flux_gallery rendering
  build_summary.txt        counts, solver status for both solve stages
"""

from __future__ import annotations

import sys
from pathlib import Path

PKG_SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(PKG_SRC))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from corneto.backend import VarType  # noqa: E402
from corneto.methods import MultiSampleFBA  # noqa: E402

from corneto_mcmt.flux_benchmark import build_core_metabolic_network  # noqa: E402
from corneto_mcmt.model_example import plot_task_flux_gallery  # noqa: E402

PKG_ROOT = Path(__file__).resolve().parents[1]
OUTDIR = PKG_ROOT / "result" / "option1_tutorial_tasks"

# Exactly the notebook's Option 1 cell -- same reaction names, same bounds.
INPUTS = ("EX_GLC_IN", "EX_FA_IN", "EX_LAC_IN", "EX_ALA_IN")
OUTPUTS = (
    "EX_GLC_OUT", "EX_LAC_OUT", "EX_LIPID_OUT", "EX_KET_OUT",
    "EX_GLYCOGEN_OUT", "EX_HEXOSAMINE_OUT",
)
# Keep this first example focused on production rather than substrate cycles.
BLOCKED_CYCLES = ("FUTILE_ATPASE", "GLYG_DEG", "G6PC", "FBP1")


def build_task_bounds() -> dict[str, dict[str, tuple[float, float]]]:
    base_bounds = {reaction: (0.0, 0.0) for reaction in (*INPUTS, *OUTPUTS, *BLOCKED_CYCLES)}
    base_bounds["ATP_MAINT"] = (1.0, 1.0)
    return {
        "Glucose to lactate": {**base_bounds, "EX_GLC_IN": (0.0, 10.0), "EX_LAC_OUT": (1.0, 10.0)},
        "Glucose to glycogen": {**base_bounds, "EX_GLC_IN": (0.0, 10.0), "EX_GLYCOGEN_OUT": (1.0, 10.0)},
        "Glucose to lipid": {**base_bounds, "EX_GLC_IN": (0.0, 10.0), "EX_LIPID_OUT": (1.0, 10.0)},
        "Fatty acid to ketone": {**base_bounds, "EX_FA_IN": (0.0, 10.0), "EX_KET_OUT": (1.0, 10.0)},
        "Glucose to hexosamine": {**base_bounds, "EX_GLC_IN": (0.0, 10.0), "EX_HEXOSAMINE_OUT": (1.0, 10.0)},
        "Lactate to lipid": {**base_bounds, "EX_LAC_IN": (0.0, 10.0), "EX_LIPID_OUT": (1.0, 10.0)},
    }


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    graph, _default_tasks = build_core_metabolic_network()
    print(f"Tutorial network: {graph.num_edges} reactions")

    task_bounds = build_task_bounds()
    print("Tasks:", ", ".join(task_bounds))

    task_method = MultiSampleFBA(lambda_reg=1.0, beta_reg=0.0)
    task_problem = task_method.build_many(graph, reaction_bounds=task_bounds)
    task_problem.add_objective(
        task_problem.expr.edge_has_flux_OR.sum(), name="explicit_minimum_network_size"
    )
    first_solve = task_problem.solve(solver="GUROBI", verbosity=1)
    status_1 = str(first_solve.status).lower()
    print(f"Stage 1 (minimize union) solve status: {status_1}")
    if status_1 != "optimal":
        raise RuntimeError(f"The six tasks could not be satisfied: {first_solve.status}")

    minimum_size = int(np.count_nonzero(np.asarray(task_problem.expr.edge_has_flux_OR.value) > 0.5))
    flow = task_problem.expr.flow
    absolute_flow = task_method.backend.Variable(
        name="task_absolute_flux", shape=flow.shape, lb=0.0, vartype=VarType.CONTINUOUS
    )
    task_problem.add_constraints([
        task_problem.expr.edge_has_flux_OR.sum() == minimum_size,
        absolute_flow >= flow,
        absolute_flow >= -flow,
    ])
    task_problem.add_objective(absolute_flow.sum(), name="parsimonious_task_flux")
    second_solve = task_problem.solve(solver="GUROBI", verbosity=1)
    status_2 = str(second_solve.status).lower()
    print(f"Stage 2 (minimize flux at fixed size) solve status: {status_2}")
    if status_2 != "optimal":
        raise RuntimeError(f"Flux refinement failed: {second_solve.status}")

    reaction_ids = list(graph.get_attr_from_edges("id"))
    task_flux = pd.DataFrame(
        np.asarray(task_problem.expr.flow.value), index=reaction_ids, columns=task_bounds
    )
    task_flux.to_csv(OUTDIR / "task_flux.csv")

    selected_reactions = np.asarray(task_problem.expr.edge_has_flux_OR.value).ravel() > 0.5
    assert selected_reactions.sum() == minimum_size
    assert np.array_equal(selected_reactions, (np.abs(task_flux.to_numpy()) > 1e-6).any(axis=1))

    print(f"Smallest joint network: {minimum_size} of {graph.num_edges} reactions")
    exchange_summary = task_flux.loc[
        ["EX_GLC_IN", "EX_FA_IN", "EX_LAC_IN", "EX_LAC_OUT", "EX_GLYCOGEN_OUT",
         "EX_LIPID_OUT", "EX_KET_OUT", "EX_HEXOSAMINE_OUT", "ATP_MAINT"]
    ].T.round(3)
    print(exchange_summary)

    selected_df = pd.DataFrame({
        "reaction": np.asarray(reaction_ids)[selected_reactions],
    }).set_index("reaction")
    selected_df.to_csv(OUTDIR / "selected_reactions.csv")

    html = plot_task_flux_gallery(
        graph, task_flux.to_numpy(), list(task_bounds), selected_reactions,
    )
    (OUTDIR / "task_flux_gallery.html").write_text(f"<html><body>{html}</body></html>")

    summary = (
        f"Tutorial network: {graph.num_edges} reactions total\n"
        f"Tasks: {', '.join(task_bounds)}\n"
        f"Stage 1 (minimize union) solve status: {status_1}\n"
        f"Stage 2 (minimize flux at fixed size) solve status: {status_2}\n"
        f"Smallest joint network: {minimum_size}/{graph.num_edges} reactions\n"
    )
    (OUTDIR / "build_summary.txt").write_text(summary)
    print("\n" + summary)


if __name__ == "__main__":
    main()
