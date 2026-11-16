"""Option 1 (metabolic-task network curation) with our OWN starvation-
relevant tasks, replacing the tutorial's generic 6 tasks (Script 19:
glucose-to-lactate/glycogen/lipid/hexosamine, lactate-to-lipid -- mostly
FED-state carbohydrate/lipid handling, not what this dataset is about).

Same network as Script 19 (the tutorial's 58-reaction hand-curated hepatic
core model, `build_core_metabolic_network()`) and the exact same two-stage
MultiSampleFBA(lambda_reg=1.0) mechanism -- only the task set changes.

Tasks, each representing one established liver function during fasting:
  Gluconeogenesis from lactate   -- Cori cycle (muscle lactate -> liver glucose)
  Ketogenesis from fatty acid    -- liver's main fasting fuel export (kept
                                     from Script 19 -- already fasting-relevant)
  Fatty acid oxidation for energy -- fat alone meeting an elevated ATP
                                     demand (ATP_MAINT forced to (5,10)
                                     instead of the baseline (1,1)), with
                                     no other carbon source allowed
  Combined fasted state          -- FA + lactate + alanine simultaneously
                                     available, must export BOTH glucose and
                                     ketone at once (the network's own
                                     built-in `fasted_gluconeogenic_ketogenic`
                                     MetabolicTask from build_core_metabolic_
                                     network(), merged with our blocked-cycle
                                     convention for consistency with the
                                     other three tasks)

Tried but dropped -- NOT representable in this 58-reaction toy network
(confirmed by testing each task in isolation, not a bounds-design error):
  Gluconeogenesis from alanine -- ALT_IN (Ala_c -> Pyr_c + NADH_c) doesn't
    consume NAD_c the way LDH_R (Lac_c + NAD_c -> Pyr_c + NADH_c) does on
    the lactate route, so the NAD_c that GAPDH_PGK_R produces downstream
    (shared by both routes) has no consumer when alanine is the sole input
    -- NAD_c isn't an exchanged metabolite, so it must balance to exactly
    zero internally. A stoichiometric gap in this simplified model, not a
    biological impossibility.
  Glycogenolysis -- no pre-existing glycogen POOL to deplete, only
    glucose<->glycogen cycling reactions with no exchange for a starting store.
  Ureagenesis -- no urea cycle / nitrogen disposal reactions modeled at all.
  Both would need the full Mouse-GEM, not this toy model.

Outputs (corneto_starvation/result/option1_starvation_tasks/): same set as
Script 19.
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
OUTDIR = PKG_ROOT / "result" / "option1_starvation_tasks"

INPUTS = ("EX_GLC_IN", "EX_FA_IN", "EX_LAC_IN", "EX_ALA_IN")
OUTPUTS = (
    "EX_GLC_OUT", "EX_LAC_OUT", "EX_LIPID_OUT", "EX_KET_OUT",
    "EX_GLYCOGEN_OUT", "EX_HEXOSAMINE_OUT",
)
# G6PC/FBP1 (the committed gluconeogenesis dephosphorylation steps) are
# deliberately NOT blocked here, unlike Script 19's tutorial tasks -- those
# never needed to produce free glucose (glucose was always an input), so
# blocking them was harmless there. Our gluconeogenesis-from-lactate/alanine
# tasks need exactly these two steps to reach free glucose; blocking them
# (a first attempt copied Script 19's list unchanged) made every
# gluconeogenesis task infeasible.
BLOCKED_CYCLES = ("FUTILE_ATPASE", "GLYG_DEG")


def build_task_bounds() -> dict[str, dict[str, tuple[float, float]]]:
    base_bounds = {reaction: (0.0, 0.0) for reaction in (*INPUTS, *OUTPUTS, *BLOCKED_CYCLES)}
    base_bounds["ATP_MAINT"] = (1.0, 1.0)
    return {
        "Gluconeogenesis from lactate": {
            **base_bounds, "EX_LAC_IN": (0.0, 10.0), "EX_GLC_OUT": (1.0, 10.0),
        },
        # "Gluconeogenesis from alanine" dropped: infeasible in this toy
        # network specifically (not a bounds-design error). ALT_IN
        # (Ala_c -> Pyr_c + NADH_c) doesn't consume NAD_c the way LDH_R
        # (Lac_c + NAD_c -> Pyr_c + NADH_c) does, so the NAD_c that
        # GAPDH_PGK_R produces downstream (part of the shared gluconeogenic
        # pathway) has no consumer on the alanine route -- NAD_c isn't an
        # exchanged metabolite, so it must balance to zero internally. A
        # stoichiometric gap in this simplified 58-reaction model, confirmed
        # by testing each task in isolation.
        "Ketogenesis from fatty acid": {
            **base_bounds, "EX_FA_IN": (0.0, 10.0), "EX_KET_OUT": (1.0, 10.0),
        },
        "Fatty acid oxidation for energy": {
            **base_bounds, "EX_FA_IN": (0.0, 10.0), "ATP_MAINT": (5.0, 10.0),
        },
        "Combined fasted state": {
            **base_bounds,
            "EX_FA_IN": (0.0, 10.0), "EX_LAC_IN": (0.0, 10.0), "EX_ALA_IN": (0.0, 10.0),
            "EX_GLC_OUT": (1.0, 1.0), "EX_KET_OUT": (1.0, 1.0),
            "ATP_MAINT": (2.0, 2.0),
        },
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
        raise RuntimeError(f"The tasks could not be satisfied: {first_solve.status}")

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
        ["EX_GLC_IN", "EX_FA_IN", "EX_LAC_IN", "EX_ALA_IN", "EX_GLC_OUT",
         "EX_KET_OUT", "ATP_MAINT"]
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
