"""Extract a biologically grounded "core network" from the FULL Mouse-GEM,
following the colleague's corrected methodology (he explicitly rejected both
our earlier hand-picked 8-pathway subgraph with ad hoc synthetic sinks, and
his own tutorial's 58-reaction toy network, as unsuitable for real analysis):

  proteomics -> low/high threshold -> map to reactions via GPR rules ->
  reaction scores -> zero-evidence reactions penalized for removal ->
  CORNETO as gap-filling (keep evidence-supported reactions + whatever extra
  reactions are structurally required to support them + a constraint
  guaranteeing a functional minimum like biomass production) -> core network.

Implementation: this is standard iMAT (confirmed with the colleague):
+1/-1/0 reaction scores from GPR evaluation, keep what's scored +1, drop
what's scored -1, and 0 (unknown) is neither rewarded nor penalized when
lambda_reg=0 -- the colleague's recommended starting point, since a
positive lambda_reg additionally penalizes unknown reactions for being
kept, giving an artificially smaller network than plain iMAT would. A
biomass reaction lower bound (reaction_bounds on the model's own
fba_objective reaction, MAR00021, mirroring the feasibility check from the
upstream corneto-metabo-mct tutorial notebook) is still applied. Mass
balance (Sv=0) plus the biomass floor forces the solver to keep any
unscored "connector" reaction that's structurally necessary to reach
biomass -- that is the gap-filling: it falls out of the MILP's own
feasibility requirement, not a separate algorithm.

Evidence here is POOLED across all 16 (genotype, timepoint) states (median
log2 abundance per gene across all 80 samples, globally pooled quantile
thresholds). Quantiles are 0.25/0.75 (symmetric quartiles, the colleague's
suggestion) rather than his own simulator's DEFAULT_PROTEOMICS_PARAMETERS
(0.50/0.76), which is asymmetric and was flagged as needing justification
we didn't have. This is a single "is this gene relevant anywhere in this
dataset" pass to build ONE starting network, not a per-timepoint signal.
Per-timepoint dynamics are Script 03's job, run on the network this script
produces.

IMPORTANT: this script does NOT persist the core network as an extracted
CORNETO subgraph object. Verified directly: CORNETO 1.0.0rc8's
Graph.edge_subgraph()/extract_subgraph() silently flattens any non-unit
stoichiometric coefficient to +-1 when copying edges into the new Graph
(the `stoichiometry` metadata dict attribute is preserved correctly, but the
internal per-vertex edge attributes that get_vertex_incidence_matrix_as_lists
actually reads for mass-balance are not rebuilt from it). A Mouse-GEM
reaction with a coefficient of 2 for a metabolite becomes coefficient 1 in
the extracted subgraph -- confirmed with MAR04283 (H+ coefficient 2 in the
full model, 1 in its own one-edge "subgraph"). This corrupted a first version
of this script: a reaction set validated as biomass-feasible via bounds on
the full graph came back infeasible when re-tested as its own extracted
subgraph object, with no warning. Downstream scripts must keep using the
full Mouse-GEM graph and restrict to the core network via reaction_bounds
(force every non-core reaction to (0, 0)), never via edge_subgraph/
extract_subgraph for anything that will be solved again.

Outputs (corneto_starvation/result/core_network_extraction/):
  core_reactions.csv     reaction IDs kept in the core network, with an
                          `evidence_supported` flag (kept because scored,
                          vs. gap-filled with no direct evidence)
  gene_scores.csv         pooled per-gene score (-1/0/+1) and raw abundance
  reaction_scores.csv     GPR-derived reaction scores on the full Mouse-GEM
  build_summary.txt       counts, biomass flux achieved, solver status
"""

from __future__ import annotations

import sys
from pathlib import Path

PKG_SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(PKG_SRC))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import corneto as cn  # noqa: E402
from corneto.methods.fba import MultiSampleFBA  # noqa: E402
from corneto.methods.imat import MultiSampleIMAT  # noqa: E402

from corneto_mcmt import load_mouse_gem  # noqa: E402

GENE_SCORES_PATH = Path(__file__).resolve().parents[1] / "result" / "gene_scores" / "all_states_gene_scores.csv"
RESULT_ROOT = Path(__file__).resolve().parents[1] / "result" / "core_network_extraction"
# Working directory during the run; renamed to RESULT_ROOT/lr{LAMBDA_REG}_{n_reactions}reactions
# once the final reaction count is known, so every lambda_reg sweep keeps its
# own labeled, reproducible output instead of overwriting the previous run.
OUTDIR = RESULT_ROOT / "_in_progress"

LOW_QUANTILE, HIGH_QUANTILE = 0.25, 0.75  #  symmetric
LAMBDA_REG = 0.5  # 0.01/0.1/0.15 all gave ~4460-4585 reactions (35% of the
# genome) -- barely moves across a 15x change in lambda_reg. 0.2 gave 374
# reactions -- a cliff somewhere between 0.15 and 0.2. Binary-searching:
# trying the midpoint of the untested 0.15-0.2 range.
EPS = 1e-2
BIOMASS_FLOOR = 0.01
MAX_SECONDS = 1200
_FEASIBILITY_TOL = 1e-6


def build_pooled_gene_scores(wide_scores: pd.DataFrame) -> tuple[pd.Series, float, float]:
    """Median log2 abundance per gene across all 16 states, then threshold
    with globally pooled quantiles -- mirrors flux_benchmark.py's own
    thresholding logic (np.nanquantile over the full gene x state matrix)."""
    pooled_abundance = wide_scores.median(axis=1, skipna=True)
    low_thresh, high_thresh = np.nanquantile(wide_scores.to_numpy(dtype=float), [LOW_QUANTILE, HIGH_QUANTILE])
    scores = pd.Series(0.0, index=pooled_abundance.index)
    scores[pooled_abundance >= high_thresh] = 1.0
    scores[pooled_abundance <= low_thresh] = -1.0
    return scores, float(low_thresh), float(high_thresh)


def score_reactions(graph, gene_scores: pd.Series) -> dict[str, float]:
    cdict = {
        "pooled": {
            gene: {"value": float(score), "mapping": "none", "role": "expression"}
            for gene, score in gene_scores.items()
        }
    }
    data = cn.data.Data.from_cdict(cdict)
    scorer = MultiSampleIMAT(high_expression_threshold=1.0, low_expression_threshold=-1.0, gpr_field="GPR")
    _, processed = scorer.preprocess(graph, data)
    scores: dict[str, float] = {}
    for feature in processed.samples["pooled"].features:
        if feature.mapping == "edge" and feature.data.get("role") == "expression" and feature.value is not None:
            scores[feature.id] = float(feature.value)
    return scores


def repair_and_validate(graph, reaction_ids, selected, biomass_reaction, biomass_bounds):
    """Fix the MILP's binary selection, then resolve as a plain continuous LP
    (no integer variables) to get an exact, cleanly mass-balanced flux for
    exactly that reaction set. A timed-out MILP incumbent can carry a tiny
    constraint violation (seen before: 2e-6, just over tolerance); this
    repair makes the accepted network's feasibility unambiguous rather than
    trusting the MILP's own possibly-imprecise continuous values.
    """
    forced_bounds = dict(biomass_bounds)
    for reaction_id, keep in zip(reaction_ids, selected):
        if not keep:
            forced_bounds[reaction_id] = (0.0, 0.0)
    method = MultiSampleFBA()
    problem = method.build_many(graph, objectives={"repair": {}}, reaction_bounds={"repair": forced_bounds})
    result = problem.solve(solver="GUROBI", verbosity=0)
    status = str(result.status).lower()
    if status != "optimal":
        raise RuntimeError(f"Repair LP for the selected reaction set is infeasible: status={status}")
    flux = np.asarray(problem.expr.flow.value, dtype=float).reshape(-1)
    biomass_flux = float(flux[reaction_ids.index(biomass_reaction)])
    if biomass_flux < BIOMASS_FLOOR - _FEASIBILITY_TOL:
        raise RuntimeError(f"Repaired biomass flux {biomass_flux:.6g} is below the required floor {BIOMASS_FLOOR}.")
    return flux, biomass_flux


def extract_core_network(graph, reaction_scores: dict[str, float]) -> tuple:
    graph_attrs = graph.get_graph_attributes()
    biomass_reaction = next(iter(graph_attrs["fba_objective"]))
    reaction_ids = list(graph.get_attr_from_edges("id"))
    default_ub = list(graph.get_attr_from_edges("default_ub"))
    biomass_upper = float(default_ub[reaction_ids.index(biomass_reaction)])
    biomass_bounds = {biomass_reaction: (BIOMASS_FLOOR, biomass_upper)}

    method = MultiSampleIMAT(
        lambda_reg=LAMBDA_REG, beta_reg=0.0, eps=EPS, scale=False, use_bigm_constraints=True,
    )
    problem = method.build_many(
        graph, reaction_scores={"pooled": reaction_scores}, reaction_bounds={"pooled": biomass_bounds}
    )
    result = problem.solve(
        solver="GUROBI", verbosity=1, max_seconds=MAX_SECONDS,
        MIPFocus=1, IntFeasTol=1e-9,
    )
    status = str(result.status).lower()
    print(f"MILP solve status: {status}")

    # lambda_reg>0 (even a small value) makes FlowMethod's structured-sparsity
    # code path create the edge_has_flux binary selector -- the clean
    # structural "is this reaction in the network" indicator, unlike the
    # flux-magnitude proxy needed when lambda_reg=0 (which gave an
    # unreliable selection: two prior attempts with that proxy failed the
    # repair/validation step).
    selector = np.asarray(problem.expr.edge_has_flux.value, dtype=float).reshape(-1)
    selected = np.rint(selector).astype(bool)

    flux, biomass_flux = repair_and_validate(graph, reaction_ids, selected, biomass_reaction, biomass_bounds)
    print(f"Repair LP: clean, validated biomass flux = {biomass_flux:.6g}")

    return status, selected, biomass_flux, biomass_reaction


def main() -> None:
    import shutil

    if OUTDIR.exists():
        shutil.rmtree(OUTDIR)  # stale leftover from an interrupted previous run
    OUTDIR.mkdir(parents=True, exist_ok=True)
    graph = load_mouse_gem()
    print(f"Full Mouse-GEM: {graph.nv} metabolites, {graph.ne} reactions")

    wide_scores = pd.read_csv(GENE_SCORES_PATH, index_col=0)
    gene_scores, low_thresh, high_thresh = build_pooled_gene_scores(wide_scores)
    gene_scores.rename("score").to_csv(OUTDIR / "gene_scores.csv")
    print(f"Pooled thresholds ({LOW_QUANTILE}/{HIGH_QUANTILE} quantiles): low<={low_thresh:.3g}, high>={high_thresh:.3g}")
    print(f"Gene scores: {int((gene_scores==1).sum())} high, {int((gene_scores==0).sum())} neutral, {int((gene_scores==-1).sum())} low")

    reaction_scores = score_reactions(graph, gene_scores)
    pd.Series(reaction_scores, name="score").to_csv(OUTDIR / "reaction_scores.csv")
    n_pos = sum(1 for v in reaction_scores.values() if v > 0)
    n_neg = sum(1 for v in reaction_scores.values() if v < 0)
    print(f"Reaction scores: {len(reaction_scores)} scored reactions ({n_pos} positive, {n_neg} negative) "
          f"out of {graph.ne} total")

    status, selected, biomass_flux, biomass_reaction = extract_core_network(graph, reaction_scores)
    n_selected = int(selected.sum())
    print(f"Gap-filling solve status: {status}; {biomass_reaction} flux = {biomass_flux:.4g} "
          f"(required >= {BIOMASS_FLOOR}); selected {n_selected}/{graph.ne} reactions")

    selected_indices = np.flatnonzero(selected)
    reaction_ids = list(graph.get_attr_from_edges("id"))
    core_reaction_ids = [reaction_ids[i] for i in selected_indices]

    # Count touched metabolites directly from the full graph's edge
    # endpoints (membership only, not coefficients -- unaffected by the
    # subgraph-extraction coefficient bug since we never extract a subgraph).
    touched_metabolites: set = set()
    for i in selected_indices:
        s, t = graph.get_edge(int(i))
        touched_metabolites |= set(s) | set(t)

    n_evidence_kept = sum(1 for rid in core_reaction_ids if reaction_scores.get(rid, 0.0) != 0)
    n_gap_filled = n_selected - n_evidence_kept

    core_df = pd.DataFrame({
        "reaction_id": core_reaction_ids,
        "evidence_supported": [reaction_scores.get(rid, 0.0) != 0 for rid in core_reaction_ids],
        "score": [reaction_scores.get(rid, 0.0) for rid in core_reaction_ids],
    })
    core_df.to_csv(OUTDIR / "core_reactions.csv", index=False)

    summary_lines = [
        f"full_model_reactions: {graph.ne}",
        f"full_model_metabolites: {graph.nv}",
        f"scored_reactions: {len(reaction_scores)} ({n_pos} positive, {n_neg} negative)",
        f"solver_status: {status}",
        f"biomass_reaction: {biomass_reaction}",
        f"biomass_flux: {biomass_flux:.6g} (floor {BIOMASS_FLOOR})",
        f"core_network_reactions: {n_selected}",
        f"  evidence-supported (kept because scored): {n_evidence_kept}",
        f"  gap-filled (kept with no direct evidence, structurally required): {n_gap_filled}",
        f"core_network_metabolites: {len(touched_metabolites)}",
        f"lambda_reg: {LAMBDA_REG}, eps: {EPS}",
    ]
    (OUTDIR / "build_summary.txt").write_text("\n".join(summary_lines) + "\n")
    print("\n".join(summary_lines))

    final_dir = RESULT_ROOT / f"lr{LAMBDA_REG}_{n_selected}reactions"
    if final_dir.exists():
        import shutil

        shutil.rmtree(final_dir)
    OUTDIR.rename(final_dir)
    print(f"\nSaved to: {final_dir}")


if __name__ == "__main__":
    main()
