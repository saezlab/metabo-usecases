"""Faithful-mechanism reproduction of the colleague's notebook Option 2
(notebooks/model-example.ipynb, cells 19-27: "qualitative metabolomics on
the full Mouse-GEM") -- but WITHOUT forcing exchange-reaction direction,
and using liver tissue metabolomics instead of plasma.

Why not direction, and why not plasma (recap of the conversation this
builds on): the notebook's qualitative_bounds approach (falling medium
glucose -> force uptake, rising medium lactate -> force secretion) only
works because a cultured cell's medium IS that cell's own net exchange --
nothing else in the dish is competing for it. Whole-body blood plasma
around an in vivo liver is not: a rising plasma metabolite can reflect
ANY organ's net secretion, not specifically the liver's, and can be the
opposite of what the liver itself is doing (e.g. plasma lactate can rise
from muscle release while the liver is simultaneously taking lactate up
for gluconeogenesis via the Cori cycle). A first attempt at direction-
forcing on plasma data ran into exactly this kind of infeasibility.

So this script keeps the notebook's actual OPTIMIZATION MECHANISM (plain
MultiSampleFBA, lambda_reg=0, minimize edge_has_flux.sum(), single
scenario, explicit post-hoc verification) but swaps WHAT is required:
instead of forcing a specific flux SIGN on a handful of plasma exchanges,
it requires that EVERY sufficiently-abundant liver tissue metabolite has
AT LEAST ONE of its incident reactions carrying nonzero flux (either
sign) -- the literal "a network that contains reactions that involve
those metabolites" from the colleague's original description, realized
as a flux-feasibility requirement rather than mere graph adjacency (which
is what Script 14's "touch evidence" was -- adjacency alone, no
requirement that flux actually flows).

Tissue (not plasma) metabolites also sidestep the attribution problem
entirely: these are measured inside the liver itself, not in a shared
whole-body pool.

Candidate metabolites, two filters applied before matching:
  1. Low abundance: median log2(x+1) score (Script 07's `metabolite`
     modality, all 16 states) below the pooled 0.25 quantile is dropped
     (171 -> 128 metabolites) -- the project's existing low/high-quantile
     convention (Scripts 01/02), applied here for the first time to this
     modality (Script 07 itself only filters by replicate completeness,
     which drops 0 of 171 tissue metabolites -- there was no abundance
     filter before this script).
  2. Generic cofactors: CoA, ATP, ADP, AMP, NAD(P)(H), FAD, UDP touch
     hundreds-to over a thousand reactions each (CoA alone: 1,249) --
     requiring "at least one reaction active" for these is nearly
     vacuous (trivially satisfied by any reasonable network) and would
     massively inflate the reaction union for no real evidence value.
     Carnitine (394 reactions) is deliberately KEPT despite its high
     degree -- it is fatty-acid-oxidation-specific, directly relevant to
     this starvation dataset, not a universal currency metabolite.
  128 -> 118 after cofactor exclusion -> 84 matched to a Mouse-GEM vertex
  (any compartment) by the same normalized-name matching as Script 14.

FVA precheck (same reasoning as Script 14's first attempt): before
building the "at least one active" requirement, candidate reactions'
structurally achievable flux range (Sv=0, every other reaction at Mouse-
GEM's own default bounds) is checked. A metabolite whose ALL incident
reactions are structurally blocked (FVA range exactly {0}, e.g. a dead-
end reaction with no other producer/consumer anywhere in the model) is
dropped and logged -- its coverage requirement would otherwise make the
whole MILP infeasible, not just fail for that one metabolite. Checking
is scoped to metabolites with <=10 candidate reactions (2 LPs/reaction on
the full 12,987-reaction stoichiometric matrix doesn't scale to the whole
~3k-reaction candidate union -- timed at ~1,000s for that many); a
metabolite with dozens of incident reactions essentially never has every
single one blocked, so the real risk is concentrated exactly where this
checks. Confirmed on a dry run: the known-blocked cases from Script 14's
plasma analysis (glutarate, malonate) are both in this low-degree group.

Mechanism for "at least one active" (CORNETO's NonZeroIndicator, not the
plain Indicator the notebook/Script 02/14 otherwise use): plain Indicator
is one-directional (active=1 only ALLOWS nonzero flux, doesn't require
it) -- forcing it to 1 would be vacuously satisfiable by zero flux.
NonZeroIndicator builds a true biconditional (I_pos/I_neg binaries that
require flux >= tolerance or <= -tolerance respectively when set). For
each metabolite's candidate reactions, summing (I_pos + I_neg) across
just that metabolite's own reactions and requiring the sum >= 1 directly
encodes "at least one of these reactions truly carries flux, in either
direction" -- no bulk MultiSampleFBA() "edge_has_flux" indicator used
for this part, only for the overall minimize-network-size objective.

No biomass/ATP-maintenance floor is required here, matching the
notebook's own Option 2 scope exactly (cell 27: "It has not been asked to
support growth, ATP maintenance, or every function of a living cell.").
The 84 metabolite coverage requirements are themselves what forces
nontrivial network activity.

Outputs (corneto_starvation/result/option2_tissue_coverage_network/):
  candidate_metabolites.csv        which tissue metabolites passed both filters
  metabolite_vertex_matches.csv    metabolite name -> matched Mouse-GEM vertices
  unmatched_metabolites.txt        candidates with no Mouse-GEM vertex match
  metabolite_fva_blocked.txt       metabolites dropped entirely (all candidate
                                    reactions structurally blocked)
  selected_reactions.csv           final network: reaction id, name
  metabolite_coverage.csv          per-metabolite: satisfied?, which reaction(s) active
  build_summary.txt                counts, solver status, verification results
"""

from __future__ import annotations

import re
import sys
import time
from pathlib import Path

PKG_SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(PKG_SRC))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.optimize import linprog  # noqa: E402
from corneto.methods.fba import MultiSampleFBA  # noqa: E402

from corneto_mcmt import load_mouse_gem  # noqa: E402
from corneto_mcmt.model_example import plot_selected_network  # noqa: E402

PKG_ROOT = Path(__file__).resolve().parents[1]
TISSUE_SCORES_PATH = PKG_ROOT / "result" / "metabolite_scores" / "metabolite" / "all_states_scores.csv"

LOW_ABUNDANCE_QUANTILE = 0.25  # same convention as Scripts 01/02's proteomics thresholds
NONZERO_TOLERANCE = 1e-3  # minimum |flux| to count as "this reaction is active"
MAX_SECONDS = 1200
MIP_GAP = 0.05  # matches the notebook's own solve() call exactly

# First attempt (no biomass floor, matching the notebook's literal Option 2
# scope) gave a 191-reaction network where 59% of reactions were plain
# transport/exchange -- the cheapest way to satisfy "at least one reaction
# touches this metabolite" under a pure minimize-count objective, since a
# transport reaction alone is trivially mass-balanced and doesn't drag in
# real metabolic chemistry. Unlike Script 02/14, which anchor to a
# functioning cell via a biomass/ATP-maintenance floor, Option 2's literal
# formulation has no such anchor -- so there's nothing stopping the solver
# from "passing" coverage requirements as cheaply as possible. Re-adding the
# same floor Script 02/14 use, to see whether forcing the network to also
# support real cellular function pulls in more genuine metabolism instead of
# relying on transport shortcuts.
REQUIRE_BIOMASS_FLOOR = True
BIOMASS_FLOOR = 0.01

OUTDIR = PKG_ROOT / "result" / (
    "option2_tissue_coverage_network_biomass" if REQUIRE_BIOMASS_FLOOR else "option2_tissue_coverage_network"
)

# Universal currency metabolites: touch hundreds-to-thousands of reactions
# each, so "at least one active" is nearly vacuous for them. Carnitine is
# deliberately NOT excluded (fatty-acid-oxidation-specific, relevant here).
COFACTOR_EXCLUDE = {"CoA", "ATP", "ADP", "AMP", "NAD+", "NADH", "NADP+", "NADPH", "FAD", "UDP"}

# Same normalized-name matching as Script 14.
SYNONYMS = {
    "Ala": "alanine", "Arg": "arginine", "Asn": "asparagine", "Asp": "aspartate",
    "Gln": "glutamine", "Glu": "glutamate", "Gly": "glycine", "His": "histidine",
    "Ile": "isoleucine", "Leu": "leucine", "Lys": "lysine", "Met": "methionine",
    "Phe": "phenylalanine", "Pro": "proline", "Ser": "serine", "Thr": "threonine",
    "Trp": "tryptophan", "Tyr": "tyrosine", "Val": "valine",
    "beta-Ala": "beta-alanine",
    "Lactate": "L-lactate",
    "Ornithine": "ornithine",
    "Putrescine(1,4-Butanediamine)": "putrescine",
    "3-Indoxyl sulfate": "indoxyl sulfate",
    "alpha-Aminoadipate": "2-aminoadipate",
}
_STEREO_PREFIXES = ("l-", "dl-")


def normalize(name: str) -> str:
    name = name.lower().strip()
    name = name.replace("(", " ").replace(")", " ")
    name = re.sub(r"[,']", "", name)
    name = re.sub(r"\s+", " ", name).strip()
    changed = True
    while changed:
        changed = False
        for prefix in _STEREO_PREFIXES:
            if name.startswith(prefix):
                name = name[len(prefix):]
                changed = True
    return name.strip()


def select_candidate_metabolites() -> pd.DataFrame:
    scores = pd.read_csv(TISSUE_SCORES_PATH, index_col=0)
    median_abundance = scores.median(axis=1, skipna=True)
    cutoff = float(np.nanquantile(median_abundance, LOW_ABUNDANCE_QUANTILE))
    rows = []
    for name, value in median_abundance.items():
        is_cofactor = name in COFACTOR_EXCLUDE
        passes_abundance = value >= cutoff
        rows.append({
            "metabolite_name": name, "median_log2_abundance": value,
            "passes_abundance_filter": passes_abundance, "is_excluded_cofactor": is_cofactor,
            "candidate": passes_abundance and not is_cofactor,
        })
    df = pd.DataFrame(rows)
    print(f"Tissue metabolites: {len(df)} total, {df['passes_abundance_filter'].sum()} pass "
          f"abundance >= {LOW_ABUNDANCE_QUANTILE} quantile ({cutoff:.4g}), "
          f"{df['candidate'].sum()} candidates after cofactor exclusion")
    return df


def build_vertex_name_index(model) -> dict[str, list[str]]:
    index: dict[str, list[str]] = {}
    for v in model.V:
        name = model.get_attr_vertex(v).get("name")
        if name:
            index.setdefault(normalize(str(name)), []).append(v)
    return index


def match_metabolite_vertices(
    metabolite_names: list[str], vertex_index: dict[str, list[str]],
) -> tuple[dict[str, list[str]], list[str]]:
    matches: dict[str, list[str]] = {}
    unmatched: list[str] = []
    for name in metabolite_names:
        lookup_name = SYNONYMS.get(name, name)
        vertices = vertex_index.get(normalize(lookup_name))
        if vertices:
            matches[name] = vertices
        else:
            unmatched.append(name)
    return matches, unmatched


def build_candidate_reactions(model, metabolite_vertices: dict[str, list[str]]) -> dict[str, np.ndarray]:
    """Metabolite name -> array of reaction indices incident to any of its vertices."""
    vertex_row = {v: i for i, v in enumerate(model.V)}
    S = model.vertex_incidence_matrix(values=True, sparse=True)
    per_metabolite: dict[str, np.ndarray] = {}
    for name, vertices in metabolite_vertices.items():
        rows = [vertex_row[v] for v in vertices if v in vertex_row]
        cols = np.asarray((S[rows, :] != 0).sum(axis=0)).ravel().nonzero()[0]
        per_metabolite[name] = cols
    return per_metabolite


def compute_fva_blocked(model, reaction_indices: np.ndarray) -> np.ndarray:
    """Boolean mask (same length as reaction_indices): True if that reaction's
    FVA range is exactly {0} given Sv=0 and every other reaction at default
    bounds -- structurally dead, independent of any requirement we'd add."""
    n_reactions = model.ne
    default_lb = np.asarray(model.get_attr_from_edges("default_lb"), dtype=float)
    default_ub = np.asarray(model.get_attr_from_edges("default_ub"), dtype=float)
    S = model.vertex_incidence_matrix(values=True, sparse=True)
    n_vertices = S.shape[0]
    blocked = np.zeros(len(reaction_indices), dtype=bool)
    for i, idx in enumerate(reaction_indices):
        c = np.zeros(n_reactions)
        c[idx] = 1.0
        lo = linprog(c=c, A_eq=S, b_eq=np.zeros(n_vertices), bounds=list(zip(default_lb, default_ub)), method="highs")
        hi = linprog(c=-c, A_eq=S, b_eq=np.zeros(n_vertices), bounds=list(zip(default_lb, default_ub)), method="highs")
        blocked[i] = abs(lo.fun) < 1e-9 and abs(hi.fun) < 1e-9
    return blocked


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    model = load_mouse_gem()
    print(f"Full Mouse-GEM: {model.nv} metabolites, {model.ne} reactions")

    candidates_df = select_candidate_metabolites()
    candidates_df.to_csv(OUTDIR / "candidate_metabolites.csv", index=False)
    candidate_names = candidates_df.loc[candidates_df["candidate"], "metabolite_name"].tolist()

    vertex_index = build_vertex_name_index(model)
    metabolite_vertices, unmatched = match_metabolite_vertices(candidate_names, vertex_index)
    (OUTDIR / "unmatched_metabolites.txt").write_text("\n".join(sorted(unmatched)) + "\n")
    pd.DataFrame(
        [{"metabolite_name": name, "vertices": ";".join(vs)} for name, vs in metabolite_vertices.items()]
    ).to_csv(OUTDIR / "metabolite_vertex_matches.csv", index=False)
    print(f"Matched to Mouse-GEM vertices: {len(metabolite_vertices)}/{len(candidate_names)}")

    per_metabolite_reactions = build_candidate_reactions(model, metabolite_vertices)
    all_candidate_reactions = sorted(set(np.concatenate(list(per_metabolite_reactions.values())).tolist()))
    print(f"Union of candidate reactions across {len(per_metabolite_reactions)} metabolites: "
          f"{len(all_candidate_reactions)}")

    # FVA (2 LPs/reaction on the full 12,987-reaction stoichiometric matrix)
    # doesn't scale to checking the whole ~3k-reaction union (timed at ~1,000s
    # for that many). The real risk -- ALL of a metabolite's candidates being
    # structurally dead -- only matters for metabolites with FEW candidates;
    # a metabolite with dozens of incident reactions essentially never has
    # every single one blocked. So only metabolites with a small candidate
    # set get checked (confirmed on a dry run: this is exactly where the
    # known-blocked cases, e.g. glutarate/malonate, already showed up).
    FVA_CHECK_MAX_CANDIDATES = 10
    reactions_to_check = sorted(set(
        r for reactions in per_metabolite_reactions.values() if len(reactions) <= FVA_CHECK_MAX_CANDIDATES
        for r in reactions
    ))
    print(f"Running FVA precheck on {len(reactions_to_check)} reactions "
          f"(metabolites with <= {FVA_CHECK_MAX_CANDIDATES} candidates only)...")
    blocked_mask = compute_fva_blocked(model, np.array(reactions_to_check))
    blocked_set = {reactions_to_check[i] for i in np.flatnonzero(blocked_mask)}
    print(f"FVA-blocked (structurally dead) reactions: {len(blocked_set)}/{len(reactions_to_check)}")

    fva_dropped_metabolites = []
    final_per_metabolite: dict[str, np.ndarray] = {}
    for name, reactions in per_metabolite_reactions.items():
        usable = np.array([r for r in reactions if r not in blocked_set])
        if usable.size == 0:
            fva_dropped_metabolites.append(name)
            continue
        final_per_metabolite[name] = usable
    (OUTDIR / "metabolite_fva_blocked.txt").write_text("\n".join(sorted(fva_dropped_metabolites)) + "\n")
    print(f"Metabolites dropped (all candidate reactions FVA-blocked): {len(fva_dropped_metabolites)}")
    print(f"Metabolites proceeding to the solve: {len(final_per_metabolite)}")

    final_reaction_union = np.array(sorted(set(np.concatenate(list(final_per_metabolite.values())).tolist())))
    reaction_ids = list(model.get_attr_from_edges("id"))

    reaction_bounds = {}
    if REQUIRE_BIOMASS_FLOOR:
        graph_attrs = model.get_graph_attributes()
        biomass_reaction = next(iter(graph_attrs["fba_objective"]))
        default_ub = list(model.get_attr_from_edges("default_ub"))
        biomass_upper = float(default_ub[reaction_ids.index(biomass_reaction)])
        reaction_bounds[biomass_reaction] = (BIOMASS_FLOOR, biomass_upper)
        print(f"Requiring biomass/ATP-maintenance floor: {biomass_reaction} >= {BIOMASS_FLOOR}")

    method = MultiSampleFBA(lambda_reg=0.0, beta_reg=0.0)
    problem = method.build(model, reaction_bounds=reaction_bounds)
    backend = method.backend
    problem += backend.NonZeroIndicator(problem.expr.flow, indexes=final_reaction_union, tolerance=NONZERO_TOLERANCE)
    active = problem.expr._flow_ipos + problem.expr._flow_ineg

    position_of = {reaction_idx: pos for pos, reaction_idx in enumerate(final_reaction_union)}
    for name, reactions in final_per_metabolite.items():
        local_positions = [position_of[r] for r in reactions]
        problem += active[local_positions].sum() >= 1

    problem.add_objective(problem.expr.edge_has_flux.sum(), name="network_size")

    print(f"\nSolving the single-scenario MILP on the full {model.ne}-reaction Mouse-GEM "
          f"(max {MAX_SECONDS}s, MIPGap={MIP_GAP})...")
    start = time.time()
    solution = problem.solve(
        solver="GUROBI", verbosity=1, max_seconds=MAX_SECONDS, MIPGap=MIP_GAP, IntFeasTol=1e-9,
    )
    elapsed = time.time() - start
    status = str(solution.status).lower()
    print(f"MILP solve status: {status} ({elapsed:.1f}s)")
    if status not in ("optimal", "optimal_inaccurate", "user_limit"):
        raise RuntimeError(f"Network extraction did not finish: {solution.status}")

    witness_flux = np.asarray(problem.expr.flow.value, dtype=float).ravel()
    selected = np.asarray(problem.expr.edge_has_flux.value).ravel() > 0.5
    n_mismatch = int(np.sum(selected != (np.abs(witness_flux) > 1e-6)))
    if n_mismatch:
        print(f"WARNING: {n_mismatch} reactions have edge_has_flux inconsistent with their witness flux "
              f"(expected with a time-limited, non-optimal MILP solution)")

    S = model.vertex_incidence_matrix(values=True, sparse=True)
    mass_balance_violation = float(np.max(np.abs(S @ witness_flux)))
    print(f"Mass balance check: max |S @ flux| = {mass_balance_violation:.3g} (should be ~0)")

    coverage_rows = []
    n_satisfied = 0
    for name, reactions in final_per_metabolite.items():
        active_reactions = [reaction_ids[r] for r in reactions if abs(witness_flux[r]) > NONZERO_TOLERANCE / 2]
        satisfied = len(active_reactions) > 0
        n_satisfied += int(satisfied)
        coverage_rows.append({
            "metabolite_name": name, "n_candidate_reactions": len(reactions),
            "satisfied": satisfied, "active_reactions": ";".join(active_reactions),
        })
    pd.DataFrame(coverage_rows).to_csv(OUTDIR / "metabolite_coverage.csv", index=False)
    print(f"Metabolite coverage satisfied: {n_satisfied}/{len(final_per_metabolite)}")

    selected_indices = np.flatnonzero(selected)
    selected_df = pd.DataFrame([
        {"reaction": reaction_ids[i], "name": model.get_attr_edge(int(i)).get("name") or "exchange / transport"}
        for i in selected_indices
    ]).set_index("reaction")
    selected_df.to_csv(OUTDIR / "selected_reactions.csv")
    print(f"Selected network: {len(selected_df)}/{model.ne} reactions")

    svg = plot_selected_network(model, selected)
    (OUTDIR / "selected_network.svg").write_text(svg)

    summary = (
        f"Mouse-GEM: {model.ne} reactions total\n"
        f"Biomass/ATP-maintenance floor required: {REQUIRE_BIOMASS_FLOOR}"
        + (f" ({reaction_bounds})\n" if REQUIRE_BIOMASS_FLOOR else "\n")
        + f"Candidate tissue metabolites (abundance + cofactor filters): {len(candidate_names)}\n"
        f"Matched to Mouse-GEM vertices: {len(metabolite_vertices)}\n"
        f"Dropped (all reactions FVA-blocked): {len(fva_dropped_metabolites)}\n"
        f"Metabolites with an active coverage requirement: {len(final_per_metabolite)}\n"
        f"MILP solve status: {status} ({elapsed:.1f}s, max {MAX_SECONDS}s, MIPGap={MIP_GAP})\n"
        f"Mass balance violation (max |S@flux|): {mass_balance_violation:.3g}\n"
        f"edge_has_flux/witness_flux mismatches: {n_mismatch}\n"
        f"Metabolite coverage satisfied: {n_satisfied}/{len(final_per_metabolite)}\n"
        f"Selected network: {len(selected_df)}/{model.ne} reactions\n"
    )
    (OUTDIR / "build_summary.txt").write_text(summary)
    print("\n" + summary)


if __name__ == "__main__":
    main()
