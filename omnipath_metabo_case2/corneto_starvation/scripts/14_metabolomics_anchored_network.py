"""Build an initial Mouse-GEM reaction network anchored on measured
metabolite coverage plus proteomics evidence -- "Example option 2:
qualitative metabolomics on the full Mouse-GEM" from the colleague's
corneto-metabo-mct notebook update (notebooks/model-example.ipynb, cells
19-27), adapted to this project's own starvation time course.

Why this script does NOT force exchange-reaction flux direction (the
notebook's own glucose/lactate toy example does: falling glucose -> forced
uptake, rising lactate -> forced secretion):

  The notebook's Option 2 measures metabolite concentration in the medium
  immediately surrounding CULTURED CELLS -- there, a rising/falling medium
  metabolite is essentially that cell population's own net secretion/uptake,
  since nothing else is exchanging with that medium.

  Our data is whole-body blood PLASMA around an in vivo liver. A plasma
  metabolite's time trend is the net effect of EVERY organ exchanging with
  that pool, not just the liver's. Concretely: during fasting, plasma
  lactate can RISE because muscle releases it, while the liver is
  simultaneously TAKING UP lactate as a gluconeogenic precursor (the Cori
  cycle) -- the exact opposite of the notebook's "rising -> forced
  secretion" rule. Forcing exchange-reaction sign from whole-body plasma
  trends is not just under-determined, it can be systematically backwards
  for exactly the metabolites (lactate, alanine, glycerol...) this
  starvation dataset cares most about. (A first attempt at this script did
  force directions this way; it ran into real infeasibilities from this
  exact issue before the design was corrected.)

  So instead this script uses the weaker, not direction-dependent part of
  the colleague's own description of Option 2: "select a metabolic network
  that **contains reactions that involve those metabolites**. That can be
  combined with transcriptomics or proteomics to [...] give a final
  evidence for each reaction." Coverage ("this metabolite's reactions
  should be in the network") does not have the attribution problem that
  forcing a flux sign does -- it doesn't matter which organ's exchange
  moved plasma lactate, the liver's own lactate-touching reactions are
  still relevant to include. Tissue (intracellular, liver-specific)
  metabolite changes have no attribution problem either way and are
  included on equal footing.

Methodology:
  1. "Variable" metabolites = any plasma OR tissue metabolite (Script 07's
     adjacent_direction.csv, both modalities) with at least one significant
     Up/Down call across the 14 non-T0 (genotype, timepoint) states. T0 is
     excluded (no previous timepoint, FC=NaN).
  2. Each variable metabolite's name is matched to every Mouse-GEM vertex
     sharing that name, in ANY compartment (not just extracellular) --
     capturing the metabolite's exchange reaction, its transporters, and
     every reaction it participates in once inside the cell. All reactions
     incident to any of those vertices are "touch evidence" for that
     reaction.
  3. Combined per-reaction score (the colleague's "final evidence for each
     reaction"): Script 02's pooled proteomics GPR score (+1/0/-1), with
     touch evidence overriding to +1 wherever a reaction touches a variable
     metabolite -- metabolomics coverage takes priority over a
     proteomics-only call, proteomics still disambiguates the reactions
     neither modality touches.
  4. Single-sample MultiSampleIMAT (exactly Script 02's mechanism:
     lambda_reg-penalized reaction selection, gap-filled to a required
     biomass/ATP-maintenance floor) selects the smallest reaction set
     consistent with that combined score and a working steady state.
     Single-sample, not Script 02-first-attempt's 14-sample joint union --
     that joint MILP never got past its own LP relaxation in 1200s; this is
     back to the same scale as Script 02's own already-converged run.

Outputs (corneto_starvation/result/metabolomics_anchored_network/):
  variable_metabolites.csv       which plasma/tissue metabolites qualified, and why
  metabolite_vertex_matches.csv  metabolite name -> matched Mouse-GEM vertices
  unmatched_metabolites.txt      variable metabolites with no Mouse-GEM vertex match
  combined_scores.csv            final per-reaction evidence score fed to MultiSampleIMAT
  selected_reactions.csv         final network: reaction id, name, combined score,
                                  proteomics score, touching variable metabolite(s)
  selected_network.svg           CORNETO plot of the final reaction set
  build_summary.txt              counts, solver status, timing
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
from corneto.methods.fba import MultiSampleFBA  # noqa: E402
from corneto.methods.imat import MultiSampleIMAT  # noqa: E402

from corneto_mcmt import load_mouse_gem  # noqa: E402
from corneto_mcmt.model_example import plot_selected_network  # noqa: E402

PKG_ROOT = Path(__file__).resolve().parents[1]
PLASMA_DIRECTION_PATH = PKG_ROOT / "result" / "metabolite_scores" / "plasma_metabolite" / "adjacent_direction.csv"
TISSUE_DIRECTION_PATH = PKG_ROOT / "result" / "metabolite_scores" / "metabolite" / "adjacent_direction.csv"
PROTEOMICS_SCORES_PATH = PKG_ROOT / "result" / "core_network_extraction" / "lr0.5_1492reactions" / "reaction_scores.csv"
STATE_ORDER = ("2h", "4h", "6h", "8h", "12h", "16h", "24h")  # T0 excluded: no previous timepoint
GENOTYPES = ("WT", "ob_ob")
STATE_COLUMNS = [f"{g}_{t}" for g in GENOTYPES for t in STATE_ORDER]

# lr0.5 (Script 02's own value) gave 4,807 reactions (too big for Script 16's
# downstream temporal fit -- never found a feasible solution in 75+ minutes).
# lr2.0 overshot badly: 388 reactions, and lost almost all metabolomics
# coverage (only 173/4136 touched reactions survived) -- defeats the point
# of this network. Trying lr1.0 as a midpoint between the two.
LAMBDA_REG = 1.0
# Labeled output directory (mirrors Script 02's lr{LAMBDA_REG}_{n}reactions/
# convention) so different lambda_reg attempts don't overwrite each other.
OUTDIR = PKG_ROOT / "result" / f"metabolomics_anchored_network_lr{LAMBDA_REG}"
EPS = 1e-2
BIOMASS_FLOOR = 0.01
MAX_SECONDS = 1200
_FEASIBILITY_TOL = 1e-6

# Metabolite panel name -> Mouse-GEM metabolite 'name' attribute, for cases
# the automated normalized-name match gets wrong or misses: three-letter
# amino acid codes (GEM uses full names), and a few explicit stereoisomer
# choices (GEM lists both D- and L-/unprefixed forms; physiological
# metabolomics measures the natural L-/unprefixed form unless noted).
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

# Only l-/dl- are treated as interchangeable with the unprefixed name; d- is
# kept distinct (a different, usually non-physiological stereoisomer).
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


def variable_metabolites(direction_path: Path, modality: str) -> pd.DataFrame:
    """Metabolites with >=1 significant Up/Down call across the 14 states."""
    direction = pd.read_csv(direction_path, index_col=0)
    calls = direction[STATE_COLUMNS]
    is_variable = (calls == "Up").any(axis=1) | (calls == "Down").any(axis=1)
    n_up = (calls == "Up").sum(axis=1)
    n_down = (calls == "Down").sum(axis=1)
    return pd.DataFrame({
        "metabolite_name": direction.index, "modality": modality,
        "n_states_up": n_up.to_numpy(), "n_states_down": n_down.to_numpy(),
    })[is_variable.to_numpy()]


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


def build_touch_evidence(
    model, metabolite_vertices: dict[str, list[str]],
) -> tuple[np.ndarray, dict[str, list[str]]]:
    """Boolean per-reaction touch mask, plus reaction -> metabolite names
    touching it (for reporting)."""
    reaction_ids = list(model.get_attr_from_edges("id"))
    vertex_row = {v: i for i, v in enumerate(model.V)}
    S = model.vertex_incidence_matrix(values=True, sparse=True)

    touched = np.zeros(len(reaction_ids), dtype=bool)
    reaction_to_metabolites: dict[str, list[str]] = {}
    for metabolite_name, vertices in metabolite_vertices.items():
        rows = [vertex_row[v] for v in vertices if v in vertex_row]
        if not rows:
            continue
        reaction_cols = np.asarray((S[rows, :] != 0).sum(axis=0)).ravel() > 0
        touched |= reaction_cols
        for idx in np.flatnonzero(reaction_cols):
            reaction_to_metabolites.setdefault(reaction_ids[idx], []).append(metabolite_name)
    return touched, reaction_to_metabolites


def build_combined_scores(
    reaction_ids: list[str], proteomics_scores: pd.Series, touched: np.ndarray,
) -> pd.Series:
    """Colleague's "final evidence for each reaction": proteomics GPR score,
    with metabolomics touch evidence overriding to +1 (coverage takes
    priority; proteomics still disambiguates everything metabolomics
    doesn't touch)."""
    combined = proteomics_scores.reindex(reaction_ids).fillna(0.0)
    combined[np.asarray(reaction_ids)[touched]] = 1.0
    return combined


def extract_network(model, combined_scores: pd.Series):
    graph_attrs = model.get_graph_attributes()
    biomass_reaction = next(iter(graph_attrs["fba_objective"]))
    reaction_ids = list(model.get_attr_from_edges("id"))
    default_ub = list(model.get_attr_from_edges("default_ub"))
    biomass_upper = float(default_ub[reaction_ids.index(biomass_reaction)])
    biomass_bounds = {biomass_reaction: (BIOMASS_FLOOR, biomass_upper)}

    method = MultiSampleIMAT(lambda_reg=LAMBDA_REG, beta_reg=0.0, eps=EPS, scale=False, use_bigm_constraints=True)
    problem = method.build_many(
        model, reaction_scores={"combined": combined_scores.to_dict()},
        reaction_bounds={"combined": biomass_bounds},
    )
    start = time.time()
    result = problem.solve(solver="GUROBI", verbosity=1, max_seconds=MAX_SECONDS, MIPFocus=1, IntFeasTol=1e-9)
    elapsed = time.time() - start
    status = str(result.status).lower()
    print(f"MILP solve status: {status} ({elapsed:.1f}s)")

    selector = np.asarray(problem.expr.edge_has_flux.value, dtype=float).reshape(-1)
    selected = np.rint(selector).astype(bool)
    return status, selected, biomass_reaction, biomass_bounds, elapsed


def repair_and_validate(model, reaction_ids, selected, biomass_reaction, biomass_bounds):
    """Re-solve as a plain LP with the selection fixed, exactly Script 02's
    pattern: unambiguous, cleanly mass-balanced flux for exactly that
    reaction set, rather than trusting a possibly time-limited MILP's own
    numeric slack."""
    forced_bounds = dict(biomass_bounds)
    for reaction_id, keep in zip(reaction_ids, selected):
        if not keep:
            forced_bounds[reaction_id] = (0.0, 0.0)
    method = MultiSampleFBA()
    problem = method.build_many(model, objectives={"repair": {}}, reaction_bounds={"repair": forced_bounds})
    result = problem.solve(solver="GUROBI", verbosity=0)
    status = str(result.status).lower()
    if status != "optimal":
        raise RuntimeError(f"Repair LP for the selected reaction set is infeasible: status={status}")
    flux = np.asarray(problem.expr.flow.value, dtype=float).reshape(-1)
    biomass_flux = float(flux[reaction_ids.index(biomass_reaction)])
    if biomass_flux < BIOMASS_FLOOR - _FEASIBILITY_TOL:
        raise RuntimeError(f"Repaired biomass flux {biomass_flux:.6g} is below the required floor {BIOMASS_FLOOR}.")
    return biomass_flux


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    model = load_mouse_gem()
    print(f"Full Mouse-GEM: {model.nv} metabolites, {model.ne} reactions")

    plasma_variable = variable_metabolites(PLASMA_DIRECTION_PATH, "plasma")
    tissue_variable = variable_metabolites(TISSUE_DIRECTION_PATH, "tissue")
    variable = pd.concat([plasma_variable, tissue_variable], ignore_index=True)
    variable.to_csv(OUTDIR / "variable_metabolites.csv", index=False)
    print(f"Variable metabolites: {len(plasma_variable)} plasma, {len(tissue_variable)} tissue")

    vertex_index = build_vertex_name_index(model)
    metabolite_vertices, unmatched = match_metabolite_vertices(list(variable["metabolite_name"]), vertex_index)
    (OUTDIR / "unmatched_metabolites.txt").write_text("\n".join(sorted(unmatched)) + "\n")
    pd.DataFrame(
        [{"metabolite_name": name, "vertices": ";".join(vs)} for name, vs in metabolite_vertices.items()]
    ).to_csv(OUTDIR / "metabolite_vertex_matches.csv", index=False)
    print(f"Matched to Mouse-GEM vertices: {len(metabolite_vertices)}/{len(variable)}")

    touched, reaction_to_metabolites = build_touch_evidence(model, metabolite_vertices)
    print(f"Reactions touching a variable metabolite: {int(touched.sum())}/{model.ne}")

    reaction_ids = list(model.get_attr_from_edges("id"))
    proteomics_scores = pd.read_csv(PROTEOMICS_SCORES_PATH, index_col=0)["score"]
    combined_scores = build_combined_scores(reaction_ids, proteomics_scores, touched)
    combined_scores.rename("combined_score").to_csv(OUTDIR / "combined_scores.csv")
    print(f"Combined scores: {int((combined_scores == 1.0).sum())} positive, "
          f"{int((combined_scores == -1.0).sum())} negative, {int((combined_scores == 0.0).sum())} neutral")

    print(f"\nSolving the single-sample MILP on the full {model.ne}-reaction Mouse-GEM (max {MAX_SECONDS}s)...")
    status, selected, biomass_reaction, biomass_bounds, elapsed = extract_network(model, combined_scores)
    n_selected = int(selected.sum())
    print(f"Selected network: {n_selected}/{model.ne} reactions")

    biomass_flux = repair_and_validate(model, reaction_ids, selected, biomass_reaction, biomass_bounds)
    print(f"Repair LP: validated, {biomass_reaction} flux = {biomass_flux:.4g} (required >= {BIOMASS_FLOOR})")

    selected_indices = np.flatnonzero(selected)
    rows = []
    for idx in selected_indices:
        reaction_id = reaction_ids[idx]
        rows.append({
            "reaction": reaction_id,
            "name": model.get_attr_edge(int(idx)).get("name") or "exchange / transport",
            "combined_score": combined_scores.get(reaction_id, 0.0),
            "proteomics_score": proteomics_scores.get(reaction_id, np.nan),
            "touches_variable_metabolite": ";".join(sorted(reaction_to_metabolites.get(reaction_id, []))),
        })
    selected_df = pd.DataFrame(rows).set_index("reaction")
    selected_df.to_csv(OUTDIR / "selected_reactions.csv")

    svg = plot_selected_network(model, selected)
    (OUTDIR / "selected_network.svg").write_text(svg)

    n_touching = int((selected_df["touches_variable_metabolite"] != "").sum())
    n_proteomics_supported = int((selected_df["proteomics_score"] > 0).sum())
    n_proteomics_contradicted = int((selected_df["proteomics_score"] < 0).sum())
    summary = (
        f"Mouse-GEM: {model.ne} reactions total\n"
        f"Variable metabolites: {len(plasma_variable)} plasma, {len(tissue_variable)} tissue "
        f"({len(metabolite_vertices)}/{len(variable)} matched to a Mouse-GEM vertex)\n"
        f"Reactions touching a variable metabolite: {int(touched.sum())}/{model.ne}\n"
        f"MILP solve status: {status} ({elapsed:.1f}s, max {MAX_SECONDS}s)\n"
        f"Selected network: {n_selected}/{model.ne} reactions\n"
        f"  touching a variable metabolite: {n_touching}\n"
        f"  proteomics-supported (score > 0): {n_proteomics_supported}\n"
        f"  proteomics-contradicted (score < 0): {n_proteomics_contradicted}\n"
        f"  repaired {biomass_reaction} flux: {biomass_flux:.4g} (required >= {BIOMASS_FLOOR})\n"
    )
    (OUTDIR / "build_summary.txt").write_text(summary)
    print("\n" + summary)


if __name__ == "__main__":
    main()
