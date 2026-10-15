"""Per-genotype temporal analysis restricted to the core network from
Script 02 (result/core_network_extraction/core_reactions.csv), reusing the
exact same TemporalMultiSampleIMAT settings and the colleague's own plotting
functions (corneto_mcmt.model_example) as Script 06.

IMPORTANT: this operates on the FULL Mouse-GEM graph throughout, not an
extracted subgraph object. Verified directly (see Script 02's docstring):
CORNETO 1.0.0rc8's edge_subgraph()/extract_subgraph() silently flattens any
non-unit stoichiometric coefficient to +-1 when copying edges, which
corrupted an earlier version of this pipeline (a core network that was
biomass-feasible within the full model came back infeasible once extracted
as its own subgraph object). "Restricted to the core network" is instead
done the way the biomass-feasibility repair in Script 02 already does it
correctly: every reaction NOT in the core set gets bounds forced to
(0.0, 0.0) for every timepoint, on the full, mathematically correct graph.
Gurobi's presolve fixes those variables away, so the actual solved problem
stays small despite the full 12,987-reaction graph being passed in.

Gene symbols are native mouse Title-case (e.g. "Pck1"), matching our
proteomics directly -- no case conversion needed.

Outputs (corneto_starvation/result/temporal_core_network/<genotype>/):
  full_model.html            core network structure (plot_full_model, drawn
                              from a visualization-only edge_subgraph -- safe
                              here since nothing is solved from it)
  input_proteomics.png       measured protein abundance heatmap
  input_reaction_scores.png  GPR-derived reaction score heatmap
  score_counts.csv           high/neutral/low score counts per timepoint
  inference_gallery.html     per-timepoint inferred-network panels
  activity_overview.png      reaction activity heatmap, core reactions only
  run_summary.txt            solver status, objective, switch counts
"""

from __future__ import annotations

import sys
from pathlib import Path

PKG_SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(PKG_SRC))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import corneto as cn  # noqa: E402
from corneto.methods.imat import MultiSampleIMAT  # noqa: E402

from corneto_mcmt import TemporalMultiSampleIMAT, load_mouse_gem  # noqa: E402
from corneto_mcmt.flux_benchmark import DEFAULT_PROTEOMICS_PARAMETERS  # noqa: E402
from corneto_mcmt.model_example import (  # noqa: E402
    plot_activity_overview,
    plot_full_model,
    plot_inference_gallery,
    plot_input_heatmaps,
    score_count_table,
)

# Which labeled Script 02 output to use -- update this to point at whichever
# core_network_extraction/lr{lambda_reg}_{n}reactions/ run you want to fit on.
CORE_NETWORK_RUN = "lr0.5_1492reactions"
CORE_REACTIONS_PATH = (
    Path(__file__).resolve().parents[1] / "result" / "core_network_extraction" / CORE_NETWORK_RUN / "core_reactions.csv"
)
GENE_SCORES_PATH = Path(__file__).resolve().parents[1] / "result" / "gene_scores" / "all_states_gene_scores.csv"

# Same as Script 06 / the colleague's notebook.
EPSILON_FIT = 0.03
LAMBDA_REG = 0.08
LAMBDA_TEMPORAL = 0.16
LAMBDA_FLUX = 0.05
TEMPORAL_MODE = "flux"
FLUX_PENALTY = "l1"
MAX_SECONDS = 10800
BIOMASS_FLOOR = 0.01
ZERO_BOUND = (0.0, 0.0)

# Labeled per core-network run + time budget, so different attempts don't
# silently overwrite each other's results.
OUTDIR = (
    Path(__file__).resolve().parents[1] / "result" / "temporal_core_network"
    / f"{CORE_NETWORK_RUN}_t{MAX_SECONDS}s"
)


def order_timepoints(states: list[str]) -> list[str]:
    return sorted(states, key=lambda s: int(s.rsplit("_", 1)[1].rstrip("h")))


class _Modality:
    def __init__(self, state_expression, thresholds, gene_scores, reaction_scores, reaction_score_matrix):
        self.state_expression = state_expression
        self.thresholds = thresholds
        self.gene_scores = gene_scores
        self.reaction_scores = reaction_scores
        self.reaction_score_matrix = reaction_score_matrix


class _Dataset:
    def __init__(self, genes, reactions, conditions, proteomics):
        self.genes = genes
        self.reactions = reactions
        self.conditions = conditions
        self.proteomics = proteomics


class _ActivityView:
    """Duck-typed stand-in restricting a TemporalIMATResult's activity
    heatmap to just the core reactions, for plot_activity_overview -- the
    real result carries all 12,987 full-model reactions, almost all
    permanently inactive (forced to (0,0) bounds), which would make that
    plot an unreadable, mostly-blank wall."""

    def __init__(self, activity, conditions, reactions):
        self.activity = activity
        self.conditions = conditions
        self.reactions = reactions


def core_reaction_genes(graph, core_reaction_ids: set[str]) -> list[str]:
    import re

    pattern = re.compile(r"[A-Za-z0-9_.\-]+")
    genes: set[str] = set()
    for index, reaction_id in enumerate(graph.get_attr_from_edges("id")):
        if reaction_id not in core_reaction_ids:
            continue
        gpr = graph.get_attr_edge(index).get("GPR") or ""
        if gpr.strip():
            genes.update(token for token in pattern.findall(gpr) if token not in ("and", "or"))
    return sorted(genes)


def build_gene_expression_matrix(wide_scores: pd.DataFrame, network_genes: list[str], states: list[str]):
    matrix = np.full((len(network_genes), len(states)), np.nan)
    for row, gene in enumerate(network_genes):
        if gene not in wide_scores.index:
            continue
        matrix[row, :] = wide_scores.loc[gene, states].to_numpy(dtype=float)
    return matrix


def score_genes(state_expression: np.ndarray, quantiles: tuple[float, float]) -> tuple[np.ndarray, float, float]:
    low_q, high_q = quantiles
    low_thresh, high_thresh = np.nanquantile(state_expression, [low_q, high_q])
    gene_scores = np.zeros_like(state_expression)
    finite = np.isfinite(state_expression)
    gene_scores[finite & (state_expression >= high_thresh)] = 1.0
    gene_scores[finite & (state_expression <= low_thresh)] = -1.0
    gene_scores[~finite] = np.nan
    return gene_scores, float(low_thresh), float(high_thresh)


def gene_scores_to_reaction_scores(graph, network_genes, gene_scores, states, core_reaction_ids: set[str]):
    cdict = {}
    for col, state in enumerate(states):
        cdict[state] = {
            gene: {"value": float(gene_scores[row, col]), "mapping": "none", "role": "expression"}
            for row, gene in enumerate(network_genes)
            if np.isfinite(gene_scores[row, col])
        }
    data = cn.data.Data.from_cdict(cdict)
    scorer = MultiSampleIMAT(high_expression_threshold=1.0, low_expression_threshold=-1.0, gpr_field="GPR")
    _, processed = scorer.preprocess(graph, data)

    reaction_scores: dict[str, dict[str, float]] = {}
    for state, sample in processed.samples.items():
        scores = {}
        for feature in sample.features:
            if (
                feature.mapping == "edge"
                and feature.data.get("role") == "expression"
                and feature.value is not None
                and feature.id in core_reaction_ids  # scope scoring to the core network only
            ):
                scores[feature.id] = float(feature.value)
        reaction_scores[state] = scores
    return reaction_scores


def write_image(image, path: Path) -> None:
    path.write_bytes(image.data)


def write_html(html: str, path: Path) -> None:
    path.write_text(f"<html><body>{html}</body></html>")


def run_genotype(
    genotype: str,
    wide_scores: pd.DataFrame,
    graph,
    core_reaction_ids: set[str],
    core_reaction_id_list: list[str],
    network_genes: list[str],
    biomass_reaction: str,
    biomass_upper: float,
    viz_graph,
) -> None:
    outdir = OUTDIR / genotype
    outdir.mkdir(parents=True, exist_ok=True)

    states = order_timepoints([s for s in wide_scores.columns if s.startswith(genotype)])
    state_expression = build_gene_expression_matrix(wide_scores, network_genes, states)

    gene_scores, low_thresh, high_thresh = score_genes(
        state_expression, DEFAULT_PROTEOMICS_PARAMETERS.threshold_quantiles
    )
    reaction_scores = gene_scores_to_reaction_scores(graph, network_genes, gene_scores, states, core_reaction_ids)

    reaction_score_matrix = np.zeros((len(core_reaction_id_list), len(states)))
    for col, state in enumerate(states):
        for row, reaction_id in enumerate(core_reaction_id_list):
            reaction_score_matrix[row, col] = reaction_scores[state].get(reaction_id, 0.0)

    dataset = _Dataset(
        genes=tuple(network_genes),
        reactions=tuple(core_reaction_id_list),
        conditions=tuple(states),
        proteomics=_Modality(
            state_expression=state_expression,
            thresholds=(low_thresh, high_thresh),
            gene_scores=gene_scores,
            reaction_scores=reaction_scores,
            reaction_score_matrix=reaction_score_matrix,
        ),
    )

    # layout="neato" (force-directed) instead of the default "dot": "dot"
    # hung >15h on this network size; neato renders the same ~4.5k-edge
    # network in ~9 minutes. The result is a dense, largely unreadable
    # "hairball" for the core as a whole -- confirmed by direct inspection --
    # but at least completes, and the periphery (low-degree reactions) is
    # legible.
    write_html(plot_full_model(viz_graph, layout="neato"), outdir / "full_model.html")
    proteomics_img, scores_img = plot_input_heatmaps(dataset)
    write_image(proteomics_img, outdir / "input_proteomics.png")
    write_image(scores_img, outdir / "input_reaction_scores.png")
    score_count_table(dataset).to_csv(outdir / "score_counts.csv")

    # Every non-core reaction is pinned to zero flux at every timepoint, on
    # the full graph; the biomass floor is kept so each timepoint must stay
    # a functioning network within the core set, not an idle one (an
    # all-zero solution is otherwise cheaper than activating a
    # mixed-evidence connector chain -- confirmed by direct testing).
    all_reaction_ids = tuple(graph.get_attr_from_edges("id"))
    zero_bounds = {rid: ZERO_BOUND for rid in all_reaction_ids if rid not in core_reaction_ids}
    per_state_bounds = {biomass_reaction: (BIOMASS_FLOOR, biomass_upper), **zero_bounds}
    reaction_bounds = {state: per_state_bounds for state in states}

    pairs = tuple(zip(states[:-1], states[1:]))
    model = TemporalMultiSampleIMAT(
        epsilon_fit=EPSILON_FIT,
        lambda_reg=LAMBDA_REG,
        lambda_temporal=LAMBDA_TEMPORAL,
        lambda_flux=LAMBDA_FLUX,
        temporal_mode=TEMPORAL_MODE,
        flux_penalty=FLUX_PENALTY,
    )
    result = model.fit(
        graph,
        reaction_scores=reaction_scores,
        reaction_bounds=reaction_bounds,
        temporal_pairs=pairs,
        solver="GUROBI",
        max_seconds=MAX_SECONDS,
        accept_feasible_joint_incumbent=True,
        verbosity=1,
    )

    # Skipped too: also Graphviz-based (draws every reaction active at any
    # timepoint, which on this network could still be a large, slow subset).
    # write_html(plot_inference_gallery(graph, result), outdir / "inference_gallery.html")

    reaction_index = {rid: i for i, rid in enumerate(result.reactions)}
    core_rows = [reaction_index[rid] for rid in core_reaction_id_list]
    core_activity = np.asarray(result.activity)[core_rows, :]
    activity_view = _ActivityView(core_activity, result.conditions, tuple(core_reaction_id_list))
    write_image(plot_activity_overview(activity_view), outdir / "activity_overview.png")

    summary_lines = [
        f"genotype: {genotype}",
        f"states: {', '.join(states)}",
        f"gene coverage: {int(np.isfinite(state_expression).any(axis=1).sum())}/{len(network_genes)} network genes measured",
        f"score thresholds (global pooled {DEFAULT_PROTEOMICS_PARAMETERS.threshold_quantiles} quantiles): low<={low_thresh:.3g}, high>={high_thresh:.3g}",
        f"solver_status: {result.solver_status}",
        f"optimality_proven: {result.optimality_proven}",
        f"objective_value: {result.objective_value:.6g}",
        f"total_switches: {result.total_switches}",
        f"temporal_mode: {TEMPORAL_MODE}",
        f"epsilon_fit: {EPSILON_FIT}, lambda_reg: {LAMBDA_REG}, lambda_temporal: {LAMBDA_TEMPORAL}, lambda_flux: {LAMBDA_FLUX}",
    ]
    (outdir / "run_summary.txt").write_text("\n".join(summary_lines) + "\n")
    print("\n".join(summary_lines))
    print()


def main() -> None:
    graph = load_mouse_gem()
    core_df = pd.read_csv(CORE_REACTIONS_PATH)
    core_reaction_id_list = core_df["reaction_id"].tolist()
    core_reaction_ids = set(core_reaction_id_list)
    network_genes = core_reaction_genes(graph, core_reaction_ids)
    print(f"Core network: {len(core_reaction_ids)} reactions, {len(network_genes)} genes "
          f"(full Mouse-GEM: {graph.ne} reactions)")

    graph_attrs = graph.get_graph_attributes()
    biomass_reaction = next(iter(graph_attrs["fba_objective"]))
    all_reaction_ids = list(graph.get_attr_from_edges("id"))
    biomass_upper = float(graph.get_attr_from_edges("default_ub")[all_reaction_ids.index(biomass_reaction)])

    # Visualization-only subgraph (edge_subgraph is fine here -- nothing is
    # solved from it, only drawn; the stoichiometry-flattening bug only
    # matters for solving).
    core_indices = [i for i, rid in enumerate(all_reaction_ids) if rid in core_reaction_ids]
    viz_graph = graph.edge_subgraph(core_indices)

    wide_scores = pd.read_csv(GENE_SCORES_PATH, index_col=0)

    for genotype in ("WT", "ob_ob"):
        print(f"=== {genotype} ===")
        run_genotype(
            genotype, wide_scores, graph, core_reaction_ids, core_reaction_id_list, network_genes,
            biomass_reaction, biomass_upper, viz_graph,
        )


if __name__ == "__main__":
    main()
