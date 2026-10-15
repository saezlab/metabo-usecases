"""Joint WT + ob_ob temporal analysis on the core network from Script 02.

Unlike Script 03 (which fits each genotype as a separate, independent
TemporalMultiSampleIMAT.fit() call), this script fits BOTH genotypes'
16 states (WT_0h..WT_24h, ob_ob_0h..ob_ob_24h) in a single joint MILP.
temporal_pairs is the union of each genotype's own 7 within-genotype
adjacent-timepoint pairs (14 total) -- there is no WT<->ob_ob pair, since
the two genotypes are not temporally adjacent to each other.

The practical effect: CORNETO's lambda_reg structural-sparsity term (which
penalizes the union of reactions ever selected across every condition
passed to one .fit() call) now operates over the union across BOTH
genotypes together, instead of two separate per-genotype unions. This
matches the project's original plan (before the full-Mouse-GEM pivot),
which called for a single joint fit across both genotypes.

Score thresholds: the two genotypes' expression values are pooled into one
16-state matrix and the (0.50, 0.76) global quantile thresholds are
computed over that pooled matrix, not separately per genotype -- so "high"
and "low" mean the same thing for both genotypes in this run.

IMPORTANT: this is a strictly bigger MILP than Script 03's single-genotype
fit (16 states instead of 8), so it is expected to be at least as hard to
converge, not easier. Run only with a generous time budget.

Outputs (corneto_starvation/result/temporal_core_network_joint/<genotype>/):
  full_model.html, input_proteomics.png, input_reaction_scores.png,
  score_counts.csv, activity_overview.png, run_summary.txt -- same as
  Script 03, split back out per genotype for readability.
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

GENOTYPES = ("WT", "ob_ob")

# Same as Script 06 / Script 03 / the colleague's notebook.
EPSILON_FIT = 0.03
LAMBDA_REG = 0.08
LAMBDA_TEMPORAL = 0.16
LAMBDA_FLUX = 0.05
TEMPORAL_MODE = "flux"
FLUX_PENALTY = "l1"
MAX_SECONDS = 21600  # 6h joint-solve budget: double Script 03's single-genotype 3h attempt, since this problem is ~2x the size
BIOMASS_FLOOR = 0.01
ZERO_BOUND = (0.0, 0.0)

# Labeled per core-network run + time budget, so different attempts don't
# silently overwrite each other's results.
OUTDIR = (
    Path(__file__).resolve().parents[1] / "result" / "temporal_core_network_joint"
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
                and feature.id in core_reaction_ids
            ):
                scores[feature.id] = float(feature.value)
        reaction_scores[state] = scores
    return reaction_scores


def write_image(image, path: Path) -> None:
    path.write_bytes(image.data)


def write_html(html: str, path: Path) -> None:
    path.write_text(f"<html><body>{html}</body></html>")


def write_genotype_outputs(
    genotype: str,
    genotype_states: list[str],
    gene_idx_by_state: dict[str, int],
    network_genes: list[str],
    core_reaction_id_list: list[str],
    all_states: list[str],
    state_expression: np.ndarray,
    gene_scores: np.ndarray,
    reaction_scores: dict[str, dict[str, float]],
    low_thresh: float,
    high_thresh: float,
    result,
) -> None:
    outdir = OUTDIR / genotype
    outdir.mkdir(parents=True, exist_ok=True)

    cols = [all_states.index(s) for s in genotype_states]
    g_state_expression = state_expression[:, cols]
    g_gene_scores = gene_scores[:, cols]

    reaction_score_matrix = np.zeros((len(core_reaction_id_list), len(genotype_states)))
    for col, state in enumerate(genotype_states):
        for row, reaction_id in enumerate(core_reaction_id_list):
            reaction_score_matrix[row, col] = reaction_scores[state].get(reaction_id, 0.0)

    dataset = _Dataset(
        genes=tuple(network_genes),
        reactions=tuple(core_reaction_id_list),
        conditions=tuple(genotype_states),
        proteomics=_Modality(
            state_expression=g_state_expression,
            thresholds=(low_thresh, high_thresh),
            gene_scores=g_gene_scores,
            reaction_scores={s: reaction_scores[s] for s in genotype_states},
            reaction_score_matrix=reaction_score_matrix,
        ),
    )

    proteomics_img, scores_img = plot_input_heatmaps(dataset)
    write_image(proteomics_img, outdir / "input_proteomics.png")
    write_image(scores_img, outdir / "input_reaction_scores.png")
    score_count_table(dataset).to_csv(outdir / "score_counts.csv")

    reaction_index = {rid: i for i, rid in enumerate(result.reactions)}
    condition_index = {c: i for i, c in enumerate(result.conditions)}
    core_rows = [reaction_index[rid] for rid in core_reaction_id_list]
    g_cols = [condition_index[s] for s in genotype_states]
    core_activity = np.asarray(result.activity)[np.ix_(core_rows, g_cols)]
    activity_view = _ActivityView(core_activity, tuple(genotype_states), tuple(core_reaction_id_list))
    write_image(plot_activity_overview(activity_view), outdir / "activity_overview.png")

    # Switches within just this genotype's own adjacent pairs (subset of
    # result.switching_counts, which also contains the other genotype's pairs).
    g_pairs = tuple(zip(genotype_states[:-1], genotype_states[1:]))
    g_switches = sum(result.switching_counts.get(pair, 0) for pair in g_pairs)

    summary_lines = [
        f"genotype: {genotype}",
        f"states: {', '.join(genotype_states)}",
        f"gene coverage: {int(np.isfinite(g_state_expression).any(axis=1).sum())}/{len(network_genes)} network genes measured",
        f"score thresholds (pooled across BOTH genotypes, {DEFAULT_PROTEOMICS_PARAMETERS.threshold_quantiles} quantiles): low<={low_thresh:.3g}, high>={high_thresh:.3g}",
        f"[joint fit] solver_status: {result.solver_status}",
        f"[joint fit] optimality_proven: {result.optimality_proven}",
        f"[joint fit] objective_value: {result.objective_value:.6g}",
        f"[joint fit] total_switches (both genotypes): {result.total_switches}",
        f"switches within {genotype} only: {g_switches}",
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

    OUTDIR.mkdir(parents=True, exist_ok=True)
    core_indices = [i for i, rid in enumerate(all_reaction_ids) if rid in core_reaction_ids]
    viz_graph = graph.edge_subgraph(core_indices)
    write_html(plot_full_model(viz_graph, layout="neato"), OUTDIR / "full_model.html")

    wide_scores = pd.read_csv(GENE_SCORES_PATH, index_col=0)

    per_genotype_states = {g: order_timepoints([s for s in wide_scores.columns if s.startswith(g)]) for g in GENOTYPES}
    all_states = per_genotype_states["WT"] + per_genotype_states["ob_ob"]
    print(f"Joint states ({len(all_states)}): {', '.join(all_states)}")

    state_expression = build_gene_expression_matrix(wide_scores, network_genes, all_states)
    gene_scores, low_thresh, high_thresh = score_genes(
        state_expression, DEFAULT_PROTEOMICS_PARAMETERS.threshold_quantiles
    )
    reaction_scores = gene_scores_to_reaction_scores(graph, network_genes, gene_scores, all_states, core_reaction_ids)

    all_reaction_ids_t = tuple(all_reaction_ids)
    zero_bounds = {rid: ZERO_BOUND for rid in all_reaction_ids_t if rid not in core_reaction_ids}
    per_state_bounds = {biomass_reaction: (BIOMASS_FLOOR, biomass_upper), **zero_bounds}
    reaction_bounds = {state: per_state_bounds for state in all_states}

    # 14 within-genotype adjacent pairs; no cross-genotype pair.
    pairs = tuple(
        pair
        for g in GENOTYPES
        for pair in zip(per_genotype_states[g][:-1], per_genotype_states[g][1:])
    )
    print(f"Temporal pairs ({len(pairs)}): {pairs}")

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

    for genotype in GENOTYPES:
        write_genotype_outputs(
            genotype, per_genotype_states[genotype], {}, network_genes, core_reaction_id_list,
            all_states, state_expression, gene_scores, reaction_scores, low_thresh, high_thresh, result,
        )


if __name__ == "__main__":
    main()
