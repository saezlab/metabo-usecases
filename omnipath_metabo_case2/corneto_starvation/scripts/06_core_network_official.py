"""Mirror the colleague's notebooks/model-example.ipynb workflow exactly,
swapping in our real data in place of his simulated data. Same network
(`build_core_metabolic_network()`, his hand-curated ~58-reaction hepatic
core model -- real human gene symbols in its GPRs, mass-balanced with
proper exchange reactions so it needs no synthetic sinks), same
TemporalMultiSampleIMAT hyperparameters as his notebook, and his own
plotting helpers from `corneto_mcmt.model_example` (not our own).

Per the notebook's own "Use your own data" table:
  Metabolic graph with reaction GPRs and bounds -> build_core_metabolic_network()
  Expression matrix + ordered sample labels    -> our real log2 protein abundance
                                                    (Script 01 output), restricted
                                                    to this network's 86 genes
  Score thresholds + GPR preprocessing          -> ours to prepare (the notebook
                                                    text says so explicitly) --
                                                    mirrored from his own
                                                    DEFAULT_PROTEOMICS_PARAMETERS
                                                    (flux_benchmark.py): global
                                                    pooled (0.50, 0.76) quantiles
                                                    over the gene x timepoint
                                                    matrix, not the per-gene
                                                    percentile our Script 03 used
  Condition -> reaction score mapping            -> corneto.methods.imat.MultiSampleIMAT
                                                    .preprocess(), exactly as the
                                                    notebook does
  Ordered adjacent sample pairs                  -> our real 7 within-genotype pairs

72 of the network's 86 genes are measured in our proteomics data (matched by
uppercasing our mouse Title-case symbols -- e.g. "Pck1" -> "PCK1" -- against
the network's human-style GPR symbols). The other 14 (ENO1, GPI, MPC1, NNT,
NNT_REV, PC, PCK2, SLC16A6/7, SLC2A1, SLC37A4, SLC38A2/4, UQCRC1) are kept in
the displayed gene list with NaN rows, matching how the colleague's own
simulator represents a missing measurement (gray cells, not zero evidence).

Fit twice (WT, ob_ob), each a single 8-timepoint trajectory -- exactly the
shape of the notebook's own single-trajectory example. No cross-genotype
joint fit here. `evaluate_benchmark_prediction`/`evaluation_table` are
skipped: they compare against simulated ground truth, which doesn't exist
for real biology.

Outputs (corneto_starvation/result/core_network/<genotype>/):
  full_model.html            graph structure (plot_full_model)
  input_proteomics.png       measured protein abundance heatmap
  input_reaction_scores.png  GPR-derived reaction score heatmap
  score_counts.csv           high/neutral/low score counts per timepoint
  inference_gallery.html     per-timepoint inferred-network panels
  activity_overview.png      reaction activity heatmap across timepoints
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
from corneto.methods.metabolism import get_unique_genes  # noqa: E402

from corneto_mcmt import TemporalMultiSampleIMAT  # noqa: E402
from corneto_mcmt.flux_benchmark import DEFAULT_PROTEOMICS_PARAMETERS, build_core_metabolic_network  # noqa: E402
from corneto_mcmt.model_example import (  # noqa: E402
    plot_activity_overview,
    plot_full_model,
    plot_inference_gallery,
    plot_input_heatmaps,
    score_count_table,
)

GENE_SCORES_PATH = Path(__file__).resolve().parents[1] / "result" / "gene_scores" / "all_states_gene_scores.csv"
OUTDIR = Path(__file__).resolve().parents[1] / "result" / "core_network"

# Mirrors the notebook exactly: epsilon_fit=0.03, lambda_reg=0.08,
# lambda_temporal=0.16, lambda_flux=0.05, temporal_mode="flux",
# flux_penalty="l1", eps left at the class default (1e-2).
EPSILON_FIT = 0.03
LAMBDA_REG = 0.08
LAMBDA_TEMPORAL = 0.16
LAMBDA_FLUX = 0.05
TEMPORAL_MODE = "flux"
FLUX_PENALTY = "l1"
MAX_SECONDS = 300


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
    """Minimal duck-typed stand-in for FluxBenchmarkDataset.

    plot_input_heatmaps/score_count_table only read .conditions, .genes,
    .reactions, and .proteomics.{state_expression,reaction_score_matrix} --
    this supplies exactly those, without faking any of the simulated-truth
    fields real data doesn't have.
    """

    def __init__(self, genes, reactions, conditions, proteomics):
        self.genes = genes
        self.reactions = reactions
        self.conditions = conditions
        self.proteomics = proteomics


def build_gene_expression_matrix(wide_scores: pd.DataFrame, network_genes: list[str], states: list[str]):
    """Mouse Title-case symbols -> network's human-style uppercase symbols."""
    upper_to_mouse: dict[str, str] = {}
    for symbol in wide_scores.index:
        upper_to_mouse.setdefault(symbol.upper(), symbol)

    matrix = np.full((len(network_genes), len(states)), np.nan)
    for row, gene in enumerate(network_genes):
        mouse_symbol = upper_to_mouse.get(gene)
        if mouse_symbol is None:
            continue
        matrix[row, :] = wide_scores.loc[mouse_symbol, states].to_numpy(dtype=float)
    return matrix


def score_genes(state_expression: np.ndarray, quantiles: tuple[float, float]) -> tuple[np.ndarray, float, float]:
    """Mirror flux_benchmark.py's own thresholding: global pooled quantiles
    over the full gene x timepoint matrix (not per-gene), exactly the logic
    used for DEFAULT_PROTEOMICS_PARAMETERS in his simulator."""
    low_q, high_q = quantiles
    low_thresh, high_thresh = np.nanquantile(state_expression, [low_q, high_q])
    gene_scores = np.zeros_like(state_expression)
    finite = np.isfinite(state_expression)
    gene_scores[finite & (state_expression >= high_thresh)] = 1.0
    gene_scores[finite & (state_expression <= low_thresh)] = -1.0
    gene_scores[~finite] = np.nan
    return gene_scores, float(low_thresh), float(high_thresh)


def gene_scores_to_reaction_scores(graph, network_genes, gene_scores, states):
    import corneto as cn_data

    cdict = {}
    for col, state in enumerate(states):
        cdict[state] = {
            gene: {"value": float(gene_scores[row, col]), "mapping": "none", "role": "expression"}
            for row, gene in enumerate(network_genes)
            if np.isfinite(gene_scores[row, col])
        }
    data = cn_data.data.Data.from_cdict(cdict)
    scorer = MultiSampleIMAT(high_expression_threshold=1.0, low_expression_threshold=-1.0, gpr_field="GPR")
    _, processed = scorer.preprocess(graph, data)

    reaction_scores: dict[str, dict[str, float]] = {}
    for state, sample in processed.samples.items():
        scores = {}
        for feature in sample.features:
            if feature.mapping == "edge" and feature.data.get("role") == "expression" and feature.value is not None:
                scores[feature.id] = float(feature.value)
        reaction_scores[state] = scores
    return reaction_scores


def write_image(image, path: Path) -> None:
    path.write_bytes(image.data)


def write_html(html: str, path: Path) -> None:
    path.write_text(f"<html><body>{html}</body></html>")


def run_genotype(genotype: str, wide_scores: pd.DataFrame, graph, reaction_ids, network_genes) -> None:
    outdir = OUTDIR / genotype
    outdir.mkdir(parents=True, exist_ok=True)

    states = order_timepoints([s for s in wide_scores.columns if s.startswith(genotype)])
    state_expression = build_gene_expression_matrix(wide_scores, network_genes, states)

    gene_scores, low_thresh, high_thresh = score_genes(
        state_expression, DEFAULT_PROTEOMICS_PARAMETERS.threshold_quantiles
    )
    reaction_scores = gene_scores_to_reaction_scores(graph, network_genes, gene_scores, states)

    reaction_score_matrix = np.zeros((len(reaction_ids), len(states)))
    for col, state in enumerate(states):
        for row, reaction_id in enumerate(reaction_ids):
            reaction_score_matrix[row, col] = reaction_scores[state].get(reaction_id, 0.0)

    dataset = _Dataset(
        genes=tuple(network_genes),
        reactions=tuple(reaction_ids),
        conditions=tuple(states),
        proteomics=_Modality(
            state_expression=state_expression,
            thresholds=(low_thresh, high_thresh),
            gene_scores=gene_scores,
            reaction_scores=reaction_scores,
            reaction_score_matrix=reaction_score_matrix,
        ),
    )

    write_html(plot_full_model(graph), outdir / "full_model.html")
    proteomics_img, scores_img = plot_input_heatmaps(dataset)
    write_image(proteomics_img, outdir / "input_proteomics.png")
    write_image(scores_img, outdir / "input_reaction_scores.png")
    score_count_table(dataset).to_csv(outdir / "score_counts.csv")

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
        temporal_pairs=pairs,
        solver="GUROBI",
        max_seconds=MAX_SECONDS,
        accept_feasible_joint_incumbent=True,
        verbosity=1,
    )

    write_html(plot_inference_gallery(graph, result), outdir / "inference_gallery.html")
    write_image(plot_activity_overview(result), outdir / "activity_overview.png")

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
    graph, _tasks = build_core_metabolic_network()
    reaction_ids = tuple(graph.get_attr_from_edges("id"))
    network_genes = sorted(get_unique_genes(graph, gpr_field="GPR"))
    print(f"Core network: {graph.nv} metabolites, {graph.ne} reactions, {len(network_genes)} genes")

    wide_scores = pd.read_csv(GENE_SCORES_PATH, index_col=0)

    for genotype in ("WT", "ob_ob"):
        print(f"=== {genotype} ===")
        run_genotype(genotype, wide_scores, graph, reaction_ids, network_genes)


if __name__ == "__main__":
    main()
