"""Flux-first stoichiometric benchmark simulator for temporal metabolic reconstruction.

Unlike purely statistical gene-expression simulators, this module generates
ground-truth time-varying steady-state flux trajectories ``v*(t)`` on a user-supplied
metabolic network (:class:`corneto.Graph`) first, verifies mass balance, bounds,
internal-L1 parsimony, and thresholded switch counts at ``eps``, and then maps active
reactions through GPR rules to latent enzyme states to simulate paired
**proteomics** (higher consistency, sparser coverage) and **transcriptomics**
(near-complete coverage, higher temporal/decoupling noise) modalities.
"""

from __future__ import annotations

import ast
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from typing import Any

import corneto as cn
import numpy as np
from corneto.data import Data
from corneto.methods.imat import MultiSampleIMAT
from corneto.methods.metabolism import get_genes_from_gpr
from scipy.optimize import linprog
from scipy.sparse import csc_matrix, csr_matrix, eye, hstack, vstack


_FEASIBILITY_TOL = 1e-6


@dataclass(frozen=True)
class MetabolicTask:
    """Definition of a metabolic regime/task on a CORNETO metabolic graph.

    ``bounds`` maps reaction IDs to ``(lb, ub)`` overrides that define the
    regime's active uptake, secretion, and metabolic task requirements when
    generating the ground-truth anchor flux. ``objectives`` optionally maps
    reaction IDs to linear objective weights maximized prior to parsimonious
    L1 flux minimization (pFBA). ``pathway`` optionally labels the regime.
    """

    name: str
    bounds: Mapping[str, tuple[float, float]]
    objectives: Mapping[str, float] | None = None
    opt_fraction: float = 0.99


@dataclass(frozen=True)
class ModalityParameters:
    """Generative observation parameters for one omics modality."""

    signal_separation: float = 2.2
    baseline_sd: float = 0.25
    temporal_noise_sd: float = 0.35
    temporal_ar1: float = 0.0
    replicate_noise_sd: float = 0.20
    decoupling_rate: float = 0.10
    constitutive_gene_rate: float = 0.06
    missing_rate: float = 0.0
    mnar_slope: float = 0.0
    low_quality_timepoints: tuple[int, ...] = ()
    low_quality_missingness_multiplier: float = 1.0
    corrupted_timepoints: tuple[int, ...] = ()
    corruption_multiplier: float = 2.2
    threshold_quantiles: tuple[float, float] = (0.45, 0.75)
    score_flip_rate: float = 0.0


DEFAULT_PROTEOMICS_PARAMETERS = ModalityParameters(
    signal_separation=1.8,
    baseline_sd=0.25,
    temporal_noise_sd=0.45,
    temporal_ar1=0.10,
    replicate_noise_sd=0.15,
    decoupling_rate=0.05,
    constitutive_gene_rate=0.03,
    missing_rate=0.45,
    mnar_slope=1.2,
    corrupted_timepoints=(1, 5),
    corruption_multiplier=2.2,
    threshold_quantiles=(0.50, 0.76),
)

DEFAULT_TRANSCRIPTOMICS_PARAMETERS = ModalityParameters(
    signal_separation=1.5,
    baseline_sd=0.35,
    temporal_noise_sd=0.80,
    temporal_ar1=0.0,
    replicate_noise_sd=0.25,
    decoupling_rate=0.10,
    constitutive_gene_rate=0.05,
    missing_rate=0.0,
    mnar_slope=0.0,
    corrupted_timepoints=(1, 5),
    corruption_multiplier=2.2,
    threshold_quantiles=(0.50, 0.78),
)


@dataclass(frozen=True)
class RegimeSupportEnvelope:
    """Cycle-free anchor flux and FVA-based reaction classification for one regime."""

    task_name: str
    anchor_flux: np.ndarray
    fva_min_abs: np.ndarray
    fva_max_abs: np.ndarray
    essential_active: np.ndarray
    alternate_active: np.ndarray
    inactive_or_decoy: np.ndarray


@dataclass(frozen=True)
class ModalityObservation:
    """Simulated gene-level expression and GPR-derived reaction scores for one modality."""

    name: str
    parameters: ModalityParameters
    latent_gene_state: np.ndarray
    replicate_expression: np.ndarray
    state_expression: np.ndarray
    thresholds: tuple[float, float]
    gene_scores: np.ndarray
    reaction_scores: Mapping[str, Mapping[str, float]]
    reaction_score_matrix: np.ndarray
    missing_fraction: float


@dataclass(frozen=True)
class FluxBenchmarkDataset:
    """Complete ground-truth flux trajectory and paired multi-omics observations."""

    graph: Any
    reactions: tuple[str, ...]
    pathways: tuple[str, ...]
    genes: tuple[str, ...]
    conditions: tuple[str, ...]
    temporal_pairs: tuple[tuple[str, str], ...]
    timepoints: np.ndarray
    eps: float
    transition_mode: str
    switch_timepoints: tuple[int, ...]
    regime_names_by_timepoint: tuple[str, ...]
    regime_envelopes: Mapping[str, RegimeSupportEnvelope]
    truth_flux: np.ndarray
    truth_activity: np.ndarray
    truth_union: np.ndarray
    truth_switches: Mapping[tuple[str, str], int]
    total_truth_switches: int
    truth_internal_cycle_flux: np.ndarray
    truth_gene_activity: np.ndarray
    active_complex_by_reaction: Mapping[str, tuple[str, ...]]
    timepoint_task_bounds: Mapping[str, Mapping[str, tuple[float, float]]]
    shared_inference_bounds: Mapping[str, Mapping[str, tuple[float, float]]]
    reaction_bounds: Mapping[str, Mapping[str, tuple[float, float]]]
    expose_task_bounds: bool
    proteomics: ModalityObservation
    transcriptomics: ModalityObservation


@dataclass(frozen=True)
class BenchmarkEvaluation:
    """Ground-truth evaluation metrics for a reconstructed temporal network."""

    active_tp: int
    active_tn: int
    active_fp: int
    active_fn: int
    active_precision: float
    active_recall: float
    active_f1: float
    active_jaccard: float
    per_timepoint_f1: tuple[float, ...]
    per_timepoint_precision: tuple[float, ...]
    per_timepoint_recall: tuple[float, ...]
    per_timepoint_tp: tuple[int, ...]
    per_timepoint_tn: tuple[int, ...]
    per_timepoint_fp: tuple[int, ...]
    per_timepoint_fn: tuple[int, ...]
    essential_recall: float
    non_decoy_precision: float
    unscored_active_recall: float
    predicted_switches: int
    true_switches: int
    switch_tp: int
    switch_tn: int
    switch_fp: int
    switch_fn: int
    switch_precision: float
    switch_recall: float
    switch_f1: float
    transition_macro_f1: float
    transition_class_f1: tuple[float, float, float, float]
    transition_confusion: tuple[tuple[int, int, int, int], ...]
    stationary_false_switches: int
    transition_switch_recall: float
    union_size: int
    true_union_size: int
    union_jaccard: float
    mean_internal_cycle_flux: float
    max_internal_cycle_flux: float


def _extract_network_matrices(
    graph: Any,
) -> tuple[csr_matrix, tuple[str, ...], np.ndarray, np.ndarray, np.ndarray, tuple[str, ...]]:
    """Return stoichiometric matrix S, reaction IDs, default lb/ub, exchange mask, and pathways."""
    n_edges = graph.num_edges
    reactions: list[str] = []
    pathways: list[str] = []
    default_lb = np.zeros(n_edges, dtype=float)
    default_ub = np.zeros(n_edges, dtype=float)
    is_exchange = np.zeros(n_edges, dtype=bool)

    for idx in range(n_edges):
        attr = graph.get_attr_edge(idx)
        rxn_id = attr.get("id", f"R_{idx}")
        reactions.append(str(rxn_id))
        pathways.append(str(attr.get("pathway", attr.get("subsystem", "Unassigned"))))
        default_lb[idx] = float(attr.get("default_lb", 0.0))
        default_ub[idx] = float(attr.get("default_ub", 1000.0))
        s_nodes, t_nodes = graph.get_edge(idx)
        if len(s_nodes) == 0 or len(t_nodes) == 0:
            is_exchange[idx] = True

    data, (rows, cols) = graph.get_vertex_incidence_matrix_as_lists(values=True)
    incidence = csr_matrix(
        (
            np.asarray(data, dtype=float),
            (np.asarray(rows, dtype=int), np.asarray(cols, dtype=int)),
        ),
        shape=(len(graph.V), n_edges),
    )
    return incidence, tuple(reactions), default_lb, default_ub, is_exchange, tuple(pathways)


def _solve_pfba(
    incidence: csr_matrix,
    lb: np.ndarray,
    ub: np.ndarray,
    *,
    objectives: Mapping[int, float] | None = None,
    opt_fraction: float = 1.0,
    previous_flux: np.ndarray | None = None,
    proximity_weight: float = 1e-3,
) -> tuple[np.ndarray, float, float | None]:
    """Solve parsimonious FBA (min ||v||_1) subject to Sv = 0 and bounds."""
    n_mets, n_rxns = incidence.shape
    bounds_list = list(zip(lb.tolist(), ub.tolist(), strict=True))

    lb_work = lb.copy()
    ub_work = ub.copy()
    objective_floor: float | None = None

    if objectives:
        c_fba = np.zeros(n_rxns, dtype=float)
        for rxn_idx, weight in objectives.items():
            c_fba[rxn_idx] = -float(weight)  # maximize weight * v
        res_fba = linprog(
            c_fba,
            A_eq=incidence,
            b_eq=np.zeros(n_mets, dtype=float),
            bounds=bounds_list,
            method="highs",
        )
        if not res_fba.success or res_fba.x is None:
            raise RuntimeError(f"Task FBA infeasible: {res_fba.message}")
        achieved_obj = float(-res_fba.fun)
        objective_floor = achieved_obj - (1.0 - opt_fraction) * abs(achieved_obj)

    # Minimize sum_r u_r (+ optional tiny proximity tie-breaker |v - v_prev|)
    # Variables: [v (n_rxns), u (n_rxns)] or [v, u, d] if previous_flux is provided
    I_n = eye(n_rxns, format="csc")
    zero_mets = csc_matrix((n_mets, n_rxns))

    if previous_flux is None:
        c = np.concatenate([np.zeros(n_rxns), np.ones(n_rxns)])
        A_eq = hstack([incidence, zero_mets], format="csc")
        b_eq = np.zeros(n_mets, dtype=float)
        # v - u <= 0 and -v - u <= 0
        A_ub = vstack(
            [
                hstack([I_n, -I_n]),
                hstack([-I_n, -I_n]),
            ],
            format="csc",
        )
        b_ub = np.zeros(2 * n_rxns, dtype=float)
        var_bounds = bounds_list + [(0.0, None)] * n_rxns
    else:
        c = np.concatenate(
            [
                np.zeros(n_rxns),
                np.ones(n_rxns),
                np.full(n_rxns, float(proximity_weight)),
            ]
        )
        A_eq = hstack([incidence, zero_mets, zero_mets], format="csc")
        b_eq = np.zeros(n_mets, dtype=float)
        zero_n = csc_matrix((n_rxns, n_rxns))
        # v - u <= 0, -v - u <= 0, v - d <= v_prev, -v - d <= -v_prev
        A_ub = vstack(
            [
                hstack([I_n, -I_n, zero_n]),
                hstack([-I_n, -I_n, zero_n]),
                hstack([I_n, zero_n, -I_n]),
                hstack([-I_n, zero_n, -I_n]),
            ],
            format="csc",
        )
        b_ub = np.concatenate(
            [
                np.zeros(2 * n_rxns, dtype=float),
                np.asarray(previous_flux, dtype=float),
                -np.asarray(previous_flux, dtype=float),
            ]
        )
        var_bounds = bounds_list + [(0.0, None)] * (2 * n_rxns)

    if objectives and objective_floor is not None:
        obj_row = np.zeros(c.size, dtype=float)
        for rxn_idx, weight in objectives.items():
            obj_row[rxn_idx] = -float(weight)
        A_ub = vstack([A_ub, csc_matrix(obj_row)], format="csc")
        b_ub = np.concatenate([b_ub, [-objective_floor]])

    res_pfba = linprog(
        c,
        A_ub=A_ub,
        b_ub=b_ub,
        A_eq=A_eq,
        b_eq=b_eq,
        bounds=var_bounds,
        method="highs",
    )
    if not res_pfba.success or res_pfba.x is None:
        raise RuntimeError(f"Parsimonious FBA infeasible: {res_pfba.message}")

    flux = np.asarray(res_pfba.x[:n_rxns], dtype=float)
    flux[np.abs(flux) < 1e-11] = 0.0
    l1_norm = float(np.abs(flux).sum())
    return flux, l1_norm, objective_floor


def compute_internal_circulation_flux(
    incidence: csr_matrix,
    lb: np.ndarray,
    ub: np.ndarray,
    is_exchange: np.ndarray,
    flux: np.ndarray,
    *,
    objectives: Mapping[int, float] | None = None,
    objective_floor: float | None = None,
) -> float:
    """Measure internal L1 excess ``||v_int||_1 - min ||v'_int||_1``.

    Fixes all exchange/boundary fluxes ``v[is_exchange]`` to their realized values
    in ``flux`` and solves for the minimum-L1 internal flux vector ``v'`` that
    satisfies ``S v' = 0`` and ``lb <= v' <= ub``. Any positive excess
    ``||v_int||_1 - ||v'_int||_1`` measures excess relative to parsimonious
    routing. It is not a thermodynamic loop test.
    """
    n_mets, n_rxns = incidence.shape
    internal_mask = ~is_exchange
    current_int_l1 = float(np.abs(flux[internal_mask]).sum())
    if current_int_l1 <= 1e-9:
        return 0.0

    lb_fixed = lb.copy()
    ub_fixed = ub.copy()
    lb_fixed[is_exchange] = flux[is_exchange]
    ub_fixed[is_exchange] = flux[is_exchange]
    lb_fixed = np.minimum(lb_fixed, flux)
    ub_fixed = np.maximum(ub_fixed, flux)

    I_n = eye(n_rxns, format="csc")
    zero_mets = csc_matrix((n_mets, n_rxns))
    weights = internal_mask.astype(float)
    c = np.concatenate([np.zeros(n_rxns), weights])
    A_eq = hstack([incidence, zero_mets], format="csc")
    b_eq = np.zeros(n_mets, dtype=float)
    A_ub = vstack([hstack([I_n, -I_n]), hstack([-I_n, -I_n])], format="csc")
    b_ub = np.zeros(2 * n_rxns, dtype=float)
    if objectives and objective_floor is not None:
        obj_row = np.zeros(2 * n_rxns, dtype=float)
        for rxn_idx, weight in objectives.items():
            obj_row[rxn_idx] = -float(weight)
        A_ub = vstack([A_ub, csc_matrix(obj_row)], format="csc")
        b_ub = np.concatenate([b_ub, [-objective_floor]])
    bounds_list = list(zip(lb_fixed.tolist(), ub_fixed.tolist(), strict=True)) + [
        (0.0, None)
    ] * n_rxns

    res = linprog(
        c,
        A_ub=A_ub,
        b_ub=b_ub,
        A_eq=A_eq,
        b_eq=b_eq,
        bounds=bounds_list,
        method="highs",
    )
    if not res.success or res.x is None:
        raise RuntimeError(f"Internal L1 minimization infeasible: {res.message}")
    min_int_l1 = float(np.dot(weights, np.abs(res.x[:n_rxns])))
    return max(0.0, current_int_l1 - min_int_l1)


def project_to_loopless_pfba(
    incidence: csr_matrix,
    lb: np.ndarray,
    ub: np.ndarray,
    is_exchange: np.ndarray,
    flux: np.ndarray,
    *,
    objectives: Mapping[int, float] | None = None,
    objective_floor: float | None = None,
) -> np.ndarray:
    """Minimize internal flux while preserving exchanges and an optional task objective."""
    n_mets, n_rxns = incidence.shape
    lb_fixed = lb.copy()
    ub_fixed = ub.copy()
    lb_fixed[is_exchange] = flux[is_exchange]
    ub_fixed[is_exchange] = flux[is_exchange]

    # Minimize ||v_int||_1 + 1e-2 * ||v - flux||_1, subject to any task objective.
    I_n = eye(n_rxns, format="csc")
    zero_mets = csc_matrix((n_mets, n_rxns))
    zero_n = csc_matrix((n_rxns, n_rxns))
    weights = (~is_exchange).astype(float)
    c = np.concatenate([np.zeros(n_rxns), weights, np.full(n_rxns, 1e-2)])
    A_eq = hstack([incidence, zero_mets, zero_mets], format="csc")
    b_eq = np.zeros(n_mets, dtype=float)
    A_ub = vstack(
        [
            hstack([I_n, -I_n, zero_n]),
            hstack([-I_n, -I_n, zero_n]),
            hstack([I_n, zero_n, -I_n]),
            hstack([-I_n, zero_n, -I_n]),
        ],
        format="csc",
    )
    b_ub = np.concatenate([np.zeros(2 * n_rxns), flux, -flux])
    if objectives and objective_floor is not None:
        obj_row = np.zeros(3 * n_rxns, dtype=float)
        for rxn_idx, weight in objectives.items():
            obj_row[rxn_idx] = -float(weight)
        A_ub = vstack([A_ub, csc_matrix(obj_row)], format="csc")
        b_ub = np.concatenate([b_ub, [-objective_floor]])
    bounds_list = list(zip(lb_fixed.tolist(), ub_fixed.tolist(), strict=True)) + [
        (0.0, None)
    ] * (2 * n_rxns)

    res = linprog(
        c,
        A_ub=A_ub,
        b_ub=b_ub,
        A_eq=A_eq,
        b_eq=b_eq,
        bounds=bounds_list,
        method="highs",
    )
    if not res.success or res.x is None:
        raise RuntimeError(f"Loopless projection infeasible: {res.message}")
    projected = np.asarray(res.x[:n_rxns], dtype=float)
    projected[np.abs(projected) < 1e-11] = 0.0
    return projected


def _compute_regime_envelope(
    incidence: csr_matrix,
    lb: np.ndarray,
    ub: np.ndarray,
    *,
    task_name: str,
    anchor_flux: np.ndarray,
    anchor_l1: float,
    eps: float,
    fva_l1_slack: float = 1.0 + 1e-5,
    compute_fva: bool = True,
    objectives: Mapping[int, float] | None = None,
    objective_floor: float | None = None,
) -> RegimeSupportEnvelope:
    """Classify reactions into essential, alternate-optimal, and inactive/decoy via FVA."""
    n_mets, n_rxns = incidence.shape
    anchor_active = np.abs(anchor_flux) >= (eps - _FEASIBILITY_TOL)

    if not compute_fva:
        return RegimeSupportEnvelope(
            task_name=task_name,
            anchor_flux=anchor_flux,
            fva_min_abs=np.abs(anchor_flux),
            fva_max_abs=np.abs(anchor_flux),
            essential_active=anchor_active.copy(),
            alternate_active=np.zeros(n_rxns, dtype=bool),
            inactive_or_decoy=~anchor_active,
        )

    # Constrain ||v||_1 <= fva_l1_slack * anchor_l1 so FVA explores alternate
    # parsimonious routes that achieve the regime task without futile loops.
    I_n = eye(n_rxns, format="csc")
    zero_mets = csc_matrix((n_mets, n_rxns))
    A_eq = hstack([incidence, zero_mets], format="csc")
    b_eq = np.zeros(n_mets, dtype=float)
    l1_row = csc_matrix(np.concatenate([np.zeros(n_rxns), np.ones(n_rxns)]))
    A_ub = vstack(
        [
            hstack([I_n, -I_n]),
            hstack([-I_n, -I_n]),
            l1_row,
        ],
        format="csc",
    )
    b_ub = np.concatenate([np.zeros(2 * n_rxns), [fva_l1_slack * max(anchor_l1, eps)]])
    if objectives and objective_floor is not None:
        obj_row = np.zeros(2 * n_rxns, dtype=float)
        for rxn_idx, weight in objectives.items():
            obj_row[rxn_idx] = -float(weight)
        A_ub = vstack([A_ub, csc_matrix(obj_row)], format="csc")
        b_ub = np.concatenate([b_ub, [-objective_floor]])
    var_bounds = list(zip(lb.tolist(), ub.tolist(), strict=True)) + [(0.0, None)] * n_rxns

    v_min = np.minimum(0.0, anchor_flux.copy())
    v_max = np.maximum(0.0, anchor_flux.copy())
    c = np.zeros(2 * n_rxns, dtype=float)

    for r in range(n_rxns):
        if lb[r] == ub[r]:
            v_min[r] = v_max[r] = lb[r]
            continue
        c.fill(0.0)
        c[r] = 1.0
        res_min = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=var_bounds, method="highs")
        if res_min.success and res_min.x is not None:
            v_min[r] = float(res_min.x[r])

        c.fill(0.0)
        c[r] = -1.0
        res_max = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=var_bounds, method="highs")
        if res_max.success and res_max.x is not None:
            v_max[r] = float(res_max.x[r])

    fva_max_abs = np.maximum(np.abs(v_min), np.abs(v_max))
    # A reaction can be zero if v_min <= 0 <= v_max; otherwise min |v_r| > 0
    crosses_zero = (v_min <= _FEASIBILITY_TOL) & (v_max >= -_FEASIBILITY_TOL)
    fva_min_abs = np.where(crosses_zero, 0.0, np.minimum(np.abs(v_min), np.abs(v_max)))

    essential_active = fva_min_abs >= (eps - _FEASIBILITY_TOL)
    can_be_active = (fva_max_abs >= (eps - _FEASIBILITY_TOL)) | anchor_active
    alternate_active = can_be_active & ~essential_active
    inactive_or_decoy = ~can_be_active

    return RegimeSupportEnvelope(
        task_name=task_name,
        anchor_flux=anchor_flux,
        fva_min_abs=fva_min_abs,
        fva_max_abs=fva_max_abs,
        essential_active=essential_active,
        alternate_active=alternate_active,
        inactive_or_decoy=inactive_or_decoy,
    )


def _parse_gpr_dnf(rule: str) -> list[tuple[str, ...]]:
    """Parse a boolean GPR string into Disjunctive Normal Form (list of AND-complexes)."""
    if not rule or not rule.strip():
        return []
    normalized = re.sub(r"\b(?:AND|OR)\b", lambda match: match.group().lower(), rule)
    try:
        tree = ast.parse(normalized, mode="eval")
    except SyntaxError:
        genes = sorted(get_genes_from_gpr(rule) or ())
        return [tuple(genes)] if genes else []

    def to_dnf(node: ast.AST) -> list[set[str]]:
        if isinstance(node, ast.Expression):
            return to_dnf(node.body)
        if isinstance(node, ast.Name):
            return [{node.id}]
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            return [{node.value}]
        if isinstance(node, ast.BoolOp):
            child_dnfs = [to_dnf(val) for val in node.values]
            if isinstance(node.op, ast.Or):
                merged: list[set[str]] = []
                for sub in child_dnfs:
                    merged.extend(sub)
                return merged
            if isinstance(node.op, ast.And):
                product = [set()]
                for sub in child_dnfs:
                    product = [a | b for a in product for b in sub]
                return product
        genes = sorted(get_genes_from_gpr(rule) or ())
        return [set(genes)] if genes else []

    clauses = to_dnf(tree)
    unique = sorted({tuple(sorted(clause)) for clause in clauses if clause})
    return unique


def _simulate_modality_observation(
    *,
    name: str,
    parameters: ModalityParameters,
    rng: np.random.Generator,
    graph: Any,
    genes: tuple[str, ...],
    reactions: tuple[str, ...],
    conditions: tuple[str, ...],
    truth_gene_activity: np.ndarray,
    replicates: int,
) -> ModalityObservation:
    """Simulate replicate/state gene expression and GPR reaction scores for one modality."""
    n_genes, n_times = truth_gene_activity.shape
    if not np.isfinite(parameters.low_quality_missingness_multiplier) or (
        parameters.low_quality_missingness_multiplier < 1.0
    ):
        raise ValueError("low_quality_missingness_multiplier must be finite and >= 1")
    if any(
        not isinstance(t, (int, np.integer)) or t < 0 or t >= n_times
        for t in parameters.low_quality_timepoints
    ):
        raise ValueError("low_quality_timepoints must be valid timepoint indices")
    latent_state = truth_gene_activity.astype(float).copy()

    # 1. Gene-level decoupling (allosteric / post-transcriptional decoupling of active enzymes)
    ever_active = np.any(truth_gene_activity, axis=1)
    active_candidates = np.flatnonzero(ever_active)
    if parameters.decoupling_rate > 0 and active_candidates.size > 0:
        n_decoupled = int(round(parameters.decoupling_rate * active_candidates.size))
        if n_decoupled > 0:
            chosen_dec = rng.choice(
                active_candidates, size=min(n_decoupled, active_candidates.size), replace=False
            )
            latent_state[chosen_dec, :] = 0.5

    # 2. Constitutively abundant inactive genes
    inactive_candidates = np.flatnonzero(~ever_active)
    if parameters.constitutive_gene_rate > 0 and inactive_candidates.size > 0:
        n_const = max(1, int(round(parameters.constitutive_gene_rate * inactive_candidates.size)))
        chosen = rng.choice(inactive_candidates, size=min(n_const, inactive_candidates.size), replace=False)
        latent_state[chosen, :] = 1.0

    # 3. Temporal AR(1) buffering (e.g., protein half-life smoothing)
    buffered_signal = np.zeros_like(latent_state)
    alpha = float(np.clip(parameters.temporal_ar1, 0.0, 0.95))
    buffered_signal[:, 0] = latent_state[:, 0]
    for t in range(1, n_times):
        buffered_signal[:, t] = alpha * buffered_signal[:, t - 1] + (1.0 - alpha) * latent_state[:, t]

    # 4. Static gene baseline + scaled activation signal + timepoint noise
    gene_baseline = rng.normal(0.0, parameters.baseline_sd, size=(n_genes, 1))
    centered_signal = (buffered_signal - 0.5) * parameters.signal_separation
    latent_abundance = gene_baseline + centered_signal

    time_noise_sd = np.full(n_times, parameters.temporal_noise_sd, dtype=float)
    for t_idx in parameters.corrupted_timepoints:
        if 0 <= t_idx < n_times:
            time_noise_sd[t_idx] *= parameters.corruption_multiplier

    state_noise = rng.normal(0.0, 1.0, size=(n_genes, n_times)) * time_noise_sd[None, :]
    noisy_state_expr = latent_abundance + state_noise

    replicate_expr = (
        noisy_state_expr[:, :, None]
        + rng.normal(0.0, parameters.replicate_noise_sd, size=(n_genes, n_times, replicates))
    )

    # 5. Missingness (MCAR + intensity-dependent MNAR for proteomics)
    missing_draw = None
    miss_prob = None
    if parameters.missing_rate > 0:
        # Detection depends on underlying abundance, not the noisy value that
        # will later be used to score a detected protein. Otherwise a negative
        # noise excursion is preferentially hidden as missing evidence.
        standardized = (latent_abundance - np.mean(latent_abundance)) / max(
            float(np.std(latent_abundance)), 1e-6
        )
        logits = -parameters.mnar_slope * standardized
        weights = 1.0 / (1.0 + np.exp(-logits))
        weights = weights / np.mean(weights)
        miss_prob = np.clip(parameters.missing_rate * weights, 0.0, 0.92)
        missing_draw = rng.random((n_genes, n_times))
        missing_mask = missing_draw < miss_prob
        for r_idx in range(replicates):
            replicate_expr[:, :, r_idx][missing_mask] = np.nan

    def state_mean() -> np.ndarray:
        state = np.full((n_genes, n_times), np.nan, dtype=float)
        valid_counts = np.isfinite(replicate_expr).sum(axis=2)
        valid_sums = np.nansum(replicate_expr, axis=2)
        np.divide(valid_sums, valid_counts, out=state, where=valid_counts > 0)
        return state

    # Calibrate global score thresholds before extra dropout, so paired
    # scenarios differ only at their low-quality timepoints.
    state_expr = state_mean()
    q_low, q_high = parameters.threshold_quantiles
    low_thresh, high_thresh = np.nanquantile(state_expr, [q_low, q_high])

    if parameters.low_quality_timepoints and miss_prob is not None:
        for t_idx in parameters.low_quality_timepoints:
            stressed_prob = np.clip(
                miss_prob[:, t_idx] * parameters.low_quality_missingness_multiplier,
                0.0,
                0.98,
            )
            extra_missing = missing_draw[:, t_idx] < stressed_prob
            replicate_expr[extra_missing, t_idx, :] = np.nan
        state_expr = state_mean()

    # 6. Assign scores with the shared thresholds.
    gene_scores = np.zeros((n_genes, n_times), dtype=float)
    finite_mask = np.isfinite(state_expr)
    gene_scores[finite_mask & (state_expr >= high_thresh)] = 1.0
    gene_scores[finite_mask & (state_expr <= low_thresh)] = -1.0

    # 7. Map gene expression to reaction scores using CORNETO's MultiSampleIMAT GPR rules
    cdict = {
        condition: {
            gene_id: {
                "value": float(state_expr[g_idx, t_idx]),
                "mapping": "none",
                "role": "expression",
            }
            for g_idx, gene_id in enumerate(genes)
            if np.isfinite(state_expr[g_idx, t_idx])
        }
        for t_idx, condition in enumerate(conditions)
    }
    scorer = MultiSampleIMAT(
        high_expression_threshold=float(high_thresh),
        low_expression_threshold=float(low_thresh),
        gpr_field="GPR",
    )
    _, processed_data = scorer.preprocess(graph, Data.from_cdict(cdict))

    rxn_index = {rxn_id: i for i, rxn_id in enumerate(reactions)}
    score_matrix = np.zeros((len(reactions), n_times), dtype=float)
    reaction_scores: dict[str, dict[str, float]] = {}

    for t_idx, condition in enumerate(conditions):
        sample = processed_data.samples[condition]
        cond_scores: dict[str, float] = {}
        for feat in sample.features:
            if feat.mapping == "edge" and feat.value is not None:
                val = float(feat.value)
                if val != 0.0 and feat.id in rxn_index:
                    cond_scores[feat.id] = val
                    score_matrix[rxn_index[feat.id], t_idx] = val
        reaction_scores[condition] = cond_scores

    # Optional post-GPR score flip stress test
    if parameters.score_flip_rate > 0:
        flip_mask = rng.random(score_matrix.shape) < parameters.score_flip_rate
        for r_idx, t_idx in np.argwhere(flip_mask):
            current = score_matrix[r_idx, t_idx]
            choices = [v for v in (-1.0, 0.0, 1.0) if v != current]
            new_val = float(rng.choice(choices))
            score_matrix[r_idx, t_idx] = new_val
            cond = conditions[t_idx]
            rxn_id = reactions[r_idx]
            if new_val == 0.0:
                reaction_scores[cond].pop(rxn_id, None)
            else:
                reaction_scores[cond][rxn_id] = new_val

    missing_fraction = float(np.isnan(state_expr).mean())
    return ModalityObservation(
        name=name,
        parameters=parameters,
        latent_gene_state=latent_state,
        replicate_expression=replicate_expr,
        state_expression=state_expr,
        thresholds=(float(low_thresh), float(high_thresh)),
        gene_scores=gene_scores,
        reaction_scores=reaction_scores,
        reaction_score_matrix=score_matrix,
        missing_fraction=missing_fraction,
    )


def simulate_flux_benchmark(
    graph: Any,
    *,
    tasks: Sequence[MetabolicTask] | None = None,
    timepoints: int | Sequence[float] = 8,
    transition_mode: str = "step",
    switch_timepoints: Sequence[int] = (4,),
    eps: float = 1e-2,
    enforce_loopless_truth: bool = True,
    expose_task_bounds: bool = False,
    compute_fva: bool = True,
    fva_l1_slack: float = 1.0 + 1e-5,
    isozyme_mode: str = "single_dominant",
    replicates: int = 3,
    proteomics_parameters: ModalityParameters | None = None,
    transcriptomics_parameters: ModalityParameters | None = None,
    seed: int = 17,
) -> FluxBenchmarkDataset:
    """Generate a ground-truth steady-state flux trajectory and paired multi-omics scores.

    Parameters
    ----------
    graph:
        Any :class:`corneto.Graph` metabolic network with reaction IDs, bounds,
        and ``GPR`` edge attributes.
    tasks:
        Sequence of :class:`MetabolicTask` specifications defining the metabolic
        regimes. When ``None`` and ``graph`` was created by
        :func:`build_core_metabolic_network`, uses the default hepatic-like
        Glycolytic/Lipogenic -> Gluconeogenic/Beta-Oxidation tasks.
    timepoints:
        Number of timepoints (int >= 3) or explicit increasing time values.
    transition_mode:
        ``"step"`` (piecewise constant regimes switching at ``switch_timepoints``),
        ``"one_point"`` (1-step boundary interpolation projected to a parsimonious
        pFBA state), or ``"pulse"`` (Task 0 baseline with Task 1 active between
        ``switch_timepoints[0]`` and ``switch_timepoints[1]``).
    switch_timepoints:
        Timepoint indices where regime transitions occur.
    eps:
        Threshold defining ground-truth active reactions ``|v*(t)| >= eps``.
    enforce_loopless_truth:
        Project every timepoint's ground-truth flux onto a minimum-internal-L1
        state, preserving task objectives, and verify near-zero internal L1 excess.
    expose_task_bounds:
        If ``False`` (default), ``dataset.reaction_bounds`` provides the common
        permissive bound envelope across all tasks so the solver must infer
        regime transitions from the omics data rather than boundary constraints.
        If ``True``, exposes the timepoint-specific task bounds.
    """
    if transition_mode not in {"step", "one_point", "pulse"}:
        raise ValueError("transition_mode must be 'step', 'one_point', or 'pulse'")
    if isozyme_mode not in {"single_dominant", "switching_isozyme", "all"}:
        raise ValueError("isozyme_mode must be 'single_dominant', 'switching_isozyme', or 'all'")
    if eps <= 0 or not np.isfinite(eps):
        raise ValueError("eps must be a finite positive number")

    if isinstance(timepoints, (int, np.integer)):
        if int(timepoints) < 3:
            raise ValueError("timepoints must be at least 3")
        times = np.linspace(0.0, 1.0, int(timepoints))
    else:
        times = np.asarray(timepoints, dtype=float)
        if times.ndim != 1 or times.size < 3 or np.any(np.diff(times) <= 0):
            raise ValueError("timepoints must be a 1D strictly increasing sequence of >= 3 values")
    n_times = int(times.size)
    conditions = tuple(f"t{i}" for i in range(n_times))
    temporal_pairs = tuple(zip(conditions[:-1], conditions[1:], strict=True))

    if tasks is None:
        default_tasks = graph.get_graph_attributes().get("default_benchmark_tasks")
        if not default_tasks:
            raise ValueError("tasks must be provided when graph has no default_benchmark_tasks")
        tasks = default_tasks
    tasks = tuple(tasks)
    if len(tasks) < 2:
        raise ValueError("At least two MetabolicTask regimes are required to benchmark transitions")

    switch_pts = tuple(int(x) for x in switch_timepoints)
    if any(pt <= 0 or pt >= n_times for pt in switch_pts):
        raise ValueError(f"switch_timepoints must lie strictly inside (0, {n_times})")

    rng = np.random.default_rng(seed)
    incidence, reactions, default_lb, default_ub, is_exchange, pathways = _extract_network_matrices(
        graph
    )
    rxn_to_idx = {rxn_id: idx for idx, rxn_id in enumerate(reactions)}
    n_rxns = len(reactions)

    # 1. Solve anchor pFBA and FVA envelope for each task
    task_lbs: dict[str, np.ndarray] = {}
    task_ubs: dict[str, np.ndarray] = {}
    task_objectives: dict[str, Mapping[int, float] | None] = {}
    task_objective_floors: dict[str, float | None] = {}
    regime_envelopes: dict[str, RegimeSupportEnvelope] = {}
    prev_anchor: np.ndarray | None = None

    for task in tasks:
        lb_k = default_lb.copy()
        ub_k = default_ub.copy()
        for rxn_id, (low, high) in task.bounds.items():
            if rxn_id not in rxn_to_idx:
                raise ValueError(f"Task {task.name!r} references unknown reaction {rxn_id!r}")
            r_idx = rxn_to_idx[rxn_id]
            lb_k[r_idx] = float(low)
            ub_k[r_idx] = float(high)
        obj_idx = (
            {rxn_to_idx[k]: float(v) for k, v in task.objectives.items()}
            if task.objectives
            else None
        )
        anchor_flux, anchor_l1, objective_floor = _solve_pfba(
            incidence,
            lb_k,
            ub_k,
            objectives=obj_idx,
            opt_fraction=task.opt_fraction,
            previous_flux=prev_anchor,
        )
        if enforce_loopless_truth:
            anchor_flux = project_to_loopless_pfba(
                incidence,
                lb_k,
                ub_k,
                is_exchange,
                anchor_flux,
                objectives=obj_idx,
                objective_floor=objective_floor,
            )
            anchor_l1 = float(np.abs(anchor_flux).sum())

        envelope = _compute_regime_envelope(
            incidence,
            lb_k,
            ub_k,
            task_name=task.name,
            anchor_flux=anchor_flux,
            anchor_l1=anchor_l1,
            eps=eps,
            fva_l1_slack=fva_l1_slack,
            compute_fva=compute_fva,
            objectives=obj_idx,
            objective_floor=objective_floor,
        )
        task_lbs[task.name] = lb_k
        task_ubs[task.name] = ub_k
        task_objectives[task.name] = obj_idx
        task_objective_floors[task.name] = objective_floor
        regime_envelopes[task.name] = envelope
        prev_anchor = anchor_flux

    # 2. Assign regime schedule across timepoints t = 0..T-1
    regime_indices = np.zeros(n_times, dtype=int)
    if transition_mode in {"step", "one_point"}:
        for step_idx, sw in enumerate(switch_pts):
            regime_indices[sw:] = min(step_idx + 1, len(tasks) - 1)
    elif transition_mode == "pulse":
        start_t = switch_pts[0]
        end_t = switch_pts[1] if len(switch_pts) > 1 else min(start_t + 2, n_times - 1)
        regime_indices[start_t:end_t] = 1

    truth_flux = np.zeros((n_rxns, n_times), dtype=float)
    timepoint_lbs = np.zeros((n_rxns, n_times), dtype=float)
    timepoint_ubs = np.zeros((n_rxns, n_times), dtype=float)
    regime_names_by_timepoint: list[str] = []
    timepoint_objectives: list[Mapping[int, float] | None] = []
    timepoint_objective_floors: list[float | None] = []

    for t in range(n_times):
        task_k = tasks[int(regime_indices[t])]
        regime_names_by_timepoint.append(task_k.name)
        lb_t = task_lbs[task_k.name].copy()
        ub_t = task_ubs[task_k.name].copy()
        flux_t = regime_envelopes[task_k.name].anchor_flux.copy()
        objective_t = task_objectives[task_k.name]
        objective_floor_t = task_objective_floors[task_k.name]

        if transition_mode == "one_point" and t in switch_pts and t > 0:
            prev_task = tasks[int(regime_indices[t - 1])]
            lb_t = 0.5 * task_lbs[prev_task.name] + 0.5 * task_lbs[task_k.name]
            ub_t = 0.5 * task_ubs[prev_task.name] + 0.5 * task_ubs[task_k.name]
            mix_flux = (
                0.5 * regime_envelopes[prev_task.name].anchor_flux
                + 0.5 * regime_envelopes[task_k.name].anchor_flux
            )
            if enforce_loopless_truth:
                flux_t = project_to_loopless_pfba(incidence, lb_t, ub_t, is_exchange, mix_flux)
            else:
                flux_t = mix_flux
            objective_t = None
            objective_floor_t = None

        truth_flux[:, t] = flux_t
        timepoint_lbs[:, t] = lb_t
        timepoint_ubs[:, t] = ub_t
        timepoint_objectives.append(objective_t)
        timepoint_objective_floors.append(objective_floor_t)

    # 3. Verify mass balance, bounds, and internal circulation on truth_flux
    residual = incidence @ truth_flux
    if np.max(np.abs(residual), initial=0.0) > _FEASIBILITY_TOL:
        raise RuntimeError("Generated ground-truth flux violates Sv = 0")

    shared_lb = np.min(timepoint_lbs, axis=1)
    shared_ub = np.max(timepoint_ubs, axis=1)
    if np.any(truth_flux < timepoint_lbs - _FEASIBILITY_TOL) or np.any(
        truth_flux > timepoint_ubs + _FEASIBILITY_TOL
    ):
        raise RuntimeError("Generated ground-truth flux violates timepoint task bounds")

    truth_cycle_flux = np.array(
        [
            compute_internal_circulation_flux(
                incidence,
                timepoint_lbs[:, t],
                timepoint_ubs[:, t],
                is_exchange,
                truth_flux[:, t],
                objectives=timepoint_objectives[t],
                objective_floor=timepoint_objective_floors[t],
            )
            for t in range(n_times)
        ],
        dtype=float,
    )
    if enforce_loopless_truth and np.max(truth_cycle_flux, initial=0.0) > 1e-5:
        raise RuntimeError(
            f"Ground-truth flux has excess internal L1 flux: max={np.max(truth_cycle_flux):.6g}"
        )

    truth_activity = np.abs(truth_flux) >= (eps - _FEASIBILITY_TOL)
    truth_union = np.any(truth_activity, axis=1)
    truth_switches = {
        (conditions[t], conditions[t + 1]): int(
            np.count_nonzero(truth_activity[:, t] != truth_activity[:, t + 1])
        )
        for t in range(n_times - 1)
    }
    total_truth_switches = int(sum(truth_switches.values()))
    if total_truth_switches == 0:
        raise RuntimeError(
            "Ground-truth trajectory has zero active-network switches at eps; "
            "ensure tasks have distinct active supports above eps."
        )

    # Build bounds dictionaries for inference
    timepoint_task_bounds: dict[str, dict[str, tuple[float, float]]] = {}
    shared_inference_bounds: dict[str, dict[str, tuple[float, float]]] = {}
    for t, cond in enumerate(conditions):
        t_dict: dict[str, tuple[float, float]] = {}
        s_dict: dict[str, tuple[float, float]] = {}
        for r_idx, rxn_id in enumerate(reactions):
            if (
                timepoint_lbs[r_idx, t] != default_lb[r_idx]
                or timepoint_ubs[r_idx, t] != default_ub[r_idx]
            ):
                t_dict[rxn_id] = (float(timepoint_lbs[r_idx, t]), float(timepoint_ubs[r_idx, t]))
            if shared_lb[r_idx] != default_lb[r_idx] or shared_ub[r_idx] != default_ub[r_idx]:
                s_dict[rxn_id] = (float(shared_lb[r_idx]), float(shared_ub[r_idx]))
        timepoint_task_bounds[cond] = t_dict
        shared_inference_bounds[cond] = s_dict

    # 4. Map active reactions to true functional gene activity via GPR DNF clauses
    all_genes: set[str] = set()
    rxn_dnf: dict[int, list[tuple[str, ...]]] = {}
    for r_idx in range(n_rxns):
        gpr = str(graph.get_attr_edge(r_idx).get("GPR", "") or "")
        clauses = _parse_gpr_dnf(gpr)
        rxn_dnf[r_idx] = clauses
        for clause in clauses:
            all_genes.update(clause)
    genes = tuple(sorted(all_genes))
    gene_to_idx = {g: i for i, g in enumerate(genes)}

    truth_gene_activity = np.zeros((len(genes), n_times), dtype=bool)
    active_complex_by_reaction: dict[str, tuple[str, ...]] = {}

    for r_idx, rxn_id in enumerate(reactions):
        clauses = rxn_dnf[r_idx]
        if not clauses:
            continue
        if isozyme_mode == "single_dominant":
            chosen_idx = int(rng.integers(0, len(clauses)))
            chosen_clause = clauses[chosen_idx]
            active_complex_by_reaction[rxn_id] = chosen_clause
            for t in range(n_times):
                if truth_activity[r_idx, t]:
                    for g in chosen_clause:
                        truth_gene_activity[gene_to_idx[g], t] = True
        elif isozyme_mode == "all":
            merged_genes = tuple(sorted({g for clause in clauses for g in clause}))
            active_complex_by_reaction[rxn_id] = merged_genes
            for t in range(n_times):
                if truth_activity[r_idx, t]:
                    for g in merged_genes:
                        truth_gene_activity[gene_to_idx[g], t] = True
        else:  # switching_isozyme across regime boundary
            active_complex_by_reaction[rxn_id] = clauses[0]
            for t in range(n_times):
                if truth_activity[r_idx, t]:
                    clause = clauses[int(regime_indices[t]) % len(clauses)]
                    for g in clause:
                        truth_gene_activity[gene_to_idx[g], t] = True

    # 5. Simulate paired proteomics and transcriptomics modalities
    prot_params = proteomics_parameters or DEFAULT_PROTEOMICS_PARAMETERS
    rna_params = transcriptomics_parameters or DEFAULT_TRANSCRIPTOMICS_PARAMETERS

    proteomics = _simulate_modality_observation(
        name="proteomics",
        parameters=prot_params,
        rng=rng,
        graph=graph,
        genes=genes,
        reactions=reactions,
        conditions=conditions,
        truth_gene_activity=truth_gene_activity,
        replicates=replicates,
    )
    transcriptomics = _simulate_modality_observation(
        name="transcriptomics",
        parameters=rna_params,
        rng=rng,
        graph=graph,
        genes=genes,
        reactions=reactions,
        conditions=conditions,
        truth_gene_activity=truth_gene_activity,
        replicates=replicates,
    )

    return FluxBenchmarkDataset(
        graph=graph,
        reactions=reactions,
        pathways=pathways,
        genes=genes,
        conditions=conditions,
        temporal_pairs=temporal_pairs,
        timepoints=times,
        eps=float(eps),
        transition_mode=transition_mode,
        switch_timepoints=switch_pts,
        regime_names_by_timepoint=tuple(regime_names_by_timepoint),
        regime_envelopes=regime_envelopes,
        truth_flux=truth_flux,
        truth_activity=truth_activity,
        truth_union=truth_union,
        truth_switches=truth_switches,
        total_truth_switches=total_truth_switches,
        truth_internal_cycle_flux=truth_cycle_flux,
        truth_gene_activity=truth_gene_activity,
        active_complex_by_reaction=active_complex_by_reaction,
        timepoint_task_bounds=timepoint_task_bounds,
        shared_inference_bounds=shared_inference_bounds,
        reaction_bounds=timepoint_task_bounds if expose_task_bounds else shared_inference_bounds,
        expose_task_bounds=bool(expose_task_bounds),
        proteomics=proteomics,
        transcriptomics=transcriptomics,
    )


def evaluate_benchmark_prediction(
    dataset: FluxBenchmarkDataset,
    *,
    flux: np.ndarray,
    activity: np.ndarray | None = None,
    modality: str | ModalityObservation | None = None,
) -> BenchmarkEvaluation:
    """Evaluate predicted fluxes and thresholded activity against ground truth."""
    pred_flux = np.asarray(flux, dtype=float)
    if pred_flux.shape != dataset.truth_flux.shape:
        raise ValueError(
            f"flux shape {pred_flux.shape} does not match truth {dataset.truth_flux.shape}"
        )
    if activity is None:
        pred_act = np.abs(pred_flux) >= (dataset.eps - _FEASIBILITY_TOL)
    else:
        pred_act = np.asarray(activity, dtype=bool)
        if pred_act.shape != dataset.truth_activity.shape:
            raise ValueError("activity shape does not match truth_activity")

    true_act = dataset.truth_activity
    tp = int(np.count_nonzero(pred_act & true_act))
    tn = int(np.count_nonzero(~pred_act & ~true_act))
    fp = int(np.count_nonzero(pred_act & ~true_act))
    fn = int(np.count_nonzero(~pred_act & true_act))
    precision = tp / max(1, tp + fp)
    recall = tp / max(1, tp + fn)
    f1 = 2.0 * precision * recall / max(1e-12, precision + recall)
    jaccard = tp / max(1, tp + fp + fn)

    n_times = len(dataset.conditions)
    per_t_f1: list[float] = []
    per_t_precision: list[float] = []
    per_t_recall: list[float] = []
    per_t_tp: list[int] = []
    per_t_tn: list[int] = []
    per_t_fp: list[int] = []
    per_t_fn: list[int] = []
    ess_hits = 0
    ess_total = 0
    non_decoy_hits = 0
    pred_total = int(np.count_nonzero(pred_act))

    for t in range(n_times):
        tp_t = int(np.count_nonzero(pred_act[:, t] & true_act[:, t]))
        tn_t = int(np.count_nonzero(~pred_act[:, t] & ~true_act[:, t]))
        fp_t = int(np.count_nonzero(pred_act[:, t] & ~true_act[:, t]))
        fn_t = int(np.count_nonzero(~pred_act[:, t] & true_act[:, t]))
        per_t_tp.append(tp_t)
        per_t_tn.append(tn_t)
        per_t_fp.append(fp_t)
        per_t_fn.append(fn_t)
        p_t = tp_t / max(1, tp_t + fp_t)
        r_t = tp_t / max(1, tp_t + fn_t)
        per_t_precision.append(float(p_t))
        per_t_recall.append(float(r_t))
        per_t_f1.append(float(2.0 * p_t * r_t / max(1e-12, p_t + r_t)))

        env = dataset.regime_envelopes[dataset.regime_names_by_timepoint[t]]
        ess_hits += int(np.count_nonzero(pred_act[:, t] & env.essential_active))
        ess_total += int(np.count_nonzero(env.essential_active))
        valid_support = env.essential_active | env.alternate_active | true_act[:, t]
        non_decoy_hits += int(np.count_nonzero(pred_act[:, t] & valid_support))

    essential_recall = ess_hits / max(1, ess_total)
    non_decoy_precision = non_decoy_hits / max(1, pred_total)

    # Recall of active reactions without a positive or negative reaction score.
    # Zero can result from missing data, mid-range expression, or absent GPRs;
    # it must not be interpreted as missing-protein recall alone.
    if isinstance(modality, str):
        mod_obs = (
            dataset.proteomics if modality == "proteomics" else dataset.transcriptomics
        )
    else:
        mod_obs = modality
    if mod_obs is not None:
        unscored_true = true_act & (mod_obs.reaction_score_matrix == 0.0)
        unscored_hits = int(np.count_nonzero(pred_act & unscored_true))
        unscored_total = int(np.count_nonzero(unscored_true))
        unscored_active_recall = unscored_hits / max(1, unscored_total)
    else:
        unscored_active_recall = recall

    # Switch evaluation across adjacent timepoints
    pred_sw_mat = pred_act[:, 1:] != pred_act[:, :-1]
    true_sw_mat = true_act[:, 1:] != true_act[:, :-1]
    pred_switches = int(np.count_nonzero(pred_sw_mat))
    true_switches = int(np.count_nonzero(true_sw_mat))
    sw_tp = int(np.count_nonzero(pred_sw_mat & true_sw_mat))
    sw_tn = int(np.count_nonzero(~pred_sw_mat & ~true_sw_mat))
    sw_fp = int(np.count_nonzero(pred_sw_mat & ~true_sw_mat))
    sw_fn = int(np.count_nonzero(~pred_sw_mat & true_sw_mat))
    sw_prec = sw_tp / max(1, sw_tp + sw_fp)
    sw_rec = sw_tp / max(1, sw_tp + sw_fn)
    sw_f1 = 2.0 * sw_prec * sw_rec / max(1e-12, sw_prec + sw_rec)

    # One four-class decision per reaction and adjacent pair: 00 stays off,
    # 01 turns on, 10 turns off, and 11 stays on. This jointly evaluates
    # activity and transition direction without mixing two separate scores.
    true_transition = 2 * true_act[:, :-1].astype(np.int8) + true_act[:, 1:].astype(np.int8)
    pred_transition = 2 * pred_act[:, :-1].astype(np.int8) + pred_act[:, 1:].astype(np.int8)
    transition_confusion = np.zeros((4, 4), dtype=int)
    np.add.at(transition_confusion, (true_transition.ravel(), pred_transition.ravel()), 1)
    transition_denominator = transition_confusion.sum(axis=1) + transition_confusion.sum(axis=0)
    transition_class_f1 = np.zeros(4, dtype=float)
    np.divide(
        2 * np.diag(transition_confusion),
        transition_denominator,
        out=transition_class_f1,
        where=transition_denominator > 0,
    )
    # A class absent from both truth and prediction has undefined F1; omit it.
    transition_macro_f1 = float(np.mean(transition_class_f1[transition_denominator > 0]))

    stationary_pairs = [
        t
        for t in range(n_times - 1)
        if dataset.regime_names_by_timepoint[t] == dataset.regime_names_by_timepoint[t + 1]
    ]
    transition_pairs = [
        t
        for t in range(n_times - 1)
        if dataset.regime_names_by_timepoint[t] != dataset.regime_names_by_timepoint[t + 1]
    ]
    stationary_false_switches = (
        int(np.count_nonzero(pred_sw_mat[:, stationary_pairs])) if stationary_pairs else 0
    )
    if transition_pairs:
        trans_tp = int(
            np.count_nonzero(
                pred_sw_mat[:, transition_pairs] & true_sw_mat[:, transition_pairs]
            )
        )
        trans_total = int(np.count_nonzero(true_sw_mat[:, transition_pairs]))
        transition_switch_recall = trans_tp / max(1, trans_total)
    else:
        transition_switch_recall = 1.0

    pred_union = np.any(pred_act, axis=1)
    true_union = dataset.truth_union
    union_tp = int(np.count_nonzero(pred_union & true_union))
    union_denom = int(np.count_nonzero(pred_union | true_union))
    union_jaccard = union_tp / max(1, union_denom)

    # Internal circulation excess on predicted fluxes
    incidence, _, default_lb, default_ub, is_exchange, _ = _extract_network_matrices(
        dataset.graph
    )
    rxn_to_idx = {r: i for i, r in enumerate(dataset.reactions)}
    cycle_fluxes: list[float] = []
    for t, cond in enumerate(dataset.conditions):
        lb_t = default_lb.copy()
        ub_t = default_ub.copy()
        for r_id, (low, high) in dataset.reaction_bounds.get(cond, {}).items():
            lb_t[rxn_to_idx[r_id]] = float(low)
            ub_t[rxn_to_idx[r_id]] = float(high)
        cycle_fluxes.append(
            compute_internal_circulation_flux(
                incidence, lb_t, ub_t, is_exchange, pred_flux[:, t]
            )
        )

    return BenchmarkEvaluation(
        active_tp=tp,
        active_tn=tn,
        active_fp=fp,
        active_fn=fn,
        active_precision=float(precision),
        active_recall=float(recall),
        active_f1=float(f1),
        active_jaccard=float(jaccard),
        per_timepoint_f1=tuple(per_t_f1),
        per_timepoint_precision=tuple(per_t_precision),
        per_timepoint_recall=tuple(per_t_recall),
        per_timepoint_tp=tuple(per_t_tp),
        per_timepoint_tn=tuple(per_t_tn),
        per_timepoint_fp=tuple(per_t_fp),
        per_timepoint_fn=tuple(per_t_fn),
        essential_recall=float(essential_recall),
        non_decoy_precision=float(non_decoy_precision),
        unscored_active_recall=float(unscored_active_recall),
        predicted_switches=pred_switches,
        true_switches=true_switches,
        switch_tp=sw_tp,
        switch_tn=sw_tn,
        switch_fp=sw_fp,
        switch_fn=sw_fn,
        switch_precision=float(sw_prec),
        switch_recall=float(sw_rec),
        switch_f1=float(sw_f1),
        transition_macro_f1=transition_macro_f1,
        transition_class_f1=tuple(float(value) for value in transition_class_f1),
        transition_confusion=tuple(
            tuple(int(value) for value in row) for row in transition_confusion
        ),
        stationary_false_switches=stationary_false_switches,
        transition_switch_recall=float(transition_switch_recall),
        union_size=int(np.count_nonzero(pred_union)),
        true_union_size=int(np.count_nonzero(true_union)),
        union_jaccard=float(union_jaccard),
        mean_internal_cycle_flux=float(np.mean(cycle_fluxes)),
        max_internal_cycle_flux=float(np.max(cycle_fluxes, initial=0.0)),
    )


def build_core_metabolic_network() -> tuple[cn.Graph, tuple[MetabolicTask, MetabolicTask]]:
    """Build a mass-balanced ~55-reaction hepatic core network with two metabolic regimes.

    The network models a switch between:
    1. ``"fed_glycolytic_lipogenic"``: Glucose import -> Glycolysis + Pentose
       Phosphate Pathway (NADPH) + De novo lipogenesis + ATP maintenance.
    2. ``"fasted_gluconeogenic_ketogenic"``: Fatty-acid + Lactate/Alanine import ->
       Mitochondrial beta-oxidation + Ketogenesis + Gluconeogenesis -> Glucose
       & Ketone export + ATP maintenance.

    Includes:
    - Bidirectional substrate-cycle pairs (GK/G6PC, PFK1/FBP1, PK/PC+PEPCK)
    - Parallel NADH/redox shuttles and isozyme families (``OR`` GPRs)
    - Multi-subunit enzyme complexes (``AND`` GPRs)
    - Decoy branches (glycogen cycling, hexosamine, pentose salvage, futile loops)
    """
    specs = [
        # --- Exchange & Transport ---
        ("EX_GLC_IN", {}, {"GLC_c": 1.0}, 0.0, 10.0, "Transport", "SLC2A2 or SLC2A1"),
        ("EX_GLC_OUT", {"GLC_c": 1.0}, {}, 0.0, 10.0, "Transport", "SLC2A2"),
        ("EX_LAC_IN", {}, {"LAC_c": 1.0}, 0.0, 10.0, "Transport", "SLC16A1 or SLC16A7"),
        ("EX_LAC_OUT", {"LAC_c": 1.0}, {}, 0.0, 10.0, "Transport", "SLC16A1"),
        ("EX_ALA_IN", {}, {"ALA_c": 1.0}, 0.0, 10.0, "Transport", "SLC38A2 or SLC38A4"),
        ("EX_FA_IN", {}, {"FA_c": 1.0}, 0.0, 10.0, "Transport", "CD36 or SLC27A2"),
        ("EX_LIPID_OUT", {"PALM_c": 1.0}, {}, 0.0, 10.0, "Lipogenesis", "MTTP and APOB"),
        ("EX_KET_OUT", {"BHB_m": 1.0}, {}, 0.0, 10.0, "Ketogenesis", "SLC16A6 or SLC16A1"),
        ("EX_GLYCOGEN_OUT", {"GLYG_c": 1.0}, {}, 0.0, 10.0, "Glycogen", ""),
        ("EX_HEXOSAMINE_OUT", {"UDPGNAc_c": 1.0}, {}, 0.0, 10.0, "Decoy_Hexosamine", ""),
        ("ATP_MAINT", {"ATP_c": 1.0}, {"ADP_c": 1.0}, 1.0, 10.0, "Maintenance", ""),

        # --- Upper Glycolysis & Gluconeogenesis (with substrate cycle pairs) ---
        ("GK", {"GLC_c": 1.0, "ATP_c": 1.0}, {"G6P_c": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Glycolysis", "GCK or HK1"),
        ("G6PC", {"G6P_c": 1.0}, {"GLC_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "G6PC and SLC37A4"),
        ("GPI_F", {"G6P_c": 1.0}, {"F6P_c": 1.0}, 0.0, 10.0, "Glycolysis", "GPI"),
        ("GPI_R", {"F6P_c": 1.0}, {"G6P_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "GPI"),
        ("PFK1", {"F6P_c": 1.0, "ATP_c": 1.0}, {"FBP_c": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Glycolysis", "PFKL or PFKM"),
        ("FBP1", {"FBP_c": 1.0}, {"F6P_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "FBP1"),
        ("ALDO_F", {"FBP_c": 1.0}, {"GAP_c": 2.0}, 0.0, 10.0, "Glycolysis", "ALDOB or ALDOA"),
        ("ALDO_R", {"GAP_c": 2.0}, {"FBP_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "ALDOB"),

        # --- Lower Glycolysis & Gluconeogenesis ---
        ("GAPDH_PGK_F", {"GAP_c": 1.0, "NAD_c": 1.0, "ADP_c": 1.0}, {"PGA3_c": 1.0, "NADH_c": 1.0, "ATP_c": 1.0}, 0.0, 10.0, "Glycolysis", "GAPDH and PGK1"),
        ("GAPDH_PGK_R", {"PGA3_c": 1.0, "NADH_c": 1.0, "ATP_c": 1.0}, {"GAP_c": 1.0, "NAD_c": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "GAPDH and PGK1"),
        ("ENO_F", {"PGA3_c": 1.0}, {"PEP_c": 1.0}, 0.0, 10.0, "Glycolysis", "ENO1 or ENO3"),
        ("ENO_R", {"PEP_c": 1.0}, {"PGA3_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "ENO1"),
        ("PK", {"PEP_c": 1.0, "ADP_c": 1.0}, {"PYR_c": 1.0, "ATP_c": 1.0}, 0.0, 10.0, "Glycolysis", "PKLR or PKM"),

        # --- Lactate / Alanine / Pyruvate ---
        ("LDH_F", {"PYR_c": 1.0, "NADH_c": 1.0}, {"LAC_c": 1.0, "NAD_c": 1.0}, 0.0, 10.0, "Glycolysis", "LDHA"),
        ("LDH_R", {"LAC_c": 1.0, "NAD_c": 1.0}, {"PYR_c": 1.0, "NADH_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "LDHB or LDHA"),
        ("ALT_IN", {"ALA_c": 1.0}, {"PYR_c": 1.0, "NADH_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "GPT or GPT2"),
        ("MPC", {"PYR_c": 1.0}, {"PYR_m": 1.0}, -10.0, 10.0, "Transport", "MPC1 and MPC2"),

        # --- Anaplerosis & PEPCK (Gluconeogenesis commit) ---
        ("PC", {"PYR_m": 1.0, "ATP_m": 1.0}, {"OAA_m": 1.0, "ADP_m": 1.0}, 0.0, 10.0, "Gluconeogenesis", "PC"),
        ("MDH_SHUTTLE_OUT", {"OAA_m": 1.0, "NADH_m": 1.0}, {"OAA_c": 1.0, "NADH_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "MDH2 and MDH1 and SLC25A11"),
        ("PEPCK_C", {"OAA_c": 1.0, "ATP_c": 1.0}, {"PEP_c": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "PCK1"),
        ("PEPCK_M", {"OAA_m": 1.0, "ATP_m": 1.0}, {"PEP_c": 1.0, "ADP_m": 1.0}, 0.0, 10.0, "Gluconeogenesis", "PCK2"),

        # --- Pentose Phosphate Pathway (NADPH supply for Lipogenesis) ---
        ("G6PDH", {"G6P_c": 1.0, "NADP_c": 1.0}, {"P6G_c": 1.0, "NADPH_c": 1.0}, 0.0, 10.0, "PPP", "G6PDX and PGLS"),
        ("PGD", {"P6G_c": 1.0, "NADP_c": 1.0}, {"RU5P_c": 1.0, "NADPH_c": 1.0}, 0.0, 10.0, "PPP", "PGD"),
        ("TKT_TALDO", {"RU5P_c": 1.0}, {"F6P_c": 1.0}, 0.0, 10.0, "PPP", "TKT and TALDO1"),
        ("ME1", {"PYR_c": 1.0, "NADP_c": 1.0}, {"ACCOA_c": 1.0, "NADPH_c": 1.0}, 0.0, 10.0, "Lipogenesis", "ME1"),

        # --- PDH, Citrate Shuttle & De Novo Lipogenesis ---
        ("PDH", {"PYR_m": 1.0, "NAD_m": 1.0}, {"ACCOA_m": 1.0, "NADH_m": 1.0}, 0.0, 10.0, "TCA_OxPhos", "PDHA1 and PDHB and DLD"),
        ("CS", {"ACCOA_m": 1.0, "OAA_m": 1.0}, {"CIT_m": 1.0}, 0.0, 10.0, "TCA_OxPhos", "CS"),
        ("CL_SHUTTLE", {"CIT_m": 1.0, "ATP_c": 1.0}, {"ACCOA_c": 1.0, "OAA_m": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Lipogenesis", "SLC25A1 and ACLY"),
        ("ACC", {"ACCOA_c": 1.0, "ATP_c": 1.0}, {"MALCOA_c": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Lipogenesis", "ACACA or ACACB"),
        ("FASN", {"ACCOA_c": 1.0, "MALCOA_c": 1.0, "NADPH_c": 2.0}, {"PALM_c": 1.0, "NADP_c": 2.0}, 0.0, 10.0, "Lipogenesis", "FASN"),

        # --- Mitochondrial Beta-Oxidation & Ketogenesis ---
        ("CPT1_CPT2", {"FA_c": 1.0, "ATP_c": 1.0}, {"FACOA_m": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Beta_Oxidation", "CPT1A and CPT2"),
        ("ACADM_STEP", {"FACOA_m": 1.0, "FAD_m": 1.0}, {"ENOYL_m": 1.0, "FADH2_m": 1.0}, 0.0, 10.0, "Beta_Oxidation", "ACADVL or ACADM or ACADS"),
        ("MTP_STEP", {"ENOYL_m": 1.0, "NAD_m": 1.0}, {"ACCOA_m": 2.0, "NADH_m": 1.0}, 0.0, 10.0, "Beta_Oxidation", "HADHA and HADHB"),
        ("HMGCS2", {"ACCOA_m": 2.0}, {"HMGCOA_m": 1.0}, 0.0, 10.0, "Ketogenesis", "HMGCS2"),
        ("HMGCL_BDH1", {"HMGCOA_m": 1.0, "NADH_m": 1.0}, {"BHB_m": 1.0, "NAD_m": 1.0}, 0.0, 10.0, "Ketogenesis", "HMGCL and BDH1"),

        # --- TCA Cycle & Oxidative Phosphorylation ---
        ("TCA_OX", {"CIT_m": 1.0, "NAD_m": 2.0, "FAD_m": 1.0, "ADP_m": 1.0}, {"OAA_m": 1.0, "NADH_m": 2.0, "FADH2_m": 1.0, "ATP_m": 1.0}, 0.0, 10.0, "TCA_OxPhos", "IDH2 and OGDH and SDHA"),
        ("ETC_NADH", {"NADH_m": 1.0, "ADP_m": 2.0}, {"NAD_m": 1.0, "ATP_m": 2.0}, 0.0, 15.0, "TCA_OxPhos", "NDUFS1 and UQCRC1 and COX4I1"),
        ("ETC_FADH2", {"FADH2_m": 1.0, "ADP_m": 1.0}, {"FAD_m": 1.0, "ATP_m": 1.0}, 0.0, 15.0, "TCA_OxPhos", "SDHA and SDHB and COX4I1"),
        ("ANT_ATP", {"ATP_m": 1.0, "ADP_c": 1.0}, {"ATP_c": 1.0, "ADP_m": 1.0}, -15.0, 15.0, "Transport", "SLC25A4 or SLC25A5"),
        ("MAS_SHUTTLE", {"NADH_c": 1.0, "NAD_m": 1.0}, {"NAD_c": 1.0, "NADH_m": 1.0}, 0.0, 10.0, "Transport", "GOT1 and GOT2 and SLC25A12"),
        ("G3P_SHUTTLE", {"NADH_c": 1.0, "FAD_m": 1.0}, {"NAD_c": 1.0, "FADH2_m": 1.0}, 0.0, 10.0, "Transport", "GPD1 and GPD2"),

        # --- Decoy Branches & Potential Futile Loops ---
        ("GLYG_SYN", {"G6P_c": 1.0, "ATP_c": 1.0}, {"GLYG_c": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Glycogen", "GYS2 and UGP2"),
        ("GLYG_DEG", {"GLYG_c": 1.0}, {"G6P_c": 1.0}, 0.0, 10.0, "Glycogen", "PYGL"),
        ("GFAT_HBP", {"F6P_c": 1.0, "ATP_c": 1.0}, {"UDPGNAc_c": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Decoy_Hexosamine", "GFPT1 and UAP1"),
        ("FUTILE_ATPASE", {"ATP_c": 1.0}, {"ADP_c": 1.0}, 0.0, 10.0, "Decoy_Futile", "ENTPD5"),
        ("NNT_FORWARD", {"NADH_m": 1.0, "NADP_c": 1.0}, {"NAD_m": 1.0, "NADPH_c": 1.0}, 0.0, 10.0, "Decoy_Futile", "NNT"),
        ("NNT_REVERSE", {"NADPH_c": 1.0, "NAD_m": 1.0}, {"NADP_c": 1.0, "NADH_m": 1.0}, 0.0, 10.0, "Decoy_Futile", "NNT_REV"),
    ]

    graph = cn.Graph()
    for rxn_id, reactants, products, lb, ub, pathway, gpr in specs:
        graph.add_edge(
            reactants,
            products,
            id=rxn_id,
            name=rxn_id,
            default_lb=float(lb),
            default_ub=float(ub),
            pathway=pathway,
            GPR=gpr,
            reversible=bool(lb < 0),
        )

    fed_task = MetabolicTask(
        name="fed_glycolytic_lipogenic",
        bounds={
            "EX_GLC_IN": (0.0, 10.0),
            "EX_GLC_OUT": (0.0, 0.0),
            "EX_FA_IN": (0.0, 0.0),
            "EX_LAC_IN": (0.0, 0.0),
            "EX_ALA_IN": (0.0, 0.0),
            "EX_LIPID_OUT": (1.0, 1.0),
            "EX_KET_OUT": (0.0, 0.0),
            "EX_GLYCOGEN_OUT": (0.0, 0.0),
            "EX_HEXOSAMINE_OUT": (0.0, 0.0),
            "ATP_MAINT": (2.0, 2.0),
        },
    )
    fasted_task = MetabolicTask(
        name="fasted_gluconeogenic_ketogenic",
        bounds={
            "EX_GLC_IN": (0.0, 0.0),
            "EX_GLC_OUT": (1.0, 1.0),
            "EX_FA_IN": (0.0, 10.0),
            "EX_LAC_IN": (0.0, 10.0),
            "EX_ALA_IN": (0.0, 10.0),
            "EX_LIPID_OUT": (0.0, 0.0),
            "EX_KET_OUT": (1.0, 1.0),
            "EX_GLYCOGEN_OUT": (0.0, 0.0),
            "EX_HEXOSAMINE_OUT": (0.0, 0.0),
            "ATP_MAINT": (2.0, 2.0),
        },
    )
    tasks = (fed_task, fasted_task)
    graph.get_graph_attributes()["default_benchmark_tasks"] = tasks
    return graph, tasks
