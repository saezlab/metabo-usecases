"""Small, reproducible simulator for paired RNA-seq and proteomics data.

The simulator creates a shared tissue scaffold plus independent modality
programs for two tissues, two conditions, and an ordered set of timepoints.
Arrays use genes on rows and states or samples on columns.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass(frozen=True)
class _ScientificAssumptions:
    """Named defaults for fixed scientific assumptions in simulator-v1."""

    active_loading_sd: float = 1.0
    background_loading_sd: float = 0.04
    rna_baseline_sd: float = 0.7
    protein_baseline_sd: float = 1.5
    temporal_bump_width: float = 2.0
    library_depth_log_sd: float = 0.08
    dispersion_log_sd: float = 0.35


_COMPONENTS = (
    "tissue",
    "condition",
    "time",
    "tissue_condition",
    "tissue_time",
    "condition_time",
    "tissue_condition_time",
)
_DEFAULT_VARIANCE_WEIGHTS = {
    "tissue": 0.70,
    "condition": 0.03,
    "time": 0.15,
    "tissue_condition": 0.02,
    "tissue_time": 0.05,
    "condition_time": 0.03,
    "tissue_condition_time": 0.02,
}
# Each row maps a latent program name to its variance component.
# The table also defines the returned program order and component membership.
_PROGRAM_DESIGN = (
    ("tissue", "tissue"),
    ("condition", "condition"),
    ("time_linear", "time"),
    ("time_quadratic", "time"),
    ("time_bump", "time"),
    ("tissue_condition", "tissue_condition"),
    ("tissue_time_linear", "tissue_time"),
    ("tissue_time_bump", "tissue_time"),
    ("condition_time_linear", "condition_time"),
    ("condition_time_bump", "condition_time"),
    ("tissue_condition_time", "tissue_condition_time"),
)
PROGRAM_NAMES = tuple(row[0] for row in _PROGRAM_DESIGN)
# Keep the internal spelling as an alias so existing implementation references
# remain easy to compare with the original simulator.
_PROGRAM_NAMES = PROGRAM_NAMES
_COMPONENT_PROGRAMS = {
    component: tuple(name for name, parent in _PROGRAM_DESIGN if parent == component)
    for component in _COMPONENTS
}

def _finite_nonnegative(name: str, value: float) -> float:
    result = float(value)
    if not np.isfinite(result) or result < 0:
        raise ValueError(f"{name} must be a finite non-negative number")
    return result


def _validate_gene_group_options(
    gene_groups: Mapping[str, Sequence[str]] | None,
    group_program_weights: Mapping[str, Mapping[str, float]] | None,
    group_strength: float,
    gene_ids: Sequence[str],
) -> dict[str, Any]:
    """Validate and canonicalize optional overlapping gene groups."""
    try:
        strength = float(group_strength)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("group_strength must be a finite number between zero and one") from exc
    if (
        isinstance(group_strength, (bool, np.bool_))
        or not np.isfinite(strength)
        or not 0 <= strength <= 1
    ):
        raise ValueError("group_strength must be a finite number between zero and one")

    if gene_groups is None:
        groups: dict[str, tuple[str, ...]] = {}
    elif not isinstance(gene_groups, Mapping):
        raise ValueError("gene_groups must be a mapping from group names to gene ID sequences")
    else:
        groups = {}
        known_genes = set(gene_ids)
        for name, members in gene_groups.items():
            if not isinstance(name, str) or not name.strip() or name != name.strip():
                raise ValueError(
                    "gene_groups keys must be non-empty names without surrounding whitespace"
                )
            if isinstance(members, (str, bytes)) or not isinstance(members, Sequence):
                raise ValueError(f"gene group {name!r} must contain a sequence of gene IDs")
            member_ids = tuple(members)
            if not member_ids:
                raise ValueError(f"gene group {name!r} must not be empty")
            if any(not isinstance(gene, str) or not gene.strip() for gene in member_ids):
                raise ValueError(f"gene group {name!r} contains an invalid gene ID")
            if len(set(member_ids)) != len(member_ids):
                raise ValueError(f"gene group {name!r} contains duplicate gene IDs")
            unknown = set(member_ids) - known_genes
            if unknown:
                raise ValueError(
                    f"gene group {name!r} contains unknown gene ID(s): {', '.join(sorted(unknown))}"
                )
            groups[name] = tuple(sorted(member_ids))
        groups = dict(sorted(groups.items()))

    if group_program_weights is not None and not groups:
        raise ValueError("group_program_weights requires at least one gene group")
    if group_program_weights is None:
        explicit_weights: dict[str, dict[str, float]] = {}
    elif not isinstance(group_program_weights, Mapping):
        raise ValueError("group_program_weights must map group names to program weight mappings")
    else:
        unknown_groups = set(group_program_weights) - set(groups)
        if unknown_groups:
            raise ValueError(
                "group_program_weights references unknown group(s): "
                + ", ".join(map(str, sorted(unknown_groups, key=str)))
            )
        explicit_weights = {}
        for group_name, weights in group_program_weights.items():
            if not isinstance(group_name, str) or not group_name.strip():
                raise ValueError("group_program_weights keys must be non-empty group names")
            if not isinstance(weights, Mapping):
                raise ValueError(
                    f"group_program_weights[{group_name!r}] must map program names to finite weights"
                )
            unknown_programs = set(weights) - set(PROGRAM_NAMES)
            if unknown_programs:
                raise ValueError(
                    f"group {group_name!r} references unknown program(s): "
                    + ", ".join(map(str, sorted(unknown_programs, key=str)))
                )
            normalized: dict[str, float] = {}
            for program_name, value in weights.items():
                if not isinstance(program_name, str):
                    raise ValueError(f"group {group_name!r} contains an invalid program name")
                try:
                    numeric = float(value)
                except (TypeError, ValueError, OverflowError) as exc:
                    raise ValueError(
                        f"group weight for {group_name!r}/{program_name!r} must be finite"
                    ) from exc
                if isinstance(value, (bool, np.bool_)) or not np.isfinite(numeric):
                    raise ValueError(
                        f"group weight for {group_name!r}/{program_name!r} must be finite"
                    )
                normalized[program_name] = numeric
            explicit_weights[group_name] = dict(sorted(normalized.items()))
        explicit_weights = dict(sorted(explicit_weights.items()))

    return {
        "gene_groups": groups or None,
        "group_program_weights": explicit_weights or None,
        "group_strength": strength,
    }


def _make_group_profiles(
    rng: np.random.Generator,
    design: dict[str, Any],
    parameters: dict[str, Any],
    assumptions: _ScientificAssumptions,
    options: dict[str, Any],
) -> dict[str, Any] | None:
    """Build normalized membership and modality-specific group coefficients."""
    groups = options["gene_groups"]
    if not groups:
        return None

    group_names = tuple(groups)
    n_genes = design["gene_count"]
    n_programs = len(design["program_names"])
    name_to_column = {name: i for i, name in enumerate(design["program_names"])}
    gene_to_row = {name: i for i, name in enumerate(design["gene_ids"])}
    membership = np.zeros((n_genes, len(group_names)), dtype=bool)
    for group_index, group_name in enumerate(group_names):
        membership[[gene_to_row[gene] for gene in groups[group_name]], group_index] = True
    row_counts = membership.sum(axis=1)
    basis = np.zeros(membership.shape, dtype=float)
    np.divide(
        membership,
        np.sqrt(np.maximum(row_counts, 1))[:, None],
        out=basis,
        where=row_counts[:, None] > 0,
    )

    explicit = options["group_program_weights"] or {}
    shared_coefficients = np.zeros((len(group_names), n_programs), dtype=float)
    rna_coefficients = np.zeros_like(shared_coefficients)
    protein_coefficients = np.zeros_like(shared_coefficients)
    shared_active = np.zeros_like(shared_coefficients, dtype=bool)
    rna_active = np.zeros_like(shared_coefficients, dtype=bool)
    protein_active = np.zeros_like(shared_coefficients, dtype=bool)
    tissue_column = name_to_column["tissue"]

    def random_component(columns: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        coefficients = np.zeros((len(group_names), n_programs), dtype=float)
        active = np.zeros_like(coefficients, dtype=bool)
        active[:, columns] = rng.random((len(group_names), len(columns))) < parameters[
            "module_probability"
        ]
        coefficients[:, columns] = rng.normal(
            0.0, assumptions.background_loading_sd, size=(len(group_names), len(columns))
        )
        coefficients[active] = rng.normal(
            0.0, assumptions.active_loading_sd, size=int(active.sum())
        )
        return coefficients, active

    non_tissue_columns = np.arange(n_programs)[np.arange(n_programs) != tissue_column]
    tissue_columns = np.asarray([tissue_column])
    random_shared, random_shared_active = random_component(tissue_columns)
    random_rna, random_rna_active = random_component(non_tissue_columns)
    random_protein, random_protein_active = random_component(non_tissue_columns)

    for group_index, group_name in enumerate(group_names):
        if group_name in explicit:
            for program_name, weight in explicit[group_name].items():
                column = name_to_column[program_name]
                if column == tissue_column:
                    shared_coefficients[group_index, column] = weight
                    shared_active[group_index, column] = weight != 0
                else:
                    rna_coefficients[group_index, column] = weight
                    protein_coefficients[group_index, column] = weight
                    rna_active[group_index, column] = weight != 0
                    protein_active[group_index, column] = weight != 0
        else:
            shared_coefficients[group_index] = random_shared[group_index]
            shared_active[group_index] = random_shared_active[group_index]
            rna_coefficients[group_index] = random_rna[group_index]
            rna_active[group_index] = random_rna_active[group_index]
            protein_coefficients[group_index] = random_protein[group_index]
            protein_active[group_index] = random_protein_active[group_index]

    return {
        "group_names": group_names,
        "membership": membership,
        "basis": basis,
        "grouped_genes": row_counts > 0,
        "coefficients": {
            "shared": shared_coefficients,
            "rna": rna_coefficients,
            "protein": protein_coefficients,
        },
        "active_profiles": {
            "shared": shared_active,
            "rna": rna_active,
            "protein": protein_active,
        },
    }


def _blend_group_loadings(
    individual: np.ndarray,
    individual_membership: np.ndarray,
    group_data: dict[str, Any] | None,
    component: str,
    strength: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Blend group and individual loadings before component normalization."""
    group_effective = np.zeros_like(individual)
    if group_data is None:
        return individual, group_effective, individual_membership
    coefficients = group_data["coefficients"][component]
    group_effective = group_data["basis"] @ coefficients
    if strength == 0:
        return individual, np.zeros_like(individual), individual_membership
    grouped = group_data["grouped_genes"]
    blended = individual.copy()
    blended[grouped] = (
        np.sqrt(1.0 - strength) * individual[grouped]
        + np.sqrt(strength) * group_effective[grouped]
    )
    active_profiles = group_data["active_profiles"][component]
    program_participation = (
        group_data["membership"].astype(np.int64)
        @ active_profiles.astype(np.int64)
        > 0
    )
    membership = individual_membership.copy()
    if strength == 1.0:
        membership[grouped] = False
    membership |= program_participation
    return blended, group_effective * np.sqrt(strength), membership


def _two_labels(name: str, labels: Sequence[str]) -> tuple[str, str]:
    if isinstance(labels, (str, bytes)):
        raise ValueError(f"{name} must contain exactly two distinct labels")
    try:
        values = tuple(labels)
    except TypeError as exc:
        raise ValueError(f"{name} must contain exactly two distinct labels") from exc
    if len(values) != 2 or any(not isinstance(v, str) or not v.strip() for v in values):
        raise ValueError(f"{name} must contain exactly two non-empty string labels")
    if values[0] == values[1]:
        raise ValueError(f"{name} labels must be distinct")
    return values  # type: ignore[return-value]


def _time_basis(
    timepoints: np.ndarray, bump_width: float
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    tau = 2.0 * (timepoints - timepoints[0]) / (timepoints[-1] - timepoints[0]) - 1.0
    bases = {
        "linear": tau,
        "quadratic": tau**2,
        "bump": np.exp(-bump_width * tau**2),
    }
    for key, basis in bases.items():
        basis = basis - basis.mean()
        sd = float(np.sqrt(np.mean(basis**2)))
        if sd <= np.finfo(float).eps:
            raise ValueError("timepoints must support non-constant smooth temporal bases")
        bases[key] = basis / sd
    return tau, bases


def _weight_values(
    variance_weights: Mapping[str, float] | Sequence[float] | None,
) -> dict[str, float]:
    if variance_weights is None:
        values = dict(_DEFAULT_VARIANCE_WEIGHTS)
    elif isinstance(variance_weights, Mapping):
        unknown = set(variance_weights) - set(_COMPONENTS)
        if unknown:
            raise ValueError(f"unknown variance component(s): {', '.join(sorted(unknown))}")
        values = dict(_DEFAULT_VARIANCE_WEIGHTS)
        values.update({key: float(value) for key, value in variance_weights.items()})
    else:
        sequence = tuple(variance_weights)
        if len(sequence) != len(_COMPONENTS):
            raise ValueError(f"variance_weights must have {len(_COMPONENTS)} values")
        values = dict(zip(_COMPONENTS, map(float, sequence), strict=True))

    for key, value in values.items():
        _finite_nonnegative(f"variance_weights[{key!r}]", value)
    total = sum(values.values())
    if not np.isclose(total, 1.0, rtol=0.0, atol=1e-8):
        raise ValueError("variance_weights must be non-negative and sum to one")
    return values


def _correlations_by_gene(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    valid = np.isfinite(left) & np.isfinite(right)
    counts = valid.sum(axis=1)
    safe_counts = np.maximum(counts, 1)
    left_mean = np.sum(np.where(valid, left, 0.0), axis=1) / safe_counts
    right_mean = np.sum(np.where(valid, right, 0.0), axis=1) / safe_counts
    left_centered = np.where(valid, left - left_mean[:, None], 0.0)
    right_centered = np.where(valid, right - right_mean[:, None], 0.0)
    numerator = np.sum(left_centered * right_centered, axis=1)
    denominator = np.sqrt(np.sum(left_centered**2, axis=1) * np.sum(right_centered**2, axis=1))
    result = np.full(left.shape[0], np.nan, dtype=float)
    np.divide(numerator, denominator, out=result, where=(denominator > 0) & (counts >= 2))
    return np.clip(result, -1.0, 1.0)


def _correlation_summary(values: np.ndarray) -> dict[str, Any]:
    finite = values[np.isfinite(values)]
    return {
        "per_gene": values,
        "n_finite": int(finite.size),
        "mean": float(np.mean(finite)) if finite.size else float("nan"),
        "median": float(np.median(finite)) if finite.size else float("nan"),
        "10th_percentile": float(np.quantile(finite, 0.1)) if finite.size else float("nan"),
        "90th_percentile": float(np.quantile(finite, 0.9)) if finite.size else float("nan"),
    }


def _state_means(values: np.ndarray, sample_state_indices: np.ndarray, n_states: int) -> np.ndarray:
    """Average gene-by-sample values within each state, ignoring missing cells."""
    means = []
    for state_index in range(n_states):
        state_values = values[:, sample_state_indices == state_index]
        valid = np.isfinite(state_values)
        count = valid.sum(axis=1)
        total = np.sum(np.where(valid, state_values, 0.0), axis=1)
        state_mean = np.full(values.shape[0], np.nan)
        np.divide(total, count, out=state_mean, where=count > 0)
        means.append(state_mean)
    return np.stack(means, axis=1)


def _replicate_variance_summary(
    values: np.ndarray, sample_state_indices: np.ndarray, n_states: int
) -> dict[str, float]:
    by_state = [values[:, sample_state_indices == state] for state in range(n_states)]
    replicate_stack = np.stack(by_state, axis=1)
    means = _state_means(values, sample_state_indices, n_states)

    def row_variance(matrix: np.ndarray, axis: int, ddof: int) -> np.ndarray:
        valid = np.isfinite(matrix)
        counts = valid.sum(axis=axis)
        safe_counts = np.maximum(counts, 1)
        average = np.sum(np.where(valid, matrix, 0.0), axis=axis) / safe_counts
        expanded_average = np.expand_dims(average, axis=axis)
        squared = np.sum(
            np.where(valid, (matrix - expanded_average) ** 2, 0.0), axis=axis
        )
        variance = np.full(average.shape, np.nan, dtype=float)
        np.divide(squared, counts - ddof, out=variance, where=counts > ddof)
        return variance

    within_variance = row_variance(
        replicate_stack, axis=2, ddof=1 if replicate_stack.shape[2] > 1 else 0
    )
    between_variance = row_variance(means, axis=1, ddof=0)
    finite_within = within_variance[np.isfinite(within_variance)]
    finite_between = between_variance[np.isfinite(between_variance)]
    within = float(np.mean(finite_within)) if finite_within.size else float("nan")
    between = float(np.mean(finite_between)) if finite_between.size else float("nan")
    return {
        "mean_within_state_variance": within,
        "mean_between_state_variance": between,
        "within_fraction": within / (within + between) if within + between > 0 else float("nan"),
    }


def _standardize_vector(values: np.ndarray) -> np.ndarray:
    centered = values - values.mean()
    scale = float(np.sqrt(np.mean(centered**2)))
    return centered / scale if scale > np.finfo(float).eps else np.zeros_like(centered)


def _state_variance_fractions(values: np.ndarray, n_times: int) -> dict[str, float]:
    """Balanced tissue-condition-time ANOVA fractions, excluding gene baselines."""
    centered = values - values.mean(axis=1, keepdims=True)
    grid = centered.reshape(values.shape[0], 2, 2, n_times)
    tissue = grid.mean(axis=(2, 3), keepdims=True)
    condition = grid.mean(axis=(1, 3), keepdims=True)
    time = grid.mean(axis=(1, 2), keepdims=True)
    tissue_condition = grid.mean(axis=3, keepdims=True) - tissue - condition
    tissue_time = grid.mean(axis=2, keepdims=True) - tissue - time
    condition_time = grid.mean(axis=1, keepdims=True) - condition - time
    three_way = centered.reshape(values.shape[0], 2, 2, n_times) - (
        tissue + condition + time + tissue_condition + tissue_time + condition_time
    )
    components = {
        "tissue": tissue,
        "condition": condition,
        "time": time,
        "tissue_condition": tissue_condition,
        "tissue_time": tissue_time,
        "condition_time": condition_time,
        "tissue_condition_time": three_way,
    }
    total = float(np.mean(centered**2))
    if total <= np.finfo(float).eps:
        return {name: 0.0 for name in _COMPONENTS}
    return {name: float(np.mean(component**2) / total) for name, component in components.items()}


def _make_contrasts(
    shared: np.ndarray,
    tissues: tuple[str, str],
    conditions: tuple[str, str],
    n_times: int,
) -> dict[str, np.ndarray]:
    by_state = shared.reshape(shared.shape[0], 2, 2, n_times)
    condition_diff = by_state[:, :, 1, :] - by_state[:, :, 0, :]
    tissue_diff = by_state[:, 1, :, :] - by_state[:, 0, :, :]
    adjacent_time_diff = np.diff(by_state, axis=3)
    condition_by_time = condition_diff - condition_diff[:, :, [0]]
    return {
        "condition_within_tissue_time": condition_diff,
        "tissue_within_condition_time": tissue_diff,
        "adjacent_time": adjacent_time_diff,
        "condition_by_time_from_first": condition_by_time,
        "tissue_labels": tissues,
        "condition_labels": conditions,
    }


def _make_design(
    *,
    gene_count: int | None,
    gene_ids: Sequence[str] | None,
    tissues: Sequence[str],
    conditions: Sequence[str],
    timepoints: Sequence[float] | None,
    replicates: int,
    assumptions: _ScientificAssumptions,
) -> dict[str, Any]:
    """Validate and expand the tissue → condition → time sample design."""
    tissue_labels = _two_labels("tissues", tissues)
    condition_labels = _two_labels("conditions", conditions)

    if gene_count is None:
        n_genes = len(gene_ids) if gene_ids is not None else 1000
    else:
        try:
            gene_count_value = float(gene_count)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("gene_count must be a positive integer") from exc
        if (
            isinstance(gene_count, (bool, np.bool_))
            or not np.isfinite(gene_count_value)
            or not gene_count_value.is_integer()
        ):
            raise ValueError("gene_count must be a positive integer")
        n_genes = int(gene_count_value)
    if n_genes < 1:
        raise ValueError("gene_count must be a positive integer")
    if gene_ids is None:
        identifiers = tuple(f"gene_{i:04d}" for i in range(n_genes))
    else:
        identifiers = tuple(gene_ids)
        if len(identifiers) != n_genes:
            raise ValueError("gene_ids length must match gene_count")
        if any(not isinstance(gene, str) or not gene.strip() for gene in identifiers):
            raise ValueError("gene_ids must contain non-empty strings")
        if len(set(identifiers)) != len(identifiers):
            raise ValueError("gene_ids must be unique")

    try:
        replicate_value = float(replicates)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("replicates must be a positive integer") from exc
    if (
        isinstance(replicates, (bool, np.bool_))
        or not np.isfinite(replicate_value)
        or not replicate_value.is_integer()
        or replicate_value < 1
    ):
        raise ValueError("replicates must be a positive integer")
    n_replicates = int(replicate_value)

    times = np.linspace(0.0, 1.0, 8) if timepoints is None else np.asarray(timepoints, dtype=float)
    if times.ndim != 1 or len(times) < 3 or not np.all(np.isfinite(times)):
        raise ValueError("timepoints must be a finite one-dimensional sequence of at least 3 values")
    if np.any(np.diff(times) <= 0):
        raise ValueError("timepoints must be strictly increasing and non-constant")
    _, temporal_bases = _time_basis(times, assumptions.temporal_bump_width)

    n_times = len(times)
    n_states = 2 * 2 * n_times
    n_samples = n_states * n_replicates
    state_metadata = [
        {
            "state_index": state_index,
            "tissue": tissue,
            "condition": condition,
            "timepoint": float(time),
            "time_index": time_index,
        }
        for state_index, (tissue, condition, time_index, time) in enumerate(
            (tissue, condition, time_index, time)
            for tissue in tissue_labels
            for condition in condition_labels
            for time_index, time in enumerate(times)
        )
    ]

    # Build the design vectors directly in the same ordered state grid.
    tissue_values = np.repeat(np.asarray([-1.0, 1.0]), 2 * n_times)
    condition_values = np.tile(np.repeat(np.asarray([-1.0, 1.0]), n_times), 2)
    time_indices = np.tile(np.arange(n_times), 4)
    activity_vectors = {
        "tissue": tissue_values,
        "condition": condition_values,
        "time_linear": temporal_bases["linear"][time_indices],
        "time_quadratic": temporal_bases["quadratic"][time_indices],
        "time_bump": temporal_bases["bump"][time_indices],
        "tissue_condition": tissue_values * condition_values,
        "tissue_time_linear": tissue_values * temporal_bases["linear"][time_indices],
        "tissue_time_bump": tissue_values * temporal_bases["bump"][time_indices],
        "condition_time_linear": condition_values * temporal_bases["linear"][time_indices],
        "condition_time_bump": condition_values * temporal_bases["bump"][time_indices],
        "tissue_condition_time": (
            tissue_values * condition_values * temporal_bases["bump"][time_indices]
        ),
    }
    activities = np.vstack([activity_vectors[name] for name in _PROGRAM_NAMES])
    activities -= activities.mean(axis=1, keepdims=True)

    return {
        "gene_count": n_genes,
        "gene_ids": identifiers,
        "tissues": tissue_labels,
        "conditions": condition_labels,
        "timepoints": times,
        "n_times": n_times,
        "n_replicates": n_replicates,
        "n_states": n_states,
        "n_samples": n_samples,
        "state_metadata": state_metadata,
        "temporal_bases": temporal_bases,
        "time_indices": time_indices,
        "activity_vectors": activity_vectors,
        "program_names": _PROGRAM_NAMES,
        "program_activities": activities,
    }


def _validate_parameters(
    *,
    design: dict[str, Any],
    seed: int,
    variance_weights: Mapping[str, float] | Sequence[float] | None,
    effect_scale: float,
    module_probability: float,
    main_time_secondary_loading_scale: float,
    coupling_mean: float,
    coupling_sd: float,
    coupling_range: tuple[float, float],
    baseline_coupling: float,
    biological_noise: float,
    replicate_noise: float,
    modality_response_noise: float,
    replicate_activity_noise: float,
    replicate_tissue_activity_noise: float,
    library_depth: float | Sequence[float],
    dispersion: float,
    proteomics_noise: float,
    missingness: bool,
    missingness_intercept: float,
    missingness_slope: float,
) -> dict[str, Any]:
    """Normalize public controls and reject invalid simulation parameters."""
    weights = _weight_values(variance_weights)
    effect_scale = _finite_nonnegative("effect_scale", effect_scale)
    main_time_secondary_loading_scale = _finite_nonnegative(
        "main_time_secondary_loading_scale", main_time_secondary_loading_scale
    )
    if main_time_secondary_loading_scale > 1:
        raise ValueError("main_time_secondary_loading_scale must be between zero and one")
    biological_noise = _finite_nonnegative("biological_noise", biological_noise)
    replicate_noise = _finite_nonnegative("replicate_noise", replicate_noise)
    modality_response_noise = _finite_nonnegative(
        "modality_response_noise", modality_response_noise
    )
    replicate_activity_noise = _finite_nonnegative(
        "replicate_activity_noise", replicate_activity_noise
    )
    replicate_tissue_activity_noise = _finite_nonnegative(
        "replicate_tissue_activity_noise", replicate_tissue_activity_noise
    )
    proteomics_noise = _finite_nonnegative("proteomics_noise", proteomics_noise)
    dispersion = _finite_nonnegative("dispersion", dispersion)
    if not np.isfinite(module_probability) or not 0 <= module_probability <= 1:
        raise ValueError("module_probability must be between zero and one")
    if not np.isfinite(coupling_mean) or not 0 <= coupling_mean <= 1:
        raise ValueError("coupling_mean must be between zero and one")
    coupling_sd = _finite_nonnegative("coupling_sd", coupling_sd)
    try:
        coupling_range = tuple(coupling_range)
    except TypeError as exc:
        raise ValueError("coupling_range must be a (minimum, maximum) pair") from exc
    if len(coupling_range) != 2:
        raise ValueError("coupling_range must be a (minimum, maximum) pair")
    coupling_low, coupling_high = map(float, coupling_range)
    if (
        not np.isfinite(coupling_low)
        or not np.isfinite(coupling_high)
        or coupling_low < 0
        or coupling_high > 1
        or coupling_low >= coupling_high
    ):
        raise ValueError("coupling_range must satisfy 0 <= minimum < maximum <= 1")
    if not coupling_low <= coupling_mean <= coupling_high:
        raise ValueError("coupling_mean must fall within coupling_range")
    if not np.isfinite(baseline_coupling) or not 0 <= baseline_coupling <= 1:
        raise ValueError("baseline_coupling must be between zero and one")
    if not isinstance(missingness, (bool, np.bool_)):
        raise ValueError("missingness must be a boolean")
    if not np.isfinite(missingness_intercept):
        raise ValueError("missingness_intercept must be finite")
    missingness_slope = _finite_nonnegative("missingness_slope", missingness_slope)

    if np.isscalar(library_depth):
        expected_depths = np.full(design["n_samples"], float(library_depth), dtype=float)
    else:
        expected_depths = np.asarray(library_depth, dtype=float)
        if expected_depths.ndim != 1 or len(expected_depths) != design["n_samples"]:
            raise ValueError("library_depth must be a positive scalar or one value per sample")
    if not np.all(np.isfinite(expected_depths)) or np.any(expected_depths <= 0):
        raise ValueError("library_depth values must be finite and positive")

    return {
        "seed": seed,
        "variance_weights": weights,
        "effect_scale": effect_scale,
        "module_probability": module_probability,
        "main_time_secondary_loading_scale": main_time_secondary_loading_scale,
        "coupling_mean": coupling_mean,
        "coupling_sd": coupling_sd,
        "coupling_low": coupling_low,
        "coupling_high": coupling_high,
        "baseline_coupling": baseline_coupling,
        "biological_noise": biological_noise,
        "replicate_noise": replicate_noise,
        "modality_response_noise": modality_response_noise,
        "replicate_activity_noise": replicate_activity_noise,
        "replicate_tissue_activity_noise": replicate_tissue_activity_noise,
        "library_depth_expected": expected_depths,
        "dispersion": dispersion,
        "proteomics_noise": proteomics_noise,
        "missingness": missingness,
        "missingness_intercept": missingness_intercept,
        "missingness_slope": missingness_slope,
    }


def _generate_shared_biology(
    rng: np.random.Generator,
    design: dict[str, Any],
    parameters: dict[str, Any],
    assumptions: _ScientificAssumptions,
    group_data: dict[str, Any] | None,
) -> dict[str, Any]:
    """Generate the common tissue program and gene baselines."""
    n_genes = design["gene_count"]
    n_programs = len(design["program_names"])
    tissue_column = design["program_names"].index("tissue")
    module_membership = np.zeros((n_genes, n_programs), dtype=bool)
    module_membership[:, tissue_column] = (
        rng.random(n_genes) < parameters["module_probability"]
    )
    program_loadings = np.zeros((n_genes, n_programs), dtype=float)
    program_loadings[:, tissue_column] = rng.normal(
        0.0, assumptions.background_loading_sd, size=n_genes
    )
    active = module_membership[:, tissue_column]
    program_loadings[active, tissue_column] = rng.normal(
        0.0, assumptions.active_loading_sd, size=int(active.sum())
    )

    program_loadings, group_loadings, module_membership = _blend_group_loadings(
        program_loadings,
        module_membership,
        group_data,
        "shared",
        parameters["group_strength"],
    )

    tissue_activity = design["program_activities"][tissue_column]
    raw_tissue = program_loadings[:, tissue_column, None] * tissue_activity[None, :]
    raw_tissue -= raw_tissue.mean(axis=1, keepdims=True)
    raw_rms = float(np.sqrt(np.mean(raw_tissue**2)))
    requested_rms = parameters["effect_scale"] * np.sqrt(
        parameters["variance_weights"]["tissue"]
    )
    multiplier = requested_rms / raw_rms if raw_rms > np.finfo(float).eps else 0.0
    program_loadings[:, tissue_column] *= multiplier
    group_loadings[:, tissue_column] *= multiplier
    shared_dynamic = program_loadings[:, tissue_column, None] * tissue_activity[None, :]
    shared_baseline = rng.normal(0.0, assumptions.rna_baseline_sd, size=n_genes)
    shared_truth = shared_baseline[:, None] + shared_dynamic
    return {
        "shared": shared_truth,
        "shared_baseline": shared_baseline,
        "shared_dynamic": shared_dynamic,
        "program_loadings": program_loadings,
        "module_membership": module_membership,
        "group_effective_loadings": group_loadings,
        "scaled_components": {"tissue": shared_dynamic},
    }


def _modality_activities(
    rng: np.random.Generator,
    design: dict[str, Any],
    response_noise: float,
) -> np.ndarray:
    """Give one modality its own smooth temporal response for each program."""
    activities = design["program_activities"].copy()
    vectors = design["activity_vectors"]
    bases = design["temporal_bases"]
    times = design["time_indices"]
    tissue = vectors["tissue"]
    condition = vectors["condition"]
    modifiers = {
        "time": np.ones_like(tissue),
        "tissue_time": tissue,
        "condition_time": condition,
        "tissue_condition_time": tissue * condition,
    }
    for program_index, name in enumerate(design["program_names"]):
        if "time" not in name:
            continue
        if name.startswith("tissue_time"):
            prefix, basis_name = "tissue_time", name.removeprefix("tissue_time_")
        elif name.startswith("condition_time"):
            prefix, basis_name = "condition_time", name.removeprefix("condition_time_")
        elif name == "tissue_condition_time":
            prefix, basis_name = name, "bump"
        else:
            prefix, basis_name = "time", name.removeprefix("time_")
        temporal_shape = bases[basis_name].copy()
        if response_noise:
            other_bases = [basis for key, basis in bases.items() if key != basis_name]
            temporal_shape += response_noise * sum(
                rng.normal() * basis for basis in other_bases
            ) / np.sqrt(len(other_bases))
        modifier = modifiers[prefix]
        values = modifier * temporal_shape[times]
        values -= values.mean()
        scale = float(np.sqrt(np.mean(values**2)))
        if scale > np.finfo(float).eps:
            activities[program_index] = values / scale
    return activities


def _generate_modality_biology(
    rng: np.random.Generator,
    design: dict[str, Any],
    parameters: dict[str, Any],
    assumptions: _ScientificAssumptions,
    group_data: dict[str, Any] | None,
    modality: str,
) -> dict[str, Any]:
    """Generate independent overlapping programs for RNA or protein.

    For genes active in multiple main-time modules, one active loading is
    sampled as primary and the other active main-time loadings are scaled by
    ``main_time_secondary_loading_scale`` before component normalization.
    """
    n_genes = design["gene_count"]
    program_names = design["program_names"]
    activities = _modality_activities(
        rng, design, parameters["modality_response_noise"]
    )
    membership = rng.random((n_genes, len(program_names))) < parameters["module_probability"]
    loadings = rng.normal(
        0.0, assumptions.background_loading_sd, size=(n_genes, len(program_names))
    )
    loadings[membership] = rng.normal(
        0.0, assumptions.active_loading_sd, size=int(membership.sum())
    )
    main_time_columns = np.asarray(
        [program_names.index(name) for name in ("time_linear", "time_quadratic", "time_bump")],
        dtype=int,
    )
    secondary_scale = parameters["main_time_secondary_loading_scale"]
    for gene in range(n_genes):
        active_columns = main_time_columns[membership[gene, main_time_columns]]
        if len(active_columns) > 1:
            primary = int(rng.choice(active_columns))
            secondary = active_columns[active_columns != primary]
            loadings[gene, secondary] *= secondary_scale
    tissue_column = program_names.index("tissue")
    # Tissue identity is the deliberately shared program; all other columns
    # are independently drawn for each modality.
    loadings[:, tissue_column] = 0.0
    membership[:, tissue_column] = False
    loadings, group_loadings, membership = _blend_group_loadings(
        loadings,
        membership,
        group_data,
        modality,
        parameters["group_strength"],
    )
    scaled_components: dict[str, np.ndarray] = {}
    for component in _COMPONENTS:
        if component == "tissue":
            continue
        columns = np.asarray(
            [program_names.index(name) for name in _COMPONENT_PROGRAMS[component]],
            dtype=int,
        )
        raw = loadings[:, columns] @ activities[columns]
        raw -= raw.mean(axis=1, keepdims=True)
        raw_rms = float(np.sqrt(np.mean(raw**2)))
        requested_rms = parameters["effect_scale"] * np.sqrt(
            parameters["variance_weights"][component]
        )
        multiplier = requested_rms / raw_rms if raw_rms > np.finfo(float).eps else 0.0
        loadings[:, columns] *= multiplier
        group_loadings[:, columns] *= multiplier
        scaled_components[component] = loadings[:, columns] @ activities[columns]
    dynamic = np.sum(np.stack(tuple(scaled_components.values()), axis=0), axis=0)
    return {
        "dynamic": dynamic,
        "program_names": program_names,
        "program_activities": activities,
        "program_loadings": loadings,
        "module_membership": membership,
        "group_effective_loadings": group_loadings,
        "scaled_components": scaled_components,
    }


def _generate_rna(
    rng: np.random.Generator,
    shared: dict[str, Any],
    rna_biology: dict[str, Any],
    noise_sd: float,
) -> dict[str, Any]:
    """Generate latent RNA as ``R = μ + T + D_R + ε_R``."""
    epsilon_R = rng.normal(0.0, noise_sd, size=shared["shared"].shape)
    R = shared["shared"] + rna_biology["dynamic"] + epsilon_R
    return {"rna": R, "rna_state_noise": epsilon_R}


def _generate_protein(
    rng: np.random.Generator,
    rna: dict[str, Any],
    rna_biology: dict[str, Any],
    shared: dict[str, Any],
    protein_biology: dict[str, Any],
    design: dict[str, Any],
    parameters: dict[str, Any],
    assumptions: _ScientificAssumptions,
) -> dict[str, Any]:
    r"""Construct protein with shared tissue and independent program biology.

    Let ``T`` be the shared tissue dynamic, ``D_R`` and ``D_P`` be the
    independent non-tissue RNA and protein program dynamics, and ``ε_R`` be
    RNA-specific state noise. Then

    ``b ~ clipped Normal(coupling_mean, coupling_sd)``
    ``P = a + T + b * (D_R + ε_R) + sqrt(1 - b**2) * D_P``

    Protein-specific program loadings include the final per-gene
    ``sqrt(1 - b**2)`` factor, so their product with their activities
    reconstructs the independent protein term. ``b`` is an RNA mixing weight,
    not a universal correlation target; the diagnostics report observed
    within-tissue and pooled correlations separately.
    """
    n_genes = design["gene_count"]
    rna_weight = np.clip(
        rng.normal(
            parameters["coupling_mean"], parameters["coupling_sd"], size=n_genes
        ),
        parameters["coupling_low"],
        parameters["coupling_high"],
    )
    residual_weight = np.sqrt(1.0 - rna_weight**2)
    protein_program_loadings = protein_biology["program_loadings"] * residual_weight[:, None]
    protein_group_loadings = protein_biology["group_effective_loadings"] * residual_weight[:, None]
    protein_specific = protein_program_loadings @ protein_biology["program_activities"]
    protein_scaled_components = {}
    for component in _COMPONENTS:
        if component == "tissue":
            continue
        columns = np.asarray(
            [
                design["program_names"].index(name)
                for name in _COMPONENT_PROGRAMS[component]
            ],
            dtype=int,
        )
        protein_scaled_components[component] = (
            protein_program_loadings[:, columns]
            @ protein_biology["program_activities"][columns]
        )
    protein_rna_component = (
        shared["shared_dynamic"]
        + rna_weight[:, None]
        * (rna_biology["dynamic"] + rna["rna_state_noise"])
    )
    baseline_shared = _standardize_vector(shared["shared_baseline"])
    baseline_specific = _standardize_vector(rng.normal(size=n_genes))
    baseline_mix = (
        parameters["baseline_coupling"] * baseline_shared
        + np.sqrt(1.0 - parameters["baseline_coupling"]**2) * baseline_specific
    )
    protein_baseline = 20.0 + assumptions.protein_baseline_sd * _standardize_vector(baseline_mix)
    protein_truth = protein_baseline[:, None] + protein_rna_component + protein_specific
    return {
        "protein": protein_truth,
        "protein_specific": protein_specific,
        "protein_rna_component": protein_rna_component,
        "protein_baseline": protein_baseline,
        "protein_program_activities": protein_biology["program_activities"],
        "protein_program_loadings": protein_program_loadings,
        "protein_group_effective_loadings": protein_group_loadings,
        "protein_specific_variance_weights": {
            component: parameters["variance_weights"][component]
            for component in _COMPONENTS
            if component != "tissue"
        },
        "protein_program_names": protein_biology["program_names"],
        "protein_module_membership": protein_biology["module_membership"],
        "protein_scaled_components": protein_scaled_components,
        "coupling_values": rna_weight,
        "residual_weight": residual_weight,
    }


def _generate_replicates(
    rng: np.random.Generator,
    rna_truth: np.ndarray,
    rna_biology: dict[str, Any],
    protein_truth: np.ndarray,
    protein: dict[str, Any],
    shared: dict[str, Any],
    design: dict[str, Any],
    parameters: dict[str, Any],
    assumptions: _ScientificAssumptions,
) -> dict[str, Any]:
    """Expand state truths with shared and modality-specific ``δH`` shifts."""
    expected_depths = parameters["library_depth_expected"]
    depth_sigma = assumptions.library_depth_log_sd
    library_depths = expected_depths * np.exp(
        rng.normal(-0.5 * depth_sigma**2, depth_sigma, size=design["n_samples"])
    )
    sample_metadata: list[dict[str, Any]] = []
    state_indices: list[int] = []
    rna_observations: list[np.ndarray] = []
    protein_observations: list[np.ndarray] = []
    shared_activity_perturbations: list[np.ndarray] = []
    rna_activity_perturbations: list[np.ndarray] = []
    protein_activity_perturbations: list[np.ndarray] = []
    program_names = design["program_names"]
    tissue_column = program_names.index("tissue")
    activity_sd = parameters["replicate_activity_noise"]
    tissue_activity_sd = parameters["replicate_tissue_activity_noise"]
    common_tissue_loading = shared["program_loadings"][:, tissue_column]
    rna_loadings = rna_biology["program_loadings"].copy()
    rna_loadings[:, tissue_column] = 0.0
    protein_loadings = protein["protein_program_loadings"]
    coupling = protein["coupling_values"]
    for state in design["state_metadata"]:
        state_idx = state["state_index"]
        for replicate in range(design["n_replicates"]):
            sample_idx = len(sample_metadata)
            sample_metadata.append(
                {
                    "sample_index": sample_idx,
                    "sample_id": f"sample_{sample_idx:04d}",
                    "state_index": state_idx,
                    "tissue": state["tissue"],
                    "condition": state["condition"],
                    "timepoint": state["timepoint"],
                    "time_index": state["time_index"],
                    "replicate": replicate + 1,
                }
            )
            state_indices.append(state_idx)
            delta_tissue = rng.normal(0.0, tissue_activity_sd)
            delta_rna = rng.normal(0.0, activity_sd, size=len(program_names))
            delta_rna[tissue_column] = 0.0
            delta_protein = rng.normal(0.0, activity_sd, size=len(program_names))
            delta_protein[tissue_column] = 0.0
            common_shift = common_tissue_loading * delta_tissue
            rna_shift = rna_loadings @ delta_rna
            protein_shift = protein_loadings @ delta_protein
            shared_activity_perturbations.append(delta_tissue)
            rna_activity_perturbations.append(delta_rna)
            protein_activity_perturbations.append(delta_protein)
            rna_observations.append(
                rna_truth[:, state_idx]
                + common_shift
                + rna_shift
                + rng.normal(0.0, parameters["replicate_noise"], size=design["gene_count"])
            )
            protein_observations.append(
                protein_truth[:, state_idx]
                + common_shift
                + coupling * rna_shift
                + protein_shift
                + rng.normal(0.0, parameters["replicate_noise"], size=design["gene_count"])
            )
    return {
        "sample_metadata": sample_metadata,
        "sample_state_indices": np.asarray(state_indices, dtype=np.int64),
        "rna_observations": np.stack(rna_observations, axis=1),
        "protein_observations": np.stack(protein_observations, axis=1),
        "shared_activity_perturbations": np.asarray(shared_activity_perturbations),
        "rna_activity_perturbations": np.stack(rna_activity_perturbations, axis=1),
        "protein_activity_perturbations": np.stack(protein_activity_perturbations, axis=1),
        "library_depths": library_depths,
    }


def _generate_observations(
    rng: np.random.Generator,
    replicates: dict[str, Any],
    design: dict[str, Any],
    parameters: dict[str, Any],
    assumptions: _ScientificAssumptions,
) -> dict[str, Any]:
    """Generate negative-binomial RNA counts and noisy protein intensities."""
    rna_observations = replicates["rna_observations"]
    protein_latent_samples = replicates["protein_observations"]

    # Stable softmax, computed independently for every sample.
    shifted = rna_observations - np.max(rna_observations, axis=0, keepdims=True)
    abundances = np.exp(shifted)
    proportions = abundances / abundances.sum(axis=0, keepdims=True)
    expected_counts = proportions * replicates["library_depths"][None, :]

    dispersion = parameters["dispersion"]
    gene_dispersions = dispersion * rng.lognormal(
        mean=-0.5 * assumptions.dispersion_log_sd**2,
        sigma=assumptions.dispersion_log_sd,
        size=design["gene_count"],
    )
    if dispersion == 0:
        rna_counts = rng.poisson(expected_counts).astype(np.int64)
    else:
        size = 1.0 / gene_dispersions[:, None]
        probability = size / (size + expected_counts)
        rna_counts = rng.negative_binomial(size, probability).astype(np.int64)

    protein_log_intensities = protein_latent_samples + rng.normal(
        0.0,
        parameters["proteomics_noise"],
        size=(design["gene_count"], design["n_samples"]),
    )
    if parameters["missingness"]:
        standardized_abundance = (
            protein_log_intensities - float(np.mean(protein_log_intensities))
        ) / max(float(np.std(protein_log_intensities)), np.finfo(float).eps)
        missing_probability = 1.0 / (
            1.0
            + np.exp(
                -np.clip(
                    parameters["missingness_intercept"]
                    - parameters["missingness_slope"] * standardized_abundance,
                    -40,
                    40,
                )
            )
        )
        missingness_mask = (
            rng.random((design["gene_count"], design["n_samples"])) < missing_probability
        )
        protein_log_intensities[missingness_mask] = np.nan
    else:
        missingness_mask = np.zeros(
            (design["gene_count"], design["n_samples"]), dtype=bool
        )
    return {
        "rna_counts": rna_counts,
        "protein_log_intensities": protein_log_intensities,
        "gene_dispersions": gene_dispersions,
        "missingness_mask": missingness_mask,
    }


def _calculate_diagnostics(
    shared: dict[str, Any],
    rna: dict[str, Any],
    rna_biology: dict[str, Any],
    protein: dict[str, Any],
    replicates: dict[str, Any],
    observed: dict[str, Any],
    design: dict[str, Any],
    parameters: dict[str, Any],
) -> dict[str, Any]:
    """Summarize modality agreement and observed replicate variation."""
    shared_truth = shared["shared"]
    rna_truth = rna["rna"]
    protein_truth = protein["protein"]
    shared_baseline = shared["shared_baseline"]
    protein_baseline = protein["protein_baseline"]
    if np.std(shared_baseline) > 0 and np.std(protein_baseline) > 0:
        baseline_correlation = float(np.corrcoef(shared_baseline, protein_baseline)[0, 1])
    else:
        baseline_correlation = float("nan")
    n_times = design["n_times"]
    n_genes = design["gene_count"]

    def distances(values: np.ndarray) -> dict[str, float]:
        grid = values.reshape(n_genes, 2, 2, n_times)
        between = [
            np.sqrt(np.mean((grid[:, 1, condition, time] - grid[:, 0, condition, time]) ** 2))
            for condition in range(2)
            for time in range(n_times)
        ]
        within = [
            np.sqrt(np.mean((grid[:, tissue, 1, time] - grid[:, tissue, 0, time]) ** 2))
            for tissue in range(2)
            for time in range(n_times)
        ]
        return {
            "between_tissue": float(np.mean(between)),
            "within_tissue_condition": float(np.mean(within)),
        }

    state_variance_fractions = {
        "shared": _state_variance_fractions(shared_truth, n_times),
        "rna": _state_variance_fractions(rna_truth, n_times),
        "protein": _state_variance_fractions(protein_truth, n_times),
    }
    counts = observed["rna_counts"]
    library_totals = counts.sum(axis=0)
    rna_log_cpm = np.log1p(counts / np.maximum(library_totals, 1)[None, :] * 1_000_000.0)
    protein_samples = observed["protein_log_intensities"]
    sample_state_indices = replicates["sample_state_indices"]
    rna_state_means = _state_means(rna_log_cpm, sample_state_indices, design["n_states"])
    protein_state_means = _state_means(
        protein_samples, sample_state_indices, design["n_states"]
    )
    state_indices_by_tissue = {
        tissue: np.asarray(
            [state["state_index"] for state in design["state_metadata"] if state["tissue"] == tissue],
            dtype=int,
        )
        for tissue in design["tissues"]
    }
    correlation_diagnostics: dict[str, Any] = {}
    measures = {
        "latent_state": (rna_truth, protein_truth, False),
        "observed_replicates": (rna_log_cpm, protein_samples, True),
        "observed_state_means": (rna_state_means, protein_state_means, False),
    }
    for measure_name, (left, right, is_sample_level) in measures.items():
        per_tissue = {}
        for tissue, state_indices in state_indices_by_tissue.items():
            if is_sample_level:
                sample_mask = np.isin(sample_state_indices, state_indices)
                per_tissue[tissue] = _correlations_by_gene(
                    left[:, sample_mask], right[:, sample_mask]
                )
            else:
                per_tissue[tissue] = _correlations_by_gene(
                    left[:, state_indices], right[:, state_indices]
                )
        correlation_diagnostics[measure_name] = {
            "pooled": _correlation_summary(_correlations_by_gene(left, right)),
            "within_tissue": {
                tissue: _correlation_summary(values) for tissue, values in per_tissue.items()
            },
            "within_tissue_combined": _correlation_summary(
                np.concatenate(tuple(per_tissue.values()))
            ),
        }

    replicate_variance = {
        "rna_log_cpm": _replicate_variance_summary(
            rna_log_cpm, sample_state_indices, design["n_states"]
        ),
        "protein_log_intensity": _replicate_variance_summary(
            protein_samples, sample_state_indices, design["n_states"]
        ),
    }
    component_variances = {
        "shared": {
            name: float(np.mean(component**2))
            for name, component in shared["scaled_components"].items()
        },
        "rna": {
            name: float(np.mean(component**2))
            for name, component in rna_biology["scaled_components"].items()
        },
        "protein": {
            name: float(np.mean(component**2))
            for name, component in protein["protein_scaled_components"].items()
        },
    }
    rna_distances = distances(rna_truth)
    diagnostics = {
        "component_variances": component_variances,
        "state_variance_fractions": state_variance_fractions,
        "total_dynamic_variance": float(np.var(shared["shared_dynamic"])),
        "rna_total_dynamic_variance": float(
            np.var(rna_truth - shared["shared_baseline"][:, None])
        ),
        "protein_total_dynamic_variance": float(
            np.var(protein_truth - protein_baseline[:, None])
        ),
        "mean_between_tissue_distance": rna_distances["between_tissue"],
        "mean_within_tissue_condition_distance": rna_distances[
            "within_tissue_condition"
        ],
        "modality_distances": {
            "shared": distances(shared_truth),
            "rna": rna_distances,
            "protein": distances(protein_truth),
        },
        "rna_protein_correlations": correlation_diagnostics,
        "replicate_variance": replicate_variance,
        "coupling_requested_mean": float(np.mean(protein["coupling_values"])),
        "modality_response_noise": float(parameters["modality_response_noise"]),
        "replicate_activity_noise": float(parameters["replicate_activity_noise"]),
        "replicate_tissue_activity_noise": float(
            parameters["replicate_tissue_activity_noise"]
        ),
        "baseline_coupling_requested": float(parameters["baseline_coupling"]),
        "baseline_coupling_realized": baseline_correlation,
        "missing_fraction": float(np.mean(observed["missingness_mask"])),
        "mean_library_depth": float(np.mean(replicates["library_depths"])),
    }
    return diagnostics

def _assemble_result(
    design: dict[str, Any],
    parameters: dict[str, Any],
    assumptions: _ScientificAssumptions,
    shared: dict[str, Any],
    rna_biology: dict[str, Any],
    rna: dict[str, Any],
    protein: dict[str, Any],
    replicates: dict[str, Any],
    observed: dict[str, Any],
    diagnostics: dict[str, Any],
    group_data: dict[str, Any] | None,
) -> dict[str, Any]:
    """Build the public result dictionary and its reconstruction truth."""
    rna_program_loadings = rna_biology["program_loadings"] + shared["program_loadings"]
    rna_program_activities = rna_biology["program_activities"]
    group_truth = None
    if group_data is not None:
        group_truth = {
            "group_names": group_data["group_names"],
            "membership": group_data["membership"],
            "coefficients": group_data["coefficients"],
            "effective_loadings": {
                "shared": shared["group_effective_loadings"],
                "rna": rna_biology["group_effective_loadings"],
                "protein": protein["protein_group_effective_loadings"],
            },
            "program_names": design["program_names"],
        }
    truth = {
        "shared": shared["shared"],
        "shared_baseline": shared["shared_baseline"],
        "shared_dynamic": shared["shared_dynamic"],
        "rna": rna["rna"],
        "rna_specific": rna_biology["dynamic"],
        "rna_state_noise": rna["rna_state_noise"],
        "protein": protein["protein"],
        "protein_baseline": protein["protein_baseline"],
        "protein_rna_component": protein["protein_rna_component"],
        "protein_specific": protein["protein_specific"],
        "program_names": design["program_names"],
        "program_activities": rna_program_activities,
        "program_loadings": rna_program_loadings,
        "module_membership": rna_biology["module_membership"] | shared["module_membership"],
        "component_programs": _COMPONENT_PROGRAMS.copy(),
        "scaled_components": shared["scaled_components"],
        "rna_scaled_components": rna_biology["scaled_components"],
        "rna_module_membership": rna_biology["module_membership"],
        "protein_program_activities": protein["protein_program_activities"],
        "protein_program_loadings": protein["protein_program_loadings"],
        "protein_program_names": protein["protein_program_names"],
        "protein_specific_variance_weights": protein["protein_specific_variance_weights"],
        "protein_scaled_components": protein["protein_scaled_components"],
        "protein_module_membership": protein["protein_module_membership"],
        "groups": group_truth,
        "sample_activity_perturbations": {
            "shared_tissue": replicates["shared_activity_perturbations"],
            "rna": replicates["rna_activity_perturbations"],
            "protein": replicates["protein_activity_perturbations"],
        },
        "contrasts": _make_contrasts(
            shared["shared"], design["tissues"], design["conditions"], design["n_times"]
        ),
        "rna_contrasts": _make_contrasts(
            rna["rna"], design["tissues"], design["conditions"], design["n_times"]
        ),
        "protein_contrasts": _make_contrasts(
            protein["protein"], design["tissues"], design["conditions"], design["n_times"]
        ),
    }
    public_parameters = {
        "seed": int(parameters["seed"]),
        "gene_count": design["gene_count"],
        "gene_ids": design["gene_ids"],
        "tissues": design["tissues"],
        "conditions": design["conditions"],
        "timepoints": tuple(map(float, design["timepoints"])),
        "replicates": design["n_replicates"],
        "variance_weights": parameters["variance_weights"],
        "effect_scale": parameters["effect_scale"],
        "module_probability": float(parameters["module_probability"]),
        "main_time_secondary_loading_scale": parameters[
            "main_time_secondary_loading_scale"
        ],
        "coupling_mean": float(parameters["coupling_mean"]),
        "coupling_sd": parameters["coupling_sd"],
        "coupling_range": (parameters["coupling_low"], parameters["coupling_high"]),
        "baseline_coupling": float(parameters["baseline_coupling"]),
        "protein_baseline_sd": assumptions.protein_baseline_sd,
        "biological_noise": parameters["biological_noise"],
        "replicate_noise": parameters["replicate_noise"],
        "modality_response_noise": parameters["modality_response_noise"],
        "replicate_activity_noise": parameters["replicate_activity_noise"],
        "replicate_tissue_activity_noise": parameters[
            "replicate_tissue_activity_noise"
        ],
        "library_depth_expected": parameters["library_depth_expected"],
        "library_depth_log_sd": assumptions.library_depth_log_sd,
        "dispersion_mean": parameters["dispersion"],
        "proteomics_noise": parameters["proteomics_noise"],
        "missingness": bool(parameters["missingness"]),
        "missingness_intercept": float(parameters["missingness_intercept"]),
        "missingness_slope": parameters["missingness_slope"],
        "gene_groups": parameters["gene_groups"],
        "group_program_weights": parameters["group_program_weights"],
        "group_strength": parameters["group_strength"],
    }
    return {
        "gene_ids": design["gene_ids"],
        "rna_counts": observed["rna_counts"],
        "protein_log_intensities": observed["protein_log_intensities"],
        "state_metadata": design["state_metadata"],
        "sample_metadata": replicates["sample_metadata"],
        "sample_state_indices": replicates["sample_state_indices"],
        "truth": truth,
        "library_depths": replicates["library_depths"],
        "gene_dispersions": observed["gene_dispersions"],
        "coupling_values": protein["coupling_values"],
        "missingness_mask": observed["missingness_mask"],
        "parameters": public_parameters,
        "diagnostics": diagnostics,
    }

def simulate(
    seed: int = 0,
    *,
    gene_count: int | None = None,
    gene_ids: Sequence[str] | None = None,
    tissues: Sequence[str] = ("liver", "muscle"),
    conditions: Sequence[str] = ("WT", "ob/ob"),
    timepoints: Sequence[float] | None = None,
    replicates: int = 3,
    variance_weights: Mapping[str, float] | Sequence[float] | None = None,
    effect_scale: float = 1.0,
    module_probability: float = 0.20,
    main_time_secondary_loading_scale: float = 0.25,
    gene_groups: Mapping[str, Sequence[str]] | None = None,
    group_program_weights: Mapping[str, Mapping[str, float]] | None = None,
    group_strength: float = 0.8,
    coupling_mean: float = 0.70,
    coupling_sd: float = 0.15,
    coupling_range: tuple[float, float] = (0.15, 0.95),
    baseline_coupling: float = 0.70,
    biological_noise: float = 0.15,
    replicate_noise: float = 0.15,
    modality_response_noise: float = 0.35,
    replicate_activity_noise: float = 0.15,
    replicate_tissue_activity_noise: float = 0.0,
    library_depth: float | Sequence[float] = 1_000_000,
    dispersion: float = 0.10,
    proteomics_noise: float = 0.20,
    missingness: bool = False,
    missingness_intercept: float = -2.5,
    missingness_slope: float = 1.0,
) -> dict[str, Any]:
    """Generate paired RNA-seq counts and proteomics log-intensities.

    Parameters are deliberately ordinary Python values. ``gene_ids`` can be
    supplied instead of ``gene_count``; otherwise identifiers are generated as
    ``gene_0000``, etc. A mapping passed as ``variance_weights`` may override
    any subset of the seven defaults, but the final weights must sum to one.

    ``biological_noise`` controls RNA-specific state variation, while
    ``replicate_noise`` controls independent per-gene deviations for both
    modalities. ``modality_response_noise`` controls how much each modality's
    smooth temporal bases vary independently. ``replicate_activity_noise``
    perturbs the non-tissue program activities independently by modality, and
    ``replicate_tissue_activity_noise`` separately controls a shared replicate
    shift to tissue identity. ``library_depth`` is an expected total and may be a
    scalar or one expected total per sample. ``dispersion`` is the mean gene dispersion
    in ``Var(count) = mean + dispersion * mean**2``; individual gene values
    vary modestly around it. Shared and modality-specific state effects scale
    together, so ``effect_scale`` changes both modalities.
    ``baseline_coupling`` blends standardized shared-RNA gene baselines with
    independent protein-specific baselines; zero and one select independent
    and fully shared baseline ranks, respectively. ``coupling_mean``,
    ``coupling_sd``, and ``coupling_range`` control gene-specific RNA mixing
    weights in protein's non-tissue dynamics. The tissue program is shared;
    all other program loadings and temporal responses are independently
    generated per modality. The mixing weight is not a correlation target.
    ``diagnostics['rna_protein_correlations']`` reports latent and observed
    per-gene correlations, separating within-tissue genotype/time-state
    correlations from pooled-tissue values.
    ``main_time_secondary_loading_scale`` scales non-dominant active loadings
    among ``time_linear``, ``time_quadratic``, and ``time_bump`` before their
    existing variance normalization. Its default is 0.25; set it to 1.0 to
    leave overlapping main-time loadings unscaled.
    ``gene_groups`` maps nonempty group names to gene IDs and may contain
    overlapping groups. ``group_program_weights`` supplies explicit nominal
    coefficients for selected groups; omitted programs in an explicit profile
    receive zero, while groups without a profile receive random coefficients
    from the simulator's existing sparsity and loading distributions. Explicit
    non-tissue coefficients are shared between RNA and protein before their
    existing modality-specific activities and protein coupling are applied;
    random non-tissue coefficients are generated independently for each
    modality. Tissue group coefficients are shared.

    For a gene in ``k`` groups, the membership basis gives each membership a
    weight of ``1/sqrt(k)``. If ``M`` is this normalized gene-by-group matrix
    and ``W`` is the group-by-program coefficient matrix, group loadings are
    ``G = M @ W``. For grouped rows, each component blends its existing
    individual loading ``I`` with ``G`` as
    ``sqrt(1 - group_strength) * I + sqrt(group_strength) * G``. Ungrouped
    rows keep their individual loadings. The simulator then applies its usual
    per-component RMS normalization to the combined signal and the reported
    group contribution. Thus ``group_strength`` controls the blend before
    normalization; it is not an exact fraction of realized variance, since
    overlapping profiles may correlate and normalization depends on the
    resulting signal. A component with no raw signal cannot realize a positive
    variance budget and retains the existing zero-signal behavior.

    ``PROGRAM_NAMES`` gives the stable program-column order. When groups are
    enabled, ``truth["groups"]`` contains ``group_names`` with shape
    ``(n_groups,)``; boolean ``membership`` with shape
    ``(n_genes, n_groups)`` in ``gene_ids`` row order; ``program_names`` equal
    to ``PROGRAM_NAMES``; ``coefficients`` with shared, RNA, and protein
    matrices of shape ``(n_groups, n_programs)`` in that program order; and
    ``effective_loadings`` with shared, RNA, and protein matrices of shape
    ``(n_genes, n_programs)``. Effective loadings are the normalized group
    contributions. Protein effective loadings include the existing per-gene
    residual coupling factor. ``truth["groups"]`` is ``None`` when no groups
    are supplied.

    Returns a dictionary of observed matrices, metadata, latent truths,
    parameters, and compact diagnostics. All matrices have genes on rows;
    ``missingness_mask`` is true where a protein intensity is missing.
    Contrast signs are second label minus first label.
    """
    rng = np.random.default_rng(seed)
    assumptions = _ScientificAssumptions()
    design = _make_design(
        gene_count=gene_count,
        gene_ids=gene_ids,
        tissues=tissues,
        conditions=conditions,
        timepoints=timepoints,
        replicates=replicates,
        assumptions=assumptions,
    )
    group_options = _validate_gene_group_options(
        gene_groups, group_program_weights, group_strength, design["gene_ids"]
    )
    parameters = _validate_parameters(
        design=design,
        seed=seed,
        variance_weights=variance_weights,
        effect_scale=effect_scale,
        module_probability=module_probability,
        main_time_secondary_loading_scale=main_time_secondary_loading_scale,
        coupling_mean=coupling_mean,
        coupling_sd=coupling_sd,
        coupling_range=coupling_range,
        baseline_coupling=baseline_coupling,
        biological_noise=biological_noise,
        replicate_noise=replicate_noise,
        modality_response_noise=modality_response_noise,
        replicate_activity_noise=replicate_activity_noise,
        replicate_tissue_activity_noise=replicate_tissue_activity_noise,
        library_depth=library_depth,
        dispersion=dispersion,
        proteomics_noise=proteomics_noise,
        missingness=missingness,
        missingness_intercept=missingness_intercept,
        missingness_slope=missingness_slope,
    )
    parameters.update(group_options)

    group_data = None
    if group_options["gene_groups"]:
        # This dedicated branch does not consume or shift the individual draws.
        group_rng = np.random.default_rng(
            np.random.SeedSequence(seed, spawn_key=(0x475250,))
        )
        group_data = _make_group_profiles(
            group_rng, design, parameters, assumptions, group_options
        )

    # Shared tissue scaffold followed by independent non-tissue programs.
    shared = _generate_shared_biology(
        rng, design, parameters, assumptions, group_data
    )

    rna_biology = _generate_modality_biology(
        rng, design, parameters, assumptions, group_data, "rna"
    )
    rna = _generate_rna(rng, shared, rna_biology, parameters["biological_noise"])
    protein_biology = _generate_modality_biology(
        rng, design, parameters, assumptions, group_data, "protein"
    )
    protein = _generate_protein(
        rng, rna, rna_biology, shared, protein_biology, design, parameters, assumptions
    )

    # Replicate-level sample values followed by assay observations.
    replicates_data = _generate_replicates(
        rng,
        rna["rna"],
        rna_biology,
        protein["protein"],
        protein,
        shared,
        design,
        parameters,
        assumptions,
    )
    observed = _generate_observations(
        rng, replicates_data, design, parameters, assumptions
    )
    diagnostics = _calculate_diagnostics(
        shared, rna, rna_biology, protein, replicates_data, observed, design, parameters
    )
    return _assemble_result(
        design,
        parameters,
        assumptions,
        shared,
        rna_biology,
        rna,
        protein,
        replicates_data,
        observed,
        diagnostics,
        group_data,
    )

def synthetic_longitudinal(
    seed: int = 0,
    *,
    gene_count: int | None = None,
    gene_ids: Sequence[str] | None = None,
    timepoints: int | Sequence[float] = 8,
    tissues: Sequence[str] = ("liver", "muscle"),
    conditions: Sequence[str] = ("WT", "ob/ob"),
    replicates: int = 3,
    gene_groups: Mapping[str, Sequence[str]] | None = None,
    group_program_weights: Mapping[str, Mapping[str, float]] | None = None,
    group_strength: float = 0.8,
    **simulation_options: Any,
) -> dict[str, Any]:
    """Generate paired sample- and state-level longitudinal omics matrices.

    An integer ``timepoints`` gives the number of equally spaced simulation
    units between zero and one; a sequence supplies explicit, ordered numeric
    timepoints. The default design has liver and muscle, WT and ob/ob, eight
    timepoints, and three replicate samples per state. Tissue, condition, and
    timepoint order in the outputs follows the supplied order.

    Sample matrices have samples on rows and genes on columns. The two state
    matrices have one ordered tissue-condition-time row per state (32 by
    default) and genes on columns.
    ``state_rna_log_cpm`` is the mean of replicate log1p(CPM) values, while
    ``state_protein_log_intensities`` is the mean observed intensity. These are
    convenient summaries; use the replicate-level matrices for statistical
    fitting so replicate variability is not discarded.

    ``truth`` contains state truth matrices with states on rows and genes on
    columns. Its ``shared_contrasts``, ``rna_contrasts``, and
    ``protein_contrasts`` are dictionaries whose numeric arrays have genes on
    the first axis: ``condition_within_tissue_time`` and
    ``tissue_within_condition_time`` have shape (genes, 2, timepoints),
    ``adjacent_time`` has shape (genes, 2, 2, timepoints - 1), and
    ``condition_by_time_from_first`` has shape (genes, 2, timepoints).
    When groups are supplied, ``truth["groups"]`` carries ``group_names``
    (``n_groups``), a gene-first boolean ``membership`` matrix
    (``n_genes``, ``n_groups``), and gene-first effective loading matrices
    (``n_genes``, ``n_programs``) for shared, RNA, and protein contributions.
    Its ``program_names`` are the public ``PROGRAM_NAMES`` order. The
    ``coefficients`` matrices for shared, RNA, and protein use groups on rows
    and programs on columns (``n_groups``, ``n_programs``); their column order
    is also ``PROGRAM_NAMES``. Protein effective loadings include the existing
    protein residual coupling factor.

    Pass GEM-mapped ``gene_ids`` for downstream use. Other simulator options,
    such as noise parameters, effect scale, or variance weights, can be passed
    as keyword arguments.
    """
    if isinstance(timepoints, (int, np.integer)) and not isinstance(
        timepoints, (bool, np.bool_)
    ):
        if timepoints < 3:
            raise ValueError("timepoints count must be at least three")
        time_values = np.linspace(0.0, 1.0, int(timepoints))
    elif np.isscalar(timepoints):
        raise ValueError("timepoints must be an integer count or a numeric sequence")
    else:
        time_values = timepoints

    simulation_options.setdefault("variance_weights", _DEFAULT_VARIANCE_WEIGHTS)
    generated = simulate(
        seed=seed,
        gene_count=gene_count,
        gene_ids=gene_ids,
        tissues=tissues,
        conditions=conditions,
        timepoints=time_values,
        replicates=replicates,
        gene_groups=gene_groups,
        group_program_weights=group_program_weights,
        group_strength=group_strength,
        **simulation_options,
    )

    counts = generated["rna_counts"]
    sample_totals = np.maximum(counts.sum(axis=0, keepdims=True), 1)
    sample_log_cpm = np.log1p(counts / sample_totals * 1e6)
    state_indices = generated["sample_state_indices"]
    n_states = len(generated["state_metadata"])

    return {
        "gene_ids": generated["gene_ids"],
        "state_metadata": generated["state_metadata"],
        "sample_metadata": generated["sample_metadata"],
        "state_rna_log_cpm": _state_means(sample_log_cpm, state_indices, n_states).T,
        "state_protein_log_intensities": _state_means(
            generated["protein_log_intensities"], state_indices, n_states
        ).T,
        "sample_rna_counts": counts.T,
        "sample_protein_log_intensities": generated["protein_log_intensities"].T,
        "protein_missingness_mask": generated["missingness_mask"].T,
        "library_depths": generated["library_depths"],
        "gene_dispersions": generated["gene_dispersions"],
        "truth": {
            "state_rna_log_abundance": generated["truth"]["rna"].T,
            "state_protein_log_intensities": generated["truth"]["protein"].T,
            "shared_contrasts": generated["truth"]["contrasts"],
            "rna_contrasts": generated["truth"]["rna_contrasts"],
            "protein_contrasts": generated["truth"]["protein_contrasts"],
            "groups": generated["truth"]["groups"],
        },
        "parameters": generated["parameters"],
        "diagnostics": generated["diagnostics"],
    }


__all__ = ["PROGRAM_NAMES", "simulate", "synthetic_longitudinal"]
