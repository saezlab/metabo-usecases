"""Temporal regularization for CORNETO's multi-sample iMAT formulation.

The implementation composes :class:`corneto.methods.imat.MultiSampleIMAT`
problems. It first measures each sample's independently achievable iMAT
agreement, then solves a joint problem constrained to retain a requested
fraction of that agreement while penalizing the reaction union and changes in
reaction selection between explicitly connected samples.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import time
from typing import Any

import numpy as np
from corneto.backend import CvxpyBackend, VarType
from corneto.methods.imat import MultiSampleIMAT
from scipy.sparse import coo_matrix

from corneto_mcmt.activity import ActivityClassifier, build_activity_classifier


_FEASIBILITY_TOL = 1e-6
_INTEGRALITY_TOL = 1e-5
_DEFAULT_GUROBI_INT_FEAS_TOL = 1e-9


@dataclass(frozen=True)
class StateFitDiagnostics:
    """Independent and joint iMAT fit measurements for one condition."""

    evidence_weight: float
    independent_mismatch: float
    independent_agreement: float
    fit_floor: float
    joint_mismatch: float
    joint_agreement: float
    fit_margin: float
    has_evidence: bool
    reference_status: str


@dataclass(frozen=True)
class WarmStartDiagnostics:
    """Validation summary for an assembled joint MIP start."""

    decision_variables: int
    initialized_variables: int
    decision_symbols: int
    objective_value: float
    flux_penalty: str
    flux_l1: float
    flux_l2_squared: float
    flux_objective: float
    canonical_objective_kind: str
    canonical_quadratic_objective_nonzeros: int
    maximum_constraint_violation: float
    maximum_bound_violation: float
    cvxpy_variables: int
    canonical_variables: int
    canonical_auxiliary_variables: int
    uninitialized_canonical_variables: int
    repaired_reference_thresholds: int
    independent_reference_mismatch: Mapping[str, float]
    reference_mismatch_before_repair: Mapping[str, float]
    reference_mismatch_after_repair: Mapping[str, float]
    fit_floor_margins: Mapping[str, float]
    reference_repair_statuses: Mapping[str, str]
    reference_repair_seconds: Mapping[str, float]


class _WarmStartFitFloorError(RuntimeError):
    """The repaired reference assignment cannot seed the original fit floors."""


@dataclass(frozen=True)
class _ReferenceSolution:
    """Detached numeric values needed to compose a joint feasible start."""

    flux: np.ndarray
    positive_thresholds: Mapping[str, tuple[float, float]]
    mismatch: float
    repaired_thresholds: int


@dataclass(frozen=True)
class TemporalIMATResult:
    """Solved temporal iMAT model and labeled numerical diagnostics.

    ``flux`` and ``selected`` have shape ``(n_reactions, n_conditions)`` even
    for a single condition. ``selected`` records CORNETO's full
    ``edge_has_flux`` network-selection indicators; a selected reaction may
    carry zero flux when it is unscored. ``activity`` is the threshold-based
    reaction-by-condition matrix (``|flux| >= eps``; MILP-classified in
    ``"activity"`` mode, post-solve thresholded in ``"flux"`` mode, and
    ``None`` in legacy ``"selection"`` mode). Switch diagnostics describe
    ``activity`` when present and ``selected`` otherwise.
    """

    problem: Any
    solver_result: Any
    solver_status: str
    optimality_proven: bool
    reactions: tuple[str, ...]
    conditions: tuple[str, ...]
    temporal_pairs: tuple[tuple[str, str], ...]
    flux: np.ndarray
    selected: np.ndarray
    activity: np.ndarray | None
    temporal_mode: str
    union: np.ndarray
    state_fits: Mapping[str, StateFitDiagnostics]
    switching_counts: Mapping[tuple[str, str], int]
    total_switches: int
    union_size: int
    union_objective: float
    temporal_objective: float
    lambda_flux: float
    flux_penalty: str
    flux_l1: float
    flux_l2_squared: float
    flux_objective: float
    objective_value: float
    warm_start: WarmStartDiagnostics | None = None


class TemporalMultiSampleIMAT:
    """Fit a multi-condition iMAT network with explicit temporal smoothing.

    ``lambda_reg`` weights CORNETO's existing global reaction-union penalty.
    ``lambda_temporal`` weights changes in threshold activity by default,
    structural ``edge_has_flux`` selection when ``temporal_mode="selection"``,
    or continuous reaction flux differences when ``temporal_mode="flux"``.
    ``lambda_flux`` optionally weights a raw flux-magnitude penalty over all
    reaction fluxes in all conditions. ``flux_penalty="l1"`` (the default)
    uses the L1 norm and adds a continuous epigraph variable per
    reaction/condition; ``flux_penalty="l2"`` uses the squared L2 norm as a
    native convex quadratic objective and adds no variables. The weight's
    effect depends on the flux units used by the model.
    Every condition first receives an independent
    reference fit, and the joint solve must retain at least
    ``1 - epsilon_fit`` of its independently achievable agreement.

    Threshold activity uses the same ``eps`` as iMAT. Flux strictly above the
    threshold is active, flux strictly inside it is inactive, and the closed
    MILP boundary permits either label. ``result.selected`` always reports
    structural selectors; ``result.activity`` reports threshold activity in
    the default ``"activity"`` mode and in ``"flux"`` mode (computed
    post-solve as ``|flux| >= eps - tol`` without MILP classifier variables),
    and is ``None`` in ``"selection"`` mode.
    Switch counts describe ``result.activity`` when present and
    ``result.selected`` otherwise; ``result.temporal_objective`` reports the
    weighted total variation for the expression selected by ``temporal_mode``.
    Only caller-supplied condition pairs are connected.

    Parameters mirror CORNETO's iMAT settings where relevant. Reaction scores
    are the only evidence input; metabolic functions belong in reaction
    bounds. The solver is selected at :meth:`fit` time and defaults to
    SciPy's MILP interface for linear objectives. A positive squared-L2 weight
    requires an MIQP-capable solver such as GUROBI.
    """

    def __init__(
        self,
        lambda_reg: float = 1.0,
        lambda_temporal: float = 0.1,
        epsilon_fit: float = 0.05,
        eps: float = 1e-2,
        scale: bool = True,
        temporal_mode: str = "activity",
        lambda_flux: float = 0.0,
        flux_penalty: str = "l1",
    ) -> None:
        self.lambda_reg = _finite_nonnegative("lambda_reg", lambda_reg)
        self.lambda_temporal = _finite_nonnegative("lambda_temporal", lambda_temporal)
        self.lambda_flux = _finite_nonnegative("lambda_flux", lambda_flux)
        if not isinstance(flux_penalty, str) or flux_penalty not in {"l1", "l2"}:
            raise ValueError("flux_penalty must be either 'l1' or 'l2'")
        self.flux_penalty = flux_penalty
        self.epsilon_fit = float(epsilon_fit)
        if not np.isfinite(self.epsilon_fit) or not 0 <= self.epsilon_fit <= 1:
            raise ValueError("epsilon_fit must be finite and between 0 and 1")
        self.eps = float(eps)
        if not np.isfinite(self.eps) or self.eps <= 0:
            raise ValueError("eps must be a finite positive number")
        if not isinstance(scale, (bool, np.bool_)):
            raise ValueError("scale must be a boolean")
        self.scale = bool(scale)
        if not isinstance(temporal_mode, str) or temporal_mode not in {
            "activity",
            "selection",
            "flux",
        }:
            raise ValueError("temporal_mode must be 'activity', 'selection', or 'flux'")
        self.temporal_mode = temporal_mode
        self.backend = CvxpyBackend()
        self.problem = None
        self.result: TemporalIMATResult | None = None

    def fit(
        self,
        graph: Any,
        *,
        reaction_scores: Mapping[str, Mapping[str, float]],
        reaction_bounds: Mapping[str, Mapping[str, tuple[float, float]]] | None = None,
        temporal_pairs: Sequence[tuple[str, str]] = (),
        solver: str = "SCIPY",
        max_seconds: int | None = None,
        verbosity: int = 0,
        solver_options: Mapping[str, Any] | None = None,
        joint_max_seconds: int | None = None,
        joint_solver_options: Mapping[str, Any] | None = None,
        accept_feasible_joint_incumbent: bool = False,
        warm_start: bool = False,
        reference_repair_max_seconds: int = 300,
    ) -> TemporalIMATResult:
        """Build, solve, and validate independent and joint temporal iMAT fits.

        ``solver_options`` contains backend-specific options forwarded to every
        independent reference solve and the joint solve. Gurobi solves default
        to ``IntFeasTol=1e-9``; an explicit ``IntFeasTol`` entry in
        ``solver_options`` (or ``joint_solver_options`` for the joint solve)
        overrides that default, case-insensitively.
        ``joint_max_seconds`` overrides the shared time limit for only the joint
        solve, and ``joint_solver_options`` adds or overrides backend-specific
        options for only that solve. A time-limited joint incumbent is rejected
        by default; set ``accept_feasible_joint_incumbent`` to accept a verified
        Gurobi time-limit incumbent without claiming optimality.
        With ``warm_start=True``, detached numeric reference fits are composed
        into a validated joint start and passed through CVXPY's GUROBI interface.
        When exact threshold indicators make a reference violate its unchanged
        fit floor, that reference is re-solved with tighter Gurobi integrality
        tolerances and the same floor before the joint start is composed. Each
        repair uses ``reference_repair_max_seconds`` as its time limit. If no
        validated repair is found, fitting stops before the joint solve. When
        flux regularization or ``temporal_mode="flux"`` is active, each
        warm-start reference flux is also regularized under its exactized iMAT
        score state to eliminate futile cycles before seeding the joint solve.
        ``temporal_mode="activity"`` smooths threshold activity (the default);
        ``temporal_mode="selection"`` retains structural-selector smoothing;
        ``temporal_mode="flux"`` smooths continuous reaction fluxes without
        introducing binary activity-classifier variables.
        ``lambda_flux`` adds a magnitude penalty over every reaction and
        condition to the joint objective only: L1 for ``flux_penalty="l1"``
        and squared L2 for ``flux_penalty="l2"``. The native reference iMAT
        objective and mismatch-to-agreement floors are unchanged. A positive
        L1 weight adds a continuous absolute-value epigraph; L2 is a native
        convex quadratic objective with linear constraints. When
        ``lambda_flux`` is zero, neither formulation adds variables or
        constraints.
        """

        if solver_options is not None and not isinstance(solver_options, Mapping):
            raise TypeError("solver_options must be a mapping of backend-specific options")
        if joint_solver_options is not None and not isinstance(joint_solver_options, Mapping):
            raise TypeError("joint_solver_options must be a mapping of backend-specific options")
        if not isinstance(accept_feasible_joint_incumbent, (bool, np.bool_)):
            raise ValueError("accept_feasible_joint_incumbent must be a boolean")
        if not isinstance(warm_start, (bool, np.bool_)):
            raise ValueError("warm_start must be a boolean")
        if (
            not isinstance(reference_repair_max_seconds, (int, np.integer))
            or reference_repair_max_seconds <= 0
        ):
            raise ValueError("reference_repair_max_seconds must be a positive integer")
        solver_options = dict(solver_options or {})
        joint_solver_options = _merge_solver_options(
            solver_options, dict(joint_solver_options or {})
        )
        joint_time_limit = max_seconds if joint_max_seconds is None else joint_max_seconds

        conditions, scores, bounds = _normalize_inputs(reaction_scores, reaction_bounds)
        pairs = _normalize_pairs(temporal_pairs, conditions)
        reactions = _reaction_labels(graph)

        available = {name.upper() for name in self.backend.available_solvers()}
        if not isinstance(solver, str) or solver.upper() not in available:
            raise ValueError(
                f"Solver {solver!r} is unavailable through CORNETO's CVXPY backend; "
                f"available solvers are {sorted(available)}."
            )
        if warm_start and solver.upper() != "GUROBI":
            raise ValueError(
                "warm_start=True currently requires CVXPY's GUROBI interface, which "
                "accepts the assembled initial-value vector."
            )

        evidence_weights: dict[str, float] = {}
        independent_mismatch: dict[str, float] = {}
        reference_status: dict[str, str] = {}
        reference_solutions: dict[str, _ReferenceSolution] = {}
        reference_mismatch_before_repair: dict[str, float] = {}
        reference_repair_statuses: dict[str, str] = {}
        reference_repair_seconds: dict[str, float] = {}
        repaired_reference_thresholds = 0

        for condition in conditions:
            reference_method = _imat_method(
                self.backend,
                lambda_reg=0.0,
                eps=self.eps,
                scale=self.scale,
            )
            reference = reference_method.build_many(
                graph,
                reaction_scores={condition: scores[condition]},
                reaction_bounds={condition: bounds[condition]},
            )
            evidence_weights[condition] = _effective_evidence_weight(
                scores[condition], scale=self.scale
            )
            reference_fit_expr = _state_mismatch_expression(reference, condition, 0)
            status, result = _solve_checked(
                reference,
                solver,
                max_seconds=max_seconds,
                verbosity=verbosity,
                solver_options=solver_options,
                stage=f"independent reference for condition {condition!r}",
            )
            reference_status[condition] = status
            independent_mismatch[condition] = _expression_value(reference_fit_expr)
            _validate_mismatch(
                independent_mismatch[condition],
                evidence_weights[condition],
                stage=f"independent reference for condition {condition!r}",
            )
            if warm_start:
                initial_reference_solution = _capture_reference_solution(
                    reference,
                    processed_data=reference_method.processed_data,
                    condition=condition,
                    reactions=reactions,
                    scale=self.scale,
                    eps=self.eps,
                )
                repaired_reference_thresholds += initial_reference_solution.repaired_thresholds
                reference_mismatch_before_repair[condition] = (
                    initial_reference_solution.mismatch
                )
                evidence_weight = evidence_weights[condition]
                original_mismatch = independent_mismatch[condition]
                maximum_mismatch = original_mismatch + self.epsilon_fit * max(
                    0.0, evidence_weight - original_mismatch
                )
                if (
                    evidence_weight > 0
                    and initial_reference_solution.mismatch
                    > maximum_mismatch + _FEASIBILITY_TOL * max(1.0, evidence_weight)
                ):
                    repaired_reference_solution, repair_status, repair_seconds = (
                        _repair_reference_solution(
                            graph,
                            scores=scores[condition],
                            bounds=bounds[condition],
                            condition=condition,
                            original_mismatch=original_mismatch,
                            evidence_weight=evidence_weight,
                            epsilon_fit=self.epsilon_fit,
                            eps=self.eps,
                            scale=self.scale,
                            solver=solver,
                            solver_options=solver_options,
                            max_seconds=int(reference_repair_max_seconds),
                            verbosity=verbosity,
                        )
                    )
                    if repaired_reference_solution.mismatch > (
                        maximum_mismatch
                        + _FEASIBILITY_TOL * max(1.0, evidence_weight)
                    ):
                        raise RuntimeError(
                            f"The repaired independent reference for {condition!r} still "
                            f"violates the unchanged agreement floor: original objective "
                            f"mismatch={original_mismatch:.8g}, allowed mismatch="
                            f"{maximum_mismatch:.8g}, repaired mismatch="
                            f"{repaired_reference_solution.mismatch:.8g}. The joint solve "
                            "was not started."
                        )
                    reference_solutions[condition] = repaired_reference_solution
                    reference_repair_statuses[condition] = repair_status
                    reference_repair_seconds[condition] = repair_seconds
                else:
                    reference_solutions[condition] = initial_reference_solution
                    reference_repair_statuses[condition] = "not_needed"
                    reference_repair_seconds[condition] = 0.0
                if self.lambda_flux > 0 or self.temporal_mode == "flux":
                    reference_solutions[condition] = _regularize_reference_solution(
                        graph,
                        reference_problem=reference,
                        reference_solution=reference_solutions[condition],
                        processed_data=reference_method.processed_data,
                        condition=condition,
                        reactions=reactions,
                        scale=self.scale,
                        eps=self.eps,
                        lambda_flux=self.lambda_flux,
                        flux_penalty=self.flux_penalty,
                        solver=solver,
                        solver_options=solver_options,
                        max_seconds=max_seconds,
                        verbosity=verbosity,
                    )

        method = _imat_method(
            self.backend,
            lambda_reg=self.lambda_reg,
            eps=self.eps,
            scale=self.scale,
            force_structural_selectors=True,
        )
        problem = method.build_many(
            graph,
            reaction_scores=scores,
            reaction_bounds=bounds,
        )
        activity_classifier = None
        if self.temporal_mode == "activity":
            positive_reactions_by_condition = {
                condition: tuple(
                    reaction
                    for reaction, weight in _scored_reaction_weights(
                        method.processed_data, condition, scale=self.scale
                    )
                    if weight > 0
                )
                for condition in conditions
            }
            activity_classifier = build_activity_classifier(
                problem,
                reactions=reactions,
                conditions=conditions,
                positive_reactions_by_condition=positive_reactions_by_condition,
                eps=self.eps,
            )
        flux_abs_symbol = None
        if self.lambda_flux > 0:
            flow = problem.expr.flow
            if self.flux_penalty == "l1":
                flux_abs_symbol = self.backend.Variable(
                    name="flux_abs",
                    shape=flow.shape,
                    lb=0.0,
                    vartype=VarType.CONTINUOUS,
                )
                # A linear epigraph expresses |v| without introducing CVXPY's
                # hidden absolute-value auxiliaries into warm-start assembly.
                problem.add_constraints([flux_abs_symbol >= flow, flux_abs_symbol >= -flow])
                problem.add_objective(
                    self.lambda_flux * flux_abs_symbol.sum(),
                    name="flux_l1_penalty",
                )
            else:
                # CORNETO's CvxpyExpression power operator preserves the sum
                # of squares as a convex quadratic objective for GUROBI's QP
                # interface; this adds no norm cone or auxiliary variables.
                problem.add_objective(
                    self.lambda_flux * (flow**2).sum(),
                    name="flux_l2_squared_penalty",
                )
        joint_fit_expr: dict[str, Any] = {}
        for index, condition in enumerate(conditions):
            expression = _state_mismatch_expression(problem, condition, index)
            joint_fit_expr[condition] = expression
            weight = evidence_weights[condition]
            if weight <= 0:
                continue
            optimum = independent_mismatch[condition]
            agreement = max(0.0, weight - optimum)
            max_mismatch = optimum + self.epsilon_fit * agreement
            problem.add_constraints(expression <= max_mismatch)

        # Keep CORNETO's named fit expressions and constraints, but remove
        # their optimization weight. Match exact generated names and condition
        # indices so regularization objectives remain untouched.
        objective_names = {objective.name: i for i, objective in enumerate(problem.objectives)}
        for index, condition in enumerate(conditions):
            for kind in ("pos", "neg"):
                name = _fit_objective_name(kind, condition, index)
                objective_index = objective_names.get(name)
                if objective_index is not None:
                    problem.weights[objective_index] = 0.0

        if self.lambda_temporal > 0 and pairs:
            if self.temporal_mode == "activity":
                temporal_expression = activity_classifier.expression
                temporal_name = "activity"
            elif self.temporal_mode == "selection":
                temporal_expression = problem.expr.edge_has_flux
                temporal_name = "selector"
            else:
                temporal_expression = problem.expr.flow
                temporal_name = "flux"
            tv = self.backend.TotalVariation(
                temporal_expression,
                pairs=[(conditions.index(left), conditions.index(right)) for left, right in pairs],
                axis=1,
                name=f"temporal_{temporal_name}",
            )
            problem += tv
            problem.add_objective(
                self.lambda_temporal * tv.expr[f"temporal_{temporal_name}_total_variation"],
                name=f"temporal_{temporal_name}_switch_penalty",
            )

        warm_start_diagnostics = None
        if warm_start:
            warm_start_diagnostics = _assemble_and_validate_warm_start(
                problem,
                graph=graph,
                solver=solver,
                references=reference_solutions,
                conditions=conditions,
                reaction_scores=scores,
                processed_data=method.processed_data,
                independent_mismatch=independent_mismatch,
                evidence_weights=evidence_weights,
                reactions=reactions,
                temporal_pairs=pairs,
                scale=self.scale,
                lambda_reg=self.lambda_reg,
                lambda_temporal=self.lambda_temporal,
                lambda_flux=self.lambda_flux,
                flux_penalty=self.flux_penalty,
                epsilon_fit=self.epsilon_fit,
                reference_mismatch_before_repair=reference_mismatch_before_repair,
                reference_repair_statuses=reference_repair_statuses,
                reference_repair_seconds=reference_repair_seconds,
                repaired_reference_thresholds=repaired_reference_thresholds,
                temporal_mode=self.temporal_mode,
                activity_classifier=activity_classifier,
            )

        status, joint_result = _solve_checked(
            problem,
            solver,
            max_seconds=joint_time_limit,
            verbosity=verbosity,
            solver_options=joint_solver_options,
            allow_feasible_incumbent=bool(accept_feasible_joint_incumbent),
            warm_start=bool(warm_start and warm_start_diagnostics is not None),
            stage="joint temporal fit",
        )
        self.problem = problem

        flux = _as_reaction_condition_matrix(problem.expr.flow.value, len(reactions), len(conditions))
        selector_values = _as_reaction_condition_matrix(
            problem.expr.edge_has_flux.value, len(reactions), len(conditions)
        )
        if not np.all(np.isfinite(flux)) or not np.all(np.isfinite(selector_values)):
            raise RuntimeError("The optimal joint solve did not return finite fluxes and selectors.")
        if np.max(np.abs(selector_values - np.rint(selector_values)), initial=0.0) > _INTEGRALITY_TOL:
            raise RuntimeError("The joint solution violates binary selector integrality tolerance.")
        if np.min(selector_values, initial=0.0) < -_INTEGRALITY_TOL or np.max(
            selector_values, initial=0.0
        ) > 1.0 + _INTEGRALITY_TOL:
            raise RuntimeError("The joint solution violates binary selector bounds.")
        selected = np.rint(selector_values).astype(bool)
        _validate_solution(graph, problem, flux, selected)
        if activity_classifier is not None:
            activity = activity_classifier.read_and_validate(problem, flux)
        elif self.temporal_mode == "flux":
            activity = np.abs(flux) >= (self.eps - _FEASIBILITY_TOL)
        else:
            activity = None
        # Solver integrality tolerances permit a small neighborhood around
        # each discrete value. Commit the nearest exact integer values, then
        # recheck the original rows. This makes subsequent mismatch/floor
        # reporting use the same rounded assignments used for the semantic
        # flux/activity checks above.
        _round_and_validate_discrete_variables(problem)
        union = np.any(selected, axis=1)
        temporal_matrix = selected if activity is None else activity
        switching_counts = {
            pair: int(
                np.count_nonzero(
                    temporal_matrix[:, conditions.index(pair[0])]
                    != temporal_matrix[:, conditions.index(pair[1])]
                )
            )
            for pair in pairs
        }
        total_switches = int(sum(switching_counts.values()))

        state_fits: dict[str, StateFitDiagnostics] = {}
        for index, condition in enumerate(conditions):
            weight = evidence_weights[condition]
            ref_mismatch = independent_mismatch[condition]
            ref_agreement = max(0.0, weight - ref_mismatch) if weight > 0 else 0.0
            joint_mismatch = _expression_value(joint_fit_expr[condition])
            _validate_mismatch(
                joint_mismatch,
                weight,
                stage=f"joint fit for condition {condition!r}",
            )
            joint_agreement = max(0.0, weight - joint_mismatch) if weight > 0 else 0.0
            fit_floor = (1.0 - self.epsilon_fit) * ref_agreement
            margin = joint_agreement - fit_floor if weight > 0 else 0.0
            if weight > 0 and margin < -_FEASIBILITY_TOL * max(1.0, weight):
                raise RuntimeError(
                    f"Condition {condition!r} violates its iMAT fit floor by {-margin:.6g}."
                )
            state_fits[condition] = StateFitDiagnostics(
                evidence_weight=float(weight),
                independent_mismatch=float(ref_mismatch),
                independent_agreement=float(ref_agreement),
                fit_floor=float(fit_floor),
                joint_mismatch=float(joint_mismatch),
                joint_agreement=float(joint_agreement),
                fit_margin=float(margin),
                has_evidence=weight > 0,
                reference_status=reference_status[condition],
            )

        union_size = int(np.count_nonzero(union))
        union_objective = self.lambda_reg * union_size
        if self.temporal_mode == "flux":
            flux_total_variation = float(
                sum(
                    np.abs(
                        flux[:, conditions.index(right)]
                        - flux[:, conditions.index(left)]
                    ).sum()
                    for left, right in pairs
                )
            )
            temporal_objective = self.lambda_temporal * flux_total_variation
        else:
            temporal_objective = self.lambda_temporal * total_switches
        flux_l1 = float(np.abs(flux).sum())
        flux_l2_squared = float(np.square(flux).sum())
        flux_objective = self.lambda_flux * (
            flux_l1 if self.flux_penalty == "l1" else flux_l2_squared
        )
        objective_value = union_objective + temporal_objective + flux_objective
        solved_objective = float(joint_result.value)
        if not np.isfinite(solved_objective) or not np.isclose(
            solved_objective,
            objective_value,
            rtol=1e-6,
            atol=1e-5,
        ):
            raise RuntimeError(
                "Reported union, temporal, and flux objective components do not match "
                f"the solver objective ({objective_value:.8g} vs {solved_objective:.8g})."
            )
        result = TemporalIMATResult(
            problem=problem,
            solver_result=joint_result,
            solver_status=status,
            optimality_proven=status == "optimal",
            reactions=reactions,
            conditions=conditions,
            temporal_pairs=pairs,
            flux=flux,
            selected=selected,
            activity=activity,
            temporal_mode=self.temporal_mode,
            union=union,
            state_fits=state_fits,
            switching_counts=switching_counts,
            total_switches=total_switches,
            union_size=union_size,
            union_objective=float(union_objective),
            temporal_objective=float(temporal_objective),
            lambda_flux=float(self.lambda_flux),
            flux_penalty=self.flux_penalty,
            flux_l1=flux_l1,
            flux_l2_squared=flux_l2_squared,
            flux_objective=float(flux_objective),
            objective_value=float(objective_value),
            warm_start=warm_start_diagnostics,
        )
        self.result = result
        return result


def _finite_nonnegative(name: str, value: float) -> float:
    result = float(value)
    if not np.isfinite(result) or result < 0:
        raise ValueError(f"{name} must be a finite non-negative number")
    return result


def _normalize_inputs(reaction_scores, reaction_bounds):
    if not isinstance(reaction_scores, Mapping) or not reaction_scores:
        raise ValueError("reaction_scores must be a non-empty condition-to-score mapping")
    conditions = tuple(reaction_scores.keys())
    if any(not isinstance(name, str) or not name.strip() for name in conditions):
        raise ValueError("Condition names must be non-empty strings")
    scores = {}
    for condition in conditions:
        state_scores = reaction_scores[condition]
        if not isinstance(state_scores, Mapping):
            raise ValueError(f"reaction_scores[{condition!r}] must be a reaction-to-score mapping")
        scores[condition] = dict(state_scores)

    if reaction_bounds is None:
        bounds = {condition: {} for condition in conditions}
    else:
        if not isinstance(reaction_bounds, Mapping):
            raise ValueError("reaction_bounds must be a condition-to-bounds mapping")
        unknown = set(reaction_bounds) - set(conditions)
        if unknown:
            raise ValueError(f"reaction_bounds contains unknown condition(s): {sorted(unknown)!r}")
        bounds = {}
        for condition in conditions:
            state_bounds = reaction_bounds.get(condition, {})
            if not isinstance(state_bounds, Mapping):
                raise ValueError(f"reaction_bounds[{condition!r}] must be a reaction-to-bounds mapping")
            bounds[condition] = dict(state_bounds)
    return conditions, scores, bounds


def _normalize_pairs(temporal_pairs, conditions):
    if isinstance(temporal_pairs, (str, bytes)):
        raise ValueError("temporal_pairs must be a sequence of condition pairs")
    try:
        raw_pairs = tuple(temporal_pairs)
    except TypeError as exc:
        raise ValueError("temporal_pairs must be a sequence of condition pairs") from exc
    condition_set = set(conditions)
    seen = set()
    pairs = []
    for pair in raw_pairs:
        if isinstance(pair, (str, bytes)):
            raise ValueError("Each temporal pair must contain two condition names")
        try:
            left, right = tuple(pair)
        except (TypeError, ValueError) as exc:
            raise ValueError("Each temporal pair must contain exactly two condition names") from exc
        if not isinstance(left, str) or not isinstance(right, str):
            raise ValueError("Temporal pair condition names must be strings")
        if left not in condition_set or right not in condition_set:
            raise ValueError(f"Unknown condition in temporal pair {(left, right)!r}")
        if left == right:
            raise ValueError(f"Temporal self-pairs are not allowed: {(left, right)!r}")
        key = frozenset((left, right))
        if key in seen:
            raise ValueError(f"Duplicate temporal pair, including reversed duplicates: {(left, right)!r}")
        seen.add(key)
        pairs.append((left, right))
    return tuple(pairs)


def _reaction_labels(graph):
    try:
        raw = tuple(graph.get_attr_from_edges("id"))
    except Exception as exc:
        raise ValueError("The model must expose an 'id' attribute for every reaction edge") from exc
    if not raw or any(not isinstance(reaction, str) or not reaction for reaction in raw):
        raise ValueError("Every reaction edge must have a non-empty string 'id' attribute")
    if len(set(raw)) != len(raw):
        raise ValueError("Reaction edge IDs must be unique")
    return raw


class _TemporalJointIMAT(MultiSampleIMAT):
    """Keep structural selectors available to the joint model at zero weight."""

    def _requires_flux_indicators(self):
        return True


def _imat_method(backend, *, lambda_reg, eps, scale, force_structural_selectors=False):
    method_type = _TemporalJointIMAT if force_structural_selectors else MultiSampleIMAT
    return method_type(
        lambda_reg=lambda_reg,
        beta_reg=0.0,
        eps=eps,
        scale=scale,
        use_bigm_constraints=True,
        backend=backend,
    )


def _fit_objective_name(kind, condition, index):
    return f"imat_fit_{kind}_{str(condition).replace(' ', '_')}_{index}"


def _state_mismatch_expression(problem, condition, index):
    names = {
        _fit_objective_name("pos", condition, index),
        _fit_objective_name("neg", condition, index),
    }
    by_name = {objective.name: objective for objective in problem.objectives}
    expressions = [by_name[name] for name in names if name in by_name]
    if not expressions:
        return problem.backend.Constant(0.0)
    total = expressions[0]
    for expression in expressions[1:]:
        total = total + expression
    return total


def _effective_evidence_weight(scores, *, scale):
    raw_weights = np.asarray([float(value) for value in scores.values()], dtype=float)
    if scale:
        denominator = float(np.abs(raw_weights).sum())
        if denominator > 0:
            raw_weights = raw_weights / denominator * 100.0
    effective = raw_weights[~np.isclose(raw_weights, 0.0)]
    return float(np.abs(effective).sum())


def _scored_reaction_weights(processed_data, condition, *, scale):
    """Return scored reaction IDs and weights in CORNETO's actual feature order."""
    sample = processed_data.samples[condition]
    features = [
        feature
        for feature in sample.features
        if feature.mapping == "edge"
        and (feature.value is not None or "imat_score" in feature.data)
        and (feature.data.get("role") != "objective" or "imat_score" in feature.data)
    ]
    reaction_ids = tuple(feature.id for feature in features)
    weights = np.asarray(
        [float(feature.data.get("imat_score", feature.value)) for feature in features],
        dtype=float,
    )
    if scale:
        denominator = float(np.abs(weights).sum())
        if denominator > 0:
            weights = (weights / denominator) * 100.0
    return tuple(
        (reaction, float(weight))
        for reaction, weight in zip(reaction_ids, weights, strict=True)
        if not np.isclose(weight, 0.0)
    )


def _copy_finite_value(value, *, label):
    if value is None:
        raise RuntimeError(f"Cannot make a warm start: {label} has no numeric solution value.")
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise RuntimeError(f"Cannot make a warm start: {label} contains non-finite values.")
    return np.array(array, dtype=float, copy=True)


def _capture_reference_solution(
    problem, *, processed_data, condition, reactions, scale, eps
):
    """Detach just the per-reaction numeric values used to seed the joint fit."""
    flow = _copy_finite_value(problem.expr.flow.value, label="reference flux")
    flow = flow.reshape(-1)
    if flow.size != len(reactions):
        raise RuntimeError(
            f"Reference flux has {flow.size} entries for {len(reactions)} reactions."
        )

    scored = _scored_reaction_weights(processed_data, condition, scale=scale)
    positive_reactions = tuple(reaction for reaction, weight in scored if weight > 0)
    thresholds: dict[str, tuple[float, float]] = {}
    if positive_reactions:
        flow_name = problem.expr.flow.name
        # MultiSampleIMAT uses these exact semantic role names for the two
        # directional reward binaries created by ThresholdIndicator. The
        # condition index is zero in each one-condition reference problem.
        positive_name = f"{flow_name}_ipos_s0"
        negative_name = f"{flow_name}_ineg_s0"
        if positive_name not in problem.symbols or negative_name not in problem.symbols:
            raise RuntimeError(
                "CORNETO's positive iMAT threshold variables do not match the expected "
                "reference schema (_flow_ipos_s0/_flow_ineg_s0)."
            )
        positive = _copy_finite_value(
            problem.symbols[positive_name].value, label=positive_name
        ).reshape(-1)
        negative = _copy_finite_value(
            problem.symbols[negative_name].value, label=negative_name
        ).reshape(-1)
        if positive.size != len(positive_reactions) or negative.size != len(positive_reactions):
            raise RuntimeError(
                "CORNETO's positive iMAT threshold variable lengths do not match the "
                f"{len(positive_reactions)} positively scored reactions."
            )
        if np.max(np.abs(positive - np.rint(positive)), initial=0.0) > _INTEGRALITY_TOL:
            raise RuntimeError("Reference positive threshold indicators are not integral.")
        if np.max(np.abs(negative - np.rint(negative)), initial=0.0) > _INTEGRALITY_TOL:
            raise RuntimeError("Reference negative threshold indicators are not integral.")
    reaction_index = {reaction: row for row, reaction in enumerate(reactions)}
    thresholds = {}
    mismatch = 0.0
    repaired_thresholds = 0
    thresholds_raw = {}
    if positive_reactions:
        thresholds_raw = {
            reaction: (float(np.rint(pos)), float(np.rint(neg)))
            for reaction, pos, neg in zip(positive_reactions, positive, negative, strict=True)
        }
    for reaction, weight in scored:
        value = float(flow[reaction_index[reaction]])
        if weight > 0:
            old_positive, old_negative = thresholds_raw[reaction]
            positive_value, negative_value = old_positive, old_negative
            threshold_tolerance = min(_FEASIBILITY_TOL, eps / 2.0)
            if value > eps + threshold_tolerance:
                positive_value = 1.0
            elif value < eps - threshold_tolerance:
                positive_value = 0.0
            if value < -eps - threshold_tolerance:
                negative_value = 1.0
            elif value > -eps + threshold_tolerance:
                negative_value = 0.0
            if positive_value + negative_value > 1.0:
                raise RuntimeError(
                    f"Captured high-expression directions overlap for {condition}/{reaction}."
                )
            thresholds[reaction] = (positive_value, negative_value)
            repaired_thresholds += int(
                old_positive > 0.5 and positive_value < 0.5
            ) + int(old_negative > 0.5 and negative_value < 0.5)
            mismatch += weight * (1.0 - positive_value - negative_value)
        elif abs(value) > _FEASIBILITY_TOL:
            mismatch += abs(weight)
    return _ReferenceSolution(
        flux=flow,
        positive_thresholds=thresholds,
        mismatch=float(mismatch),
        repaired_thresholds=repaired_thresholds,
    )


def _repair_reference_solution(
    graph,
    *,
    scores,
    bounds,
    condition,
    original_mismatch,
    evidence_weight,
    epsilon_fit,
    eps,
    scale,
    solver,
    solver_options,
    max_seconds,
    verbosity,
):
    """Resolve a numerically inconsistent reference under its original fit floor."""
    method = _imat_method(
        CvxpyBackend(),
        lambda_reg=0.0,
        eps=eps,
        scale=scale,
    )
    problem = method.build_many(
        graph,
        reaction_scores={condition: scores},
        reaction_bounds={condition: bounds},
    )
    fit_expression = _state_mismatch_expression(problem, condition, 0)
    maximum_mismatch = original_mismatch + epsilon_fit * max(
        0.0, evidence_weight - original_mismatch
    )
    problem.add_constraints(fit_expression <= maximum_mismatch)

    started = time.perf_counter()
    status, result = _solve_checked(
        problem,
        solver,
        max_seconds=max_seconds,
        verbosity=verbosity,
        solver_options=_merge_solver_options(
            {
                "NoRelHeurTime": 0,
                "NoRelHeurWork": 0,
                "IntegralityFocus": 1,
            },
            solver_options,
        ),
        stage=f"threshold feasibility repair for condition {condition!r}",
    )
    elapsed = time.perf_counter() - started
    measured_mismatch = _expression_value(fit_expression)
    _validate_mismatch(
        measured_mismatch,
        evidence_weight,
        stage=f"threshold feasibility repair for condition {condition!r}",
    )
    flow = np.asarray(problem.expr.flow.value, dtype=float).reshape(-1)
    if not np.all(np.isfinite(flow)):
        raise RuntimeError(
            f"The threshold repair for {condition!r} returned non-finite flux values."
        )
    selected = (np.abs(flow) > _FEASIBILITY_TOL)[:, None]
    _validate_solution(graph, problem, flow[:, None], selected)
    _check_problem_constraints(problem)
    solution = _capture_reference_solution(
        problem,
        processed_data=method.processed_data,
        condition=condition,
        reactions=_reaction_labels(graph),
        scale=scale,
        eps=eps,
    )
    mismatch_tolerance = _FEASIBILITY_TOL * max(1.0, evidence_weight)
    if abs(solution.mismatch - measured_mismatch) > mismatch_tolerance:
        raise RuntimeError(
            f"The threshold repair for {condition!r} returned exactized mismatch "
            f"{solution.mismatch:.8g}, but its native iMAT mismatch is "
            f"{measured_mismatch:.8g}; difference exceeds {mismatch_tolerance:.3g}."
        )
    if solution.mismatch > maximum_mismatch + mismatch_tolerance:
        raise RuntimeError(
            f"The threshold repair for {condition!r} failed the unchanged fit floor: "
            f"original mismatch={original_mismatch:.8g}, allowed={maximum_mismatch:.8g}, "
            f"exactized repair mismatch={solution.mismatch:.8g}."
        )
    return solution, status, float(elapsed)


def _regularize_reference_solution(
    graph,
    *,
    reference_problem,
    reference_solution: _ReferenceSolution,
    processed_data,
    condition: str,
    reactions: Sequence[str],
    scale: bool,
    eps: float,
    lambda_flux: float,
    flux_penalty: str,
    solver: str,
    solver_options: Mapping[str, Any],
    max_seconds: int | None,
    verbosity: int,
) -> _ReferenceSolution:
    """Minimize reference flux magnitude while preserving exactized iMAT score states."""
    flow_symbol = reference_problem.expr.flow
    lb = np.asarray(flow_symbol.lb, dtype=float).reshape(-1).copy()
    ub = np.asarray(flow_symbol.ub, dtype=float).reshape(-1).copy()
    ref_flux = np.asarray(reference_solution.flux, dtype=float).reshape(-1)
    reaction_index = {reaction: row for row, reaction in enumerate(reactions)}
    scored = _scored_reaction_weights(processed_data, condition, scale=scale)

    for reaction, weight in scored:
        row = reaction_index[reaction]
        ref_val = float(ref_flux[row])
        if weight > 0:
            pos, neg = reference_solution.positive_thresholds[reaction]
            if pos > 0.5:
                lb[row] = max(lb[row], min(eps, ref_val))
            elif neg > 0.5:
                ub[row] = min(ub[row], max(-eps, ref_val))
            else:
                lb[row] = max(lb[row], min(0.0, ref_val))
                ub[row] = min(ub[row], max(0.0, ref_val))
        elif abs(ref_val) <= _FEASIBILITY_TOL:
            lb[row] = max(lb[row], min(0.0, ref_val))
            ub[row] = min(ub[row], max(0.0, ref_val))

    backend = CvxpyBackend()
    polish = backend.Flow(graph, lb=lb, ub=ub, values=True)
    flow = polish.expr.flow
    if lambda_flux > 0 and flux_penalty == "l2":
        polish.add_objective((flow**2).sum(), name="reference_flux_l2_polish")
    else:
        flux_abs = backend.Variable(
            name="reference_flux_abs",
            shape=flow.shape,
            lb=0.0,
            vartype=VarType.CONTINUOUS,
        )
        polish.add_constraints([flux_abs >= flow, flux_abs >= -flow])
        polish.add_objective(flux_abs.sum(), name="reference_flux_l1_polish")

    _solve_checked(
        polish,
        solver,
        max_seconds=max_seconds,
        verbosity=verbosity,
        solver_options=solver_options,
        stage=f"reference flux regularization for condition {condition!r}",
    )
    polished_flux = _copy_finite_value(
        polish.expr.flow.value, label="regularized reference flux"
    ).reshape(-1)
    mismatch = 0.0
    for reaction, weight in scored:
        value = float(polished_flux[reaction_index[reaction]])
        if weight > 0:
            pos, neg = reference_solution.positive_thresholds[reaction]
            mismatch += weight * (1.0 - pos - neg)
        elif abs(value) > _FEASIBILITY_TOL:
            mismatch += abs(weight)
    return _ReferenceSolution(
        flux=polished_flux,
        positive_thresholds=reference_solution.positive_thresholds,
        mismatch=float(mismatch),
        repaired_thresholds=reference_solution.repaired_thresholds,
    )


def _assign_symbol_start(problem, name, values):
    symbol = problem.symbols.get(name)
    if symbol is None or not getattr(symbol, "is_variable", False):
        raise RuntimeError(f"Cannot assemble a warm start: expected decision variable {name!r}.")
    if symbol._vartype == VarType.BINARY:
        numeric = np.asarray(values, dtype=float)
        if not np.all(np.isfinite(numeric)) or np.max(
            np.abs(numeric - np.rint(numeric)), initial=0.0
        ) > _INTEGRALITY_TOL:
            raise RuntimeError(f"Warm-start values for binary variable {name!r} are not integral.")
        array = np.asarray(np.rint(numeric), dtype=bool)
    elif symbol._vartype == VarType.INTEGER:
        numeric = np.asarray(values, dtype=float)
        if not np.all(np.isfinite(numeric)) or np.max(
            np.abs(numeric - np.rint(numeric)), initial=0.0
        ) > _INTEGRALITY_TOL:
            raise RuntimeError(f"Warm-start values for integer variable {name!r} are not integral.")
        array = np.asarray(np.rint(numeric), dtype=int)
    else:
        array = np.asarray(values, dtype=float)
    if array.size != int(np.prod(symbol.shape, dtype=int)):
        raise RuntimeError(
            f"Warm-start values for {name!r} have {array.size} entries, "
            f"expected shape {symbol.shape}."
        )
    symbol.value = np.array(array.reshape(symbol.shape), copy=True)


def _assemble_and_validate_warm_start(
    problem,
    *,
    graph,
    solver,
    references,
    conditions,
    reaction_scores,
    processed_data,
    independent_mismatch,
    evidence_weights,
    reactions,
    temporal_pairs,
    scale,
    lambda_reg,
    lambda_temporal,
    lambda_flux,
    flux_penalty,
    epsilon_fit,
    reference_mismatch_before_repair,
    reference_repair_statuses,
    reference_repair_seconds,
    repaired_reference_thresholds,
    temporal_mode,
    activity_classifier,
):
    """Compose independent numeric solutions into a fully checked joint start."""
    if tuple(references) != tuple(conditions):
        raise RuntimeError("Warm-start references do not match the joint condition order.")
    n_reactions = len(reactions)
    n_conditions = len(conditions)
    reaction_index = {reaction: index for index, reaction in enumerate(reactions)}
    flux = np.column_stack([references[condition].flux for condition in conditions])
    if flux.shape != (n_reactions, n_conditions):
        raise RuntimeError(
            f"Composed reference flux has shape {flux.shape}; expected "
            f"{(n_reactions, n_conditions)}."
        )

    selected = np.abs(flux) > _FEASIBILITY_TOL
    for condition_index, condition in enumerate(conditions):
        reference = references[condition]
        expected_positive = tuple(
            reaction
            for reaction, weight in _scored_reaction_weights(
                processed_data, condition, scale=scale
            )
            if weight > 0
        )
        if set(reference.positive_thresholds) != set(expected_positive):
            raise RuntimeError(
                f"Warm-start threshold mapping for condition {condition!r} does not match "
                "its positively scored reaction IDs."
            )
        for reaction, (positive, negative) in reference.positive_thresholds.items():
            row = reaction_index[reaction]
            selected[row, condition_index] |= (
                positive > 0.5 or negative > 0.5
            )

    flow_symbol = problem.expr.flow
    _assign_symbol_start(problem, flow_symbol.name, flux)
    if lambda_flux > 0 and flux_penalty == "l1":
        flux_abs_symbol = problem.symbols.get("flux_abs")
        if flux_abs_symbol is None:
            raise RuntimeError("Positive lambda_flux has no flux_abs epigraph variable.")
        _assign_symbol_start(problem, flux_abs_symbol.name, np.abs(flux))
    elif "flux_abs" in problem.symbols:
        raise RuntimeError("A non-L1 flux penalty unexpectedly added a flux_abs variable.")

    selector_symbol = problem.symbols.get("edge_has_flux")
    if selector_symbol is None:
        raise RuntimeError(
            "The joint CORNETO problem has no structural edge_has_flux selectors to seed."
        )
    for condition_index, condition in enumerate(conditions):
        positive_reactions = tuple(
            reaction
            for reaction, weight in _scored_reaction_weights(
                processed_data, condition, scale=scale
            )
            if weight > 0
        )
        if not positive_reactions:
            continue
        flow_name = flow_symbol.name
        positive_name = f"{flow_name}_ipos_s{condition_index}"
        negative_name = f"{flow_name}_ineg_s{condition_index}"
        if positive_name not in problem.symbols or negative_name not in problem.symbols:
            raise RuntimeError(
                f"Joint CORNETO threshold variable schema is missing the directional "
                f"reward variables for condition index {condition_index}."
            )
        positive_values = []
        negative_values = []
        for reaction in positive_reactions:
            try:
                positive, negative = references[condition].positive_thresholds[reaction]
            except KeyError as exc:
                raise RuntimeError(
                    f"Missing detached threshold value for {condition}/{reaction}."
                ) from exc
            positive_values.append(positive)
            negative_values.append(negative)
        _assign_symbol_start(problem, positive_name, positive_values)
        _assign_symbol_start(problem, negative_name, negative_values)

    activity_start = None
    if activity_classifier is not None:
        activity_start = activity_classifier.seed_additional_variables(
            problem, flux, references
        )
        selected |= activity_start
    _assign_symbol_start(problem, selector_symbol.name, selected)

    # CORNETO's union indicator is registered as both a solver symbol and an
    # expression. Assign the symbol handle so CVXPY receives the value.
    union_symbol = problem.symbols.get("_edge_has_flux_OR")
    if union_symbol is not None:
        _assign_symbol_start(problem, union_symbol.name, np.any(selected, axis=1))

    # TotalVariation exposes a writable absolute-value auxiliary only when its
    # penalty is enabled. Its signed differences are derived expressions.
    if temporal_mode == "selection":
        temporal_start_matrix = selected
        tv_prefix = "selector"
    elif temporal_mode == "flux":
        temporal_start_matrix = flux
        tv_prefix = "flux"
    else:
        temporal_start_matrix = activity_start
        tv_prefix = "activity"
    abs_bound_name = f"temporal_{tv_prefix}_abs_bound"
    if abs_bound_name in problem.symbols:
        if temporal_start_matrix is None:
            raise RuntimeError("Activity-mode warm start has no activity classifier values.")
        differences = np.column_stack(
            [
                np.abs(
                    np.asarray(temporal_start_matrix[:, conditions.index(right)], dtype=float)
                    - np.asarray(temporal_start_matrix[:, conditions.index(left)], dtype=float)
                )
                for left, right in temporal_pairs
            ]
        )
        _assign_symbol_start(problem, abs_bound_name, differences)
    elif temporal_pairs and lambda_temporal > 0:
        raise RuntimeError("Temporal variation has no writable absolute-bound variable.")

    decision_symbols = [
        symbol
        for symbol in problem.symbols.values()
        if getattr(symbol, "is_variable", False)
    ]
    missing = [symbol.name for symbol in decision_symbols if symbol.value is None]
    if missing:
        raise RuntimeError(
            "Warm-start assembly left decision variables uninitialized: "
            + ", ".join(missing)
        )

    maximum_bound_violation = 0.0
    for symbol in decision_symbols:
        values = _copy_finite_value(symbol.value, label=f"warm-start variable {symbol.name}")
        lower = np.asarray(symbol.lb, dtype=float)
        upper = np.asarray(symbol.ub, dtype=float)
        lower_violation = np.maximum(lower - values, 0.0)
        upper_violation = np.maximum(values - upper, 0.0)
        maximum_bound_violation = max(
            maximum_bound_violation,
            float(np.max(lower_violation, initial=0.0)),
            float(np.max(upper_violation, initial=0.0)),
        )
        if symbol._vartype == VarType.BINARY and np.max(
            np.abs(values - np.rint(values)), initial=0.0
        ) > _INTEGRALITY_TOL:
            raise RuntimeError(f"Warm-start variable {symbol.name!r} is not integral.")
    if maximum_bound_violation > _FEASIBILITY_TOL:
        raise RuntimeError(
            f"Warm-start variable bounds are violated by {maximum_bound_violation:.6g}."
        )

    selector_values = np.asarray(selector_symbol.value, dtype=float)
    selector_matrix = _as_reaction_condition_matrix(
        selector_values, n_reactions, n_conditions
    )
    flux_matrix = _as_reaction_condition_matrix(
        np.asarray(flow_symbol.value, dtype=float), n_reactions, n_conditions
    )
    _validate_solution(
        graph,
        problem,
        flux_matrix,
        np.rint(selector_matrix).astype(bool),
    )
    if activity_classifier is not None:
        validated_activity_start = activity_classifier.read_and_validate(
            problem, flux_matrix
        )
        if not np.array_equal(validated_activity_start, activity_start):
            raise RuntimeError("Activity classifier warm-start directions changed during validation.")

    fit_floor_margins = {}
    repaired_reference_mismatch = {}
    candidate_mismatch_by_condition = {}
    allowed_mismatch_by_condition = {}
    floor_violations = []
    for index, condition in enumerate(conditions):
        mismatch = _expression_value(_state_mismatch_expression(problem, condition, index))
        weight = evidence_weights[condition]
        candidate_mismatch_by_condition[condition] = float(mismatch)
        if weight <= 0:
            fit_floor_margins[condition] = float("nan")
            repaired_reference_mismatch[condition] = references[condition].mismatch
            allowed_mismatch_by_condition[condition] = float("nan")
            continue
        _validate_mismatch(mismatch, weight, stage=f"warm-start condition {condition!r}")
        reference_mismatch = independent_mismatch[condition]
        maximum_mismatch = reference_mismatch + epsilon_fit * max(
            0.0, weight - reference_mismatch
        )
        allowed_mismatch_by_condition[condition] = float(maximum_mismatch)
        margin = maximum_mismatch - mismatch
        fit_floor_margins[condition] = float(margin)
        repaired_reference_mismatch[condition] = references[condition].mismatch
        if margin < -_FEASIBILITY_TOL * max(1.0, weight):
            floor_violations.append(
                f"{condition!r}: original reference mismatch={reference_mismatch:.8g}, "
                f"exactized reference mismatch={references[condition].mismatch:.8g}, "
                f"candidate mismatch={mismatch:.8g}, allowed mismatch={maximum_mismatch:.8g}, "
                f"floor margin={margin:.8g}"
            )

    if floor_violations:
        details = "; ".join(
            f"{condition!r}: original reference mismatch="
            f"{independent_mismatch[condition]:.8g}, exactized reference mismatch="
            f"{repaired_reference_mismatch[condition]:.8g}, exactization increase="
            f"{repaired_reference_mismatch[condition] - independent_mismatch[condition]:.8g}, "
            f"candidate mismatch={candidate_mismatch_by_condition[condition]:.8g}, "
            f"allowed mismatch={allowed_mismatch_by_condition[condition]:.8g}, "
            f"floor margin={fit_floor_margins[condition]:.8g}"
            for condition in conditions
        )
        raise _WarmStartFitFloorError(
            "Repaired warm start violates unchanged iMAT fit floor(s); no start was "
            "submitted and the original independent objectives/constraints were retained. "
            "Per-condition original reference, exactized reference, candidate mismatch and "
            "floor-margin details: "
            + details
            + "."
        )

    maximum_constraint_violation = _check_problem_constraints(problem)

    union = np.any(selected, axis=1)
    if temporal_mode == "flux":
        temporal_variation = float(
            sum(
                np.abs(
                    flux[:, conditions.index(right)]
                    - flux[:, conditions.index(left)]
                ).sum()
                for left, right in temporal_pairs
            )
        )
    else:
        if temporal_start_matrix is None:
            temporal_start_matrix = selected
        temporal_variation = float(
            sum(
                int(
                    np.count_nonzero(
                        temporal_start_matrix[:, conditions.index(left)]
                        != temporal_start_matrix[:, conditions.index(right)]
                    )
                )
                for left, right in temporal_pairs
            )
        )
    flux_l1 = float(np.abs(flux).sum())
    flux_l2_squared = float(np.square(flux).sum())
    flux_objective = float(
        lambda_flux * (flux_l1 if flux_penalty == "l1" else flux_l2_squared)
    )
    objective = float(
        lambda_reg * np.count_nonzero(union)
        + lambda_temporal * temporal_variation
        + flux_objective
    )
    if not np.isfinite(objective):
        raise RuntimeError("The assembled warm-start objective is not finite.")
    cvxpy_coverage = _inspect_cvxpy_warm_start(
        problem,
        solver=solver,
        inspect_weighted_objective=(lambda_flux > 0 and flux_penalty == "l2"),
    )
    return WarmStartDiagnostics(
        decision_variables=sum(
            int(np.prod(symbol.shape, dtype=int)) for symbol in decision_symbols
        ),
        initialized_variables=sum(
            int(np.prod(symbol.shape, dtype=int)) for symbol in decision_symbols
        ),
        decision_symbols=len(decision_symbols),
        objective_value=objective,
        flux_penalty=flux_penalty,
        flux_l1=flux_l1,
        flux_l2_squared=flux_l2_squared,
        flux_objective=flux_objective,
        maximum_constraint_violation=maximum_constraint_violation,
        maximum_bound_violation=maximum_bound_violation,
        **cvxpy_coverage,
        repaired_reference_thresholds=repaired_reference_thresholds,
        independent_reference_mismatch=dict(independent_mismatch),
        reference_mismatch_before_repair=dict(reference_mismatch_before_repair),
        reference_mismatch_after_repair=repaired_reference_mismatch,
        fit_floor_margins=fit_floor_margins,
        reference_repair_statuses=dict(reference_repair_statuses),
        reference_repair_seconds=dict(reference_repair_seconds),
    )


def _inspect_cvxpy_warm_start(problem, *, solver, inspect_weighted_objective=False):
    """Check canonical warm-start values, optionally through the actual QP route."""
    import cvxpy as cp

    extras = []
    for symbol in problem.symbols.values():
        if symbol._provided_lb is not None:
            extras.append(symbol >= symbol.lb)
        if symbol._provided_ub is not None:
            extras.append(symbol <= symbol.ub)
    constraints = [constraint.e for constraint in extras + problem.constraints]
    if inspect_weighted_objective:
        # Mirror CvxpyBackend.solve's weighted objective sum. This makes the
        # start check exercise GUROBI's quadratic-objective canonicalization,
        # rather than the LP-only feasibility route used for linear starts.
        if problem.objectives is not None and len(problem.objectives) > 1:
            weighted_objective = sum(
                problem.weights[index] * problem.objectives[index]
                for index in range(len(problem.objectives))
            )
        elif problem.objectives and problem.weights[0] != 0:
            weighted_objective = problem.weights[0] * problem.objectives[0]
        else:
            weighted_objective = None
        objective = cp.Constant(0) if weighted_objective is None else weighted_objective.e
    else:
        objective = cp.Constant(0)
    cvxpy_problem = cp.Problem(cp.Minimize(objective), constraints)
    if inspect_weighted_objective and (
        not cvxpy_problem.is_qp()
        or not all(constraint.expr.is_affine() for constraint in cvxpy_problem.constraints)
    ):
        raise RuntimeError("L2 warm-start inspection did not produce a QP with linear constraints.")
    cvxpy_variables = cvxpy_problem.variables()
    missing = [variable.name() for variable in cvxpy_variables if variable.value is None]
    if missing:
        raise RuntimeError(
            "CVXPY would receive uninitialized warm-start variables: " + ", ".join(missing)
        )
    try:
        data, _, _ = cvxpy_problem.get_problem_data(solver=solver.upper())
    except Exception as exc:
        raise RuntimeError(
            f"Could not inspect CVXPY's {solver.upper()} warm-start vector: {exc}"
        ) from exc
    initial_values = data.get("init_value")
    if initial_values is None:
        raise RuntimeError(
            f"CVXPY did not produce an initial-value vector for {solver.upper()}."
        )
    initial_values = np.asarray(initial_values, dtype=float).reshape(-1)
    if not np.all(np.isfinite(initial_values)):
        raise RuntimeError("CVXPY's canonical warm-start vector contains non-finite values.")
    cvxpy_variable_count = int(sum(variable.size for variable in cvxpy_variables))
    canonical_variable_count = int(initial_values.size)
    auxiliary_count = canonical_variable_count - cvxpy_variable_count
    if auxiliary_count < 0:
        raise RuntimeError(
            "CVXPY's canonical warm-start vector is shorter than its declared variables."
        )
    if auxiliary_count:
        raise RuntimeError(
            f"CVXPY introduced {auxiliary_count} canonical auxiliary variable(s) without "
            "reaction/condition mappings. Refusing to guess their warm-start values."
        )
    quadratic_objective = data.get("P")
    quadratic_nonzeros = int(getattr(quadratic_objective, "nnz", 0))
    if inspect_weighted_objective and quadratic_nonzeros <= 0:
        raise RuntimeError(
            "L2 warm-start inspection did not produce a native quadratic objective matrix."
        )
    # Linear CORNETO iMAT constraints should not add canonical decision
    # variables. An unexpected auxiliary is measured above and fails closed.
    return {
        "cvxpy_variables": cvxpy_variable_count,
        "canonical_variables": canonical_variable_count,
        "canonical_auxiliary_variables": auxiliary_count,
        "uninitialized_canonical_variables": 0,
        "canonical_objective_kind": "quadratic" if quadratic_nonzeros else "linear",
        "canonical_quadratic_objective_nonzeros": quadratic_nonzeros,
    }


def _solve_checked(
    problem,
    solver,
    *,
    max_seconds,
    verbosity,
    solver_options,
    stage,
    allow_feasible_incumbent=False,
    warm_start=False,
):
    options = {
        key: dict(value) if isinstance(value, Mapping) else value
        for key, value in solver_options.items()
    }
    if isinstance(solver, str) and solver.upper() == "GUROBI":
        if not any(str(key).casefold() == "intfeastol".casefold() for key in options):
            options["IntFeasTol"] = _DEFAULT_GUROBI_INT_FEAS_TOL
    try:
        result = problem.solve(
            solver=solver,
            max_seconds=max_seconds,
            verbosity=verbosity,
            warm_start=warm_start,
            **options,
        )
    except Exception as exc:
        raise RuntimeError(f"Solver failed during {stage}: {exc}") from exc
    status = str(getattr(result, "status", "unknown")).lower()
    if status != "optimal":
        if (
            allow_feasible_incumbent
            and status in {"user_limit", "optimal_inaccurate", "infeasible_inaccurate"}
            and _has_gurobi_time_limit_incumbent(result, solver)
        ):
            return status, result
        details = _format_solver_diagnostics(result, solver)
        raise RuntimeError(
            f"Solver returned status {status!r} during {stage}; an optimal solution is required."
            f"{details}"
        )
    return status, result


def _merge_solver_options(*option_mappings):
    """Merge backend options without case-insensitive duplicate Gurobi keys."""
    merged = {}
    keys_by_casefold = {}
    for mapping in option_mappings:
        for key, value in mapping.items():
            normalized = str(key).casefold()
            previous = keys_by_casefold.get(normalized)
            if previous is not None:
                del merged[previous]
            merged[key] = value
            keys_by_casefold[normalized] = key
    return merged


def _gurobi_extra_stats(solver_result, solver):
    if not isinstance(solver, str) or solver.upper() != "GUROBI":
        return None
    solver_stats = getattr(solver_result, "solver_stats", None)
    return getattr(solver_stats, "extra_stats", None)


def _native_stat(extra_stats, name):
    try:
        value = getattr(extra_stats, name)
    except Exception:
        return None
    if value is None:
        return None
    try:
        if isinstance(value, (float, np.floating)) and not np.isfinite(value):
            return None
    except TypeError:
        return None
    return value


def _has_gurobi_time_limit_incumbent(solver_result, solver):
    extra_stats = _gurobi_extra_stats(solver_result, solver)
    if extra_stats is None:
        return False
    status_code = _native_stat(extra_stats, "Status")
    try:
        from gurobipy import GRB

        time_limit_code = GRB.TIME_LIMIT
    except ImportError:
        time_limit_code = 9  # Gurobi's TIME_LIMIT status code.
    if status_code != time_limit_code:
        return False
    try:
        solution_count = int(_native_stat(extra_stats, "SolCount"))
        incumbent_objective = float(_native_stat(extra_stats, "ObjVal"))
    except (TypeError, ValueError, OverflowError):
        return False
    return solution_count > 0 and np.isfinite(incumbent_objective)


def _format_solver_diagnostics(solver_result, solver):
    extra_stats = _gurobi_extra_stats(solver_result, solver)
    if extra_stats is None:
        return ""
    stats = {
        name: _native_stat(extra_stats, name)
        for name in ("Status", "SolCount", "ObjVal", "ObjBound", "MIPGap", "Runtime")
    }
    if not any(value is not None for value in stats.values()):
        return " Native Gurobi statistics were unavailable."
    formatted = ", ".join(
        f"{name}={value if value is not None else 'unavailable'}"
        for name, value in stats.items()
    )
    return f" Native Gurobi statistics: {formatted}."


def _expression_value(expression):
    value = getattr(expression, "value", None)
    if value is None:
        raise RuntimeError("A solved iMAT mismatch expression has no numerical value")
    values = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(values)):
        raise RuntimeError("A solved iMAT mismatch expression contains non-finite values")
    return float(values.sum())


def _validate_mismatch(mismatch, evidence_weight, *, stage):
    if not np.isfinite(mismatch):
        raise RuntimeError(f"The mismatch expression for {stage} is not finite: {mismatch!r}.")
    tolerance = _FEASIBILITY_TOL * max(1.0, evidence_weight)
    if mismatch < -tolerance or mismatch > evidence_weight + tolerance:
        raise RuntimeError(
            f"The mismatch expression for {stage} is outside its expected range "
            f"[0, {evidence_weight:.8g}]: {mismatch:.8g}."
        )


def _as_reaction_condition_matrix(values, n_reactions, n_conditions):
    result = np.asarray(values, dtype=float)
    if n_conditions == 1 and result.shape == (n_reactions,):
        result = result[:, None]
    if result.shape != (n_reactions, n_conditions):
        raise RuntimeError(
            f"CORNETO returned array shape {result.shape}, expected {(n_reactions, n_conditions)}"
        )
    return result


def _validate_solution(graph, problem, flux, selected):
    flow_symbol = problem.expr.flow
    lower = np.asarray(flow_symbol.lb, dtype=float)
    upper = np.asarray(flow_symbol.ub, dtype=float)
    n_reactions, n_conditions = flux.shape
    if n_conditions == 1 and lower.shape == (n_reactions,):
        lower = lower[:, None]
        upper = upper[:, None]
    if lower.shape != flux.shape or upper.shape != flux.shape:
        raise RuntimeError("Could not verify reaction bounds: CORNETO flow-bound shape mismatch")
    if np.any(flux < lower - _FEASIBILITY_TOL) or np.any(flux > upper + _FEASIBILITY_TOL):
        raise RuntimeError("The solved fluxes violate reaction bounds beyond numerical tolerance")

    data, (rows, columns) = graph.get_vertex_incidence_matrix_as_lists(values=True)
    incidence = coo_matrix(
        (
            np.asarray(data, dtype=float),
            (np.asarray(rows, dtype=int), np.asarray(columns, dtype=int)),
        ),
        shape=(len(graph.V), graph.num_edges),
    ).tocsr()
    residual = incidence @ flux
    if np.max(np.abs(residual), initial=0.0) > _FEASIBILITY_TOL:
        raise RuntimeError("The solved fluxes violate mass balance beyond numerical tolerance")
    if selected.shape != flux.shape:
        raise RuntimeError("The solved selection matrix does not cover all reactions and conditions")
    if np.any((~selected) & (np.abs(flux) > _FEASIBILITY_TOL)):
        raise RuntimeError("An unselected reaction carries flux beyond numerical tolerance")

    union_expression = getattr(problem.expr, "edge_has_flux_OR", None)
    union_values = getattr(union_expression, "value", None)
    if union_values is not None:
        union_values = np.asarray(union_values, dtype=float).reshape(-1)
        expected_union = np.any(selected, axis=1)
        if not np.all(np.isfinite(union_values)):
            raise RuntimeError("The joint reaction-union selectors contain non-finite values")
        if np.max(np.abs(union_values - np.rint(union_values)), initial=0.0) > _INTEGRALITY_TOL:
            raise RuntimeError("The joint reaction-union selectors violate integrality tolerance")
        if np.any(union_values < -_INTEGRALITY_TOL) or np.any(
            union_values > 1.0 + _INTEGRALITY_TOL
        ):
            raise RuntimeError("The joint reaction-union selectors violate binary bounds")
        if not np.array_equal(np.rint(union_values).astype(bool), expected_union):
            raise RuntimeError("The joint reaction-union selectors do not match condition selections")


def _validate_problem_constraints(problem):
    _check_problem_constraints(problem)


def _round_and_validate_discrete_variables(problem):
    """Round every solved integer symbol and validate the exact assignment."""
    for symbol in problem.symbols.values():
        if not getattr(symbol, "is_variable", False) or symbol._vartype not in (
            VarType.BINARY,
            VarType.INTEGER,
        ):
            continue
        values = _copy_finite_value(symbol.value, label=f"solved variable {symbol.name}")
        if np.max(np.abs(values - np.rint(values)), initial=0.0) > _INTEGRALITY_TOL:
            raise RuntimeError(f"Solved discrete variable {symbol.name!r} is not integral.")
        rounded = np.rint(values)
        lower = np.asarray(symbol.lb, dtype=float)
        upper = np.asarray(symbol.ub, dtype=float)
        if np.any(rounded < lower - _FEASIBILITY_TOL) or np.any(
            rounded > upper + _FEASIBILITY_TOL
        ):
            raise RuntimeError(
                f"Rounded discrete variable {symbol.name!r} violates its bounds."
            )
        symbol.value = np.array(rounded.reshape(symbol.shape), copy=True)
    return _validate_problem_constraints(problem)


def _check_problem_constraints(problem):
    maximum_violation = 0.0
    for index, constraint in enumerate(problem.constraints):
        violation_method = getattr(constraint.e, "violation", None)
        if violation_method is None:
            raise RuntimeError(f"Could not verify solved constraint {index}")
        violation = violation_method()
        if violation is None:
            raise RuntimeError(f"Solved constraint {index} has no numerical residual")
        violation = np.asarray(violation, dtype=float)
        maximum = float(np.max(violation, initial=0.0))
        maximum_violation = max(maximum_violation, maximum)
        if not np.all(np.isfinite(violation)) or maximum > _FEASIBILITY_TOL:
            raise RuntimeError(
                f"Solved constraint {index} is violated by {maximum:.6g}, above "
                f"the {_FEASIBILITY_TOL:g} tolerance. Constraint: {constraint.e}"
            )
    return maximum_violation


__all__ = ["StateFitDiagnostics", "TemporalIMATResult", "TemporalMultiSampleIMAT"]
