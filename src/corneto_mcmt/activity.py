"""Threshold activity classification for temporal iMAT models.

The CORNETO structural selector answers whether a reaction may carry flux. This
module builds a separate expression from the native positive-score iMAT reward
binaries and packed binaries for the remaining reaction/state/sign entries.
The packed representation keeps the additional variable count proportional to
feasible directions rather than allocating a dense binary matrix.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
from corneto.backend import VarType
from scipy.sparse import coo_matrix


@dataclass(frozen=True)
class _BinaryAddress:
    """Location of one direction indicator in a backend binary vector."""

    name: str
    offset: int
    native_reward: bool


@dataclass(frozen=True)
class ActivityClassifier:
    """Reusable activity expression, classifier mappings, and warm-start hooks."""

    reactions: tuple[str, ...]
    conditions: tuple[str, ...]
    eps: float
    expression: object
    positive_addresses: Mapping[tuple[int, int], _BinaryAddress]
    negative_addresses: Mapping[tuple[int, int], _BinaryAddress]
    additional_positive_cells: tuple[tuple[int, int], ...]
    additional_negative_cells: tuple[tuple[int, int], ...]
    positive_name: str | None
    negative_name: str | None

    def seed_additional_variables(self, problem, flux: np.ndarray, references: Mapping) -> np.ndarray:
        """Set new direction binaries and return the rounded activity start."""
        positive = np.zeros((len(self.reactions), len(self.conditions)), dtype=bool)
        negative = np.zeros_like(positive)

        for cells, target, direction in (
            (self.additional_positive_cells, positive, 1),
            (self.additional_negative_cells, negative, -1),
        ):
            values = []
            for reaction_index, condition_index in cells:
                value = flux[reaction_index, condition_index]
                values.append(value >= self.eps if direction > 0 else value <= -self.eps)
                target[reaction_index, condition_index] = values[-1]
            variable_name = self.positive_name if direction > 0 else self.negative_name
            if cells:
                _assign_binary_start(problem, variable_name, values)

        for direction, addresses, target in (
            (1, self.positive_addresses, positive),
            (-1, self.negative_addresses, negative),
        ):
            for (row, column), address in addresses.items():
                if not address.native_reward:
                    continue
                reaction = self.reactions[row]
                condition = self.conditions[column]
                try:
                    high_positive, high_negative = references[condition].positive_thresholds[reaction]
                except KeyError as exc:
                    raise RuntimeError(
                        f"Missing captured high-expression classifier value for {condition}/{reaction}."
                    ) from exc
                target[row, column] = high_positive if direction > 0 else high_negative

        return positive | negative

    def read_and_validate(self, problem, flux: np.ndarray) -> np.ndarray:
        """Round and validate both direction vectors against solved fluxes."""
        direction_values: dict[str, np.ndarray] = {}
        all_names = {
            address.name
            for mapping in (self.positive_addresses, self.negative_addresses)
            for address in mapping.values()
        }
        for name in all_names:
            symbol = problem.symbols.get(name)
            if symbol is None or getattr(symbol, "_vartype", None) != VarType.BINARY:
                raise RuntimeError(f"Activity classifier variable {name!r} is missing or is not binary.")
            values = _finite_array(symbol.value, f"activity classifier variable {name}").reshape(-1)
            if np.max(np.abs(values - np.rint(values)), initial=0.0) > _integrality_tolerance():
                raise RuntimeError(f"Activity classifier variable {name!r} is not integral.")
            if np.min(values, initial=0.0) < -_integrality_tolerance() or np.max(
                values, initial=0.0
            ) > 1.0 + _integrality_tolerance():
                raise RuntimeError(f"Activity classifier variable {name!r} is outside binary bounds.")
            direction_values[name] = np.rint(values).astype(int)

        active = np.zeros((len(self.reactions), len(self.conditions)), dtype=bool)
        flux_tolerance = min(1e-6, self.eps / 2.0)
        for row in range(len(self.reactions)):
            for column in range(len(self.conditions)):
                positive_address = self.positive_addresses.get((row, column))
                negative_address = self.negative_addresses.get((row, column))
                p = (
                    int(direction_values[positive_address.name][positive_address.offset])
                    if positive_address is not None
                    else 0
                )
                n = (
                    int(direction_values[negative_address.name][negative_address.offset])
                    if negative_address is not None
                    else 0
                )
                value = float(flux[row, column])
                if p + n > 1:
                    raise RuntimeError(
                        f"Activity directions overlap for {self.conditions[column]}/"
                        f"{self.reactions[row]}."
                    )
                if p and value < self.eps - flux_tolerance:
                    raise RuntimeError(
                        f"Active positive classifier for {self.conditions[column]}/"
                        f"{self.reactions[row]} has flux {value:.8g}, below eps={self.eps:.8g}."
                    )
                if n and value > -self.eps + flux_tolerance:
                    raise RuntimeError(
                        f"Active negative classifier for {self.conditions[column]}/"
                        f"{self.reactions[row]} has flux {value:.8g}, above -eps={-self.eps:.8g}."
                    )
                if not p and not n and abs(value) > self.eps + flux_tolerance:
                    raise RuntimeError(
                        f"Inactive classifier for {self.conditions[column]}/"
                        f"{self.reactions[row]} has clearly above-threshold flux {value:.8g}."
                    )
                active[row, column] = bool(p or n)
        expression_values = _finite_array(self.expression.value, "registered activity expression")
        if expression_values.shape != active.shape or not np.array_equal(
            np.rint(expression_values).astype(bool), active
        ):
            raise RuntimeError("Registered activity expression does not match its direction binaries.")
        return active


def build_activity_classifier(
    problem,
    *,
    reactions: Sequence[str],
    conditions: Sequence[str],
    positive_reactions_by_condition: Mapping[str, Sequence[str]],
    eps: float,
) -> ActivityClassifier:
    """Add finite-bound activity indicators and register their full matrix.

    Positive-score reactions reuse CORNETO's native threshold binaries. Other
    entries receive packed binaries only for directions admitted by their
    original flux bounds.
    """
    reactions = tuple(reactions)
    conditions = tuple(conditions)
    n_reactions, n_conditions = len(reactions), len(conditions)
    reaction_index = {reaction: index for index, reaction in enumerate(reactions)}
    condition_index = {condition: index for index, condition in enumerate(conditions)}
    flow = problem.expr.flow
    lower, upper = _bound_matrices(flow.lb, flow.ub, n_reactions, n_conditions)
    if not np.all(np.isfinite(lower)) or not np.all(np.isfinite(upper)):
        raise ValueError("Activity classification requires finite reaction bounds.")

    positive_addresses: dict[tuple[int, int], _BinaryAddress] = {}
    negative_addresses: dict[tuple[int, int], _BinaryAddress] = {}
    native_positive_by_condition: dict[str, tuple[str, str, tuple[int, ...]]] = {}
    native_classifier_constraints = []

    for condition in conditions:
        column = condition_index[condition]
        high_reactions = tuple(positive_reactions_by_condition.get(condition, ()))
        if len(set(high_reactions)) != len(high_reactions):
            raise ValueError(f"Positive-score reaction IDs are repeated for condition {condition!r}.")
        unknown = set(high_reactions) - set(reaction_index)
        if unknown:
            raise ValueError(
                f"Positive-score reactions for {condition!r} are absent from the graph: {sorted(unknown)!r}."
            )
        if not high_reactions:
            native_positive_by_condition[condition] = ("", "", ())
            continue
        flow_name = flow.name
        positive_name = f"{flow_name}_ipos_s{column}"
        negative_name = f"{flow_name}_ineg_s{column}"
        positive_symbol = problem.symbols.get(positive_name)
        negative_symbol = problem.symbols.get(negative_name)
        if positive_symbol is None or negative_symbol is None:
            raise RuntimeError(
                "CORNETO's positive iMAT ThresholdIndicator variables are missing; "
                f"expected {positive_name!r} and {negative_name!r}."
            )
        if positive_symbol.shape != (len(high_reactions),) or negative_symbol.shape != (
            len(high_reactions),
        ):
            raise RuntimeError(
                f"Native high-expression classifier shape for {condition!r} does not match "
                f"its {len(high_reactions)} positive-score reactions."
            )
        native_positive_by_condition[condition] = (
            positive_name,
            negative_name,
            tuple(reaction_index[reaction] for reaction in high_reactions),
        )
        native_rows = np.asarray(
            [reaction_index[reaction] for reaction in high_reactions], dtype=int
        )
        native_columns = np.full(len(high_reactions), column, dtype=int)
        native_flow = _select_cells(
            flow,
            tuple((int(row), int(col)) for row, col in zip(native_rows, native_columns, strict=True)),
            n_conditions,
        )
        native_lower = lower[native_rows, native_columns]
        native_upper = upper[native_rows, native_columns]
        native_positive = positive_symbol
        native_negative = negative_symbol
        native_classifier_constraints.extend(
            [
                native_positive + native_negative <= 1,
                native_flow
                <= np.minimum(native_upper, eps)
                + native_positive.multiply(np.maximum(native_upper - eps, 0.0)),
                native_flow
                >= np.maximum(native_lower, -eps)
                + native_negative.multiply(np.minimum(native_lower + eps, 0.0)),
            ]
        )
        for offset, reaction in enumerate(high_reactions):
            cell = (reaction_index[reaction], column)
            positive_addresses[cell] = _BinaryAddress(positive_name, offset, True)
            negative_addresses[cell] = _BinaryAddress(negative_name, offset, True)

    additional_positive_cells: list[tuple[int, int]] = []
    additional_negative_cells: list[tuple[int, int]] = []
    for row in range(n_reactions):
        for column in range(n_conditions):
            cell = (row, column)
            if cell in positive_addresses:
                continue
            if upper[row, column] >= eps:
                additional_positive_cells.append(cell)
            if lower[row, column] <= -eps:
                additional_negative_cells.append(cell)

    backend = problem.backend
    positive_name = "activity_classifier_positive" if additional_positive_cells else None
    negative_name = "activity_classifier_negative" if additional_negative_cells else None
    positive_variable = (
        backend.Variable(
            name=positive_name,
            shape=(len(additional_positive_cells),),
            vartype=VarType.BINARY,
        )
        if additional_positive_cells
        else None
    )
    negative_variable = (
        backend.Variable(
            name=negative_name,
            shape=(len(additional_negative_cells),),
            vartype=VarType.BINARY,
        )
        if additional_negative_cells
        else None
    )

    for offset, cell in enumerate(additional_positive_cells):
        positive_addresses[cell] = _BinaryAddress(positive_name, offset, False)
    for offset, cell in enumerate(additional_negative_cells):
        negative_addresses[cell] = _BinaryAddress(negative_name, offset, False)

    new_constraints = list(native_classifier_constraints)
    for direction, cells, variable in (
        (1, additional_positive_cells, positive_variable),
        (-1, additional_negative_cells, negative_variable),
    ):
        if not cells:
            continue
        selected_flow = _select_cells(flow, cells, n_conditions)
        row_indices = np.asarray([row for row, _ in cells], dtype=int)
        column_indices = np.asarray([column for _, column in cells], dtype=int)
        cell_lower = lower[row_indices, column_indices]
        cell_upper = upper[row_indices, column_indices]
        if direction > 0:
            # p=1 implies v>=eps; p=0 implies v<=eps. These two
            # implications partition the original bounded flux interval.
            new_constraints.extend(
                [
                    selected_flow >= cell_lower + variable.multiply(eps - cell_lower),
                    selected_flow
                    <= np.minimum(cell_upper, eps)
                    + variable.multiply(np.maximum(cell_upper - eps, 0.0)),
                ]
            )
        else:
            # n=1 implies v<=-eps; n=0 implies v>=-eps.
            new_constraints.extend(
                [
                    selected_flow <= cell_upper - variable.multiply(cell_upper + eps),
                    selected_flow
                    >= np.maximum(cell_lower, -eps)
                    + variable.multiply(np.minimum(cell_lower + eps, 0.0)),
                ]
            )

    # Enforce p+n<=1 for every newly classified reversible reaction. Native
    # high-expression indicators already include the same constraint.
    positive_offsets = {cell: offset for offset, cell in enumerate(additional_positive_cells)}
    negative_offsets = {cell: offset for offset, cell in enumerate(additional_negative_cells)}
    reversible_cells = tuple(cell for cell in additional_positive_cells if cell in negative_offsets)
    if reversible_cells:
        pos = np.asarray([positive_offsets[cell] for cell in reversible_cells], dtype=int)
        neg = np.asarray([negative_offsets[cell] for cell in reversible_cells], dtype=int)
        new_constraints.append(positive_variable[pos] + negative_variable[neg] <= 1)

    if new_constraints:
        problem.add_constraints(new_constraints)

    positive_expression = _scatter_direction(
        problem,
        backend,
        n_reactions * n_conditions,
        native_positive_by_condition,
        positive_variable,
        additional_positive_cells,
        condition_index,
        reactions,
        n_conditions,
        positive=True,
    )
    negative_expression = _scatter_direction(
        problem,
        backend,
        n_reactions * n_conditions,
        native_positive_by_condition,
        negative_variable,
        additional_negative_cells,
        condition_index,
        reactions,
        n_conditions,
        positive=False,
    )
    activity_flat = positive_expression + negative_expression
    activity_expression = activity_flat.reshape((n_reactions, n_conditions))

    structural_selector = problem.expr.edge_has_flux
    selector_flat = structural_selector.reshape((n_reactions * n_conditions,))
    problem.add_constraints(activity_flat <= selector_flat)
    problem.register("activity", activity_expression)

    return ActivityClassifier(
        reactions=reactions,
        conditions=conditions,
        eps=float(eps),
        expression=activity_expression,
        positive_addresses=positive_addresses,
        negative_addresses=negative_addresses,
        additional_positive_cells=tuple(additional_positive_cells),
        additional_negative_cells=tuple(additional_negative_cells),
        positive_name=positive_name,
        negative_name=negative_name,
    )


def _scatter_direction(
    problem,
    backend,
    size,
    native_by_condition,
    additional_variable,
    additional_cells,
    condition_index,
    reactions,
    n_conditions,
    *,
    positive,
):
    expressions = []
    for condition, (positive_name, negative_name, rows) in native_by_condition.items():
        name = positive_name if positive else negative_name
        if not name:
            continue
        column = condition_index[condition]
        values = np.asarray(rows, dtype=int)
        columns = np.arange(len(rows), dtype=int)
        scatter = coo_matrix(
            (np.ones(len(rows)), (values + len(reactions) * column, columns)),
            shape=(size, len(rows)),
        ).tocsr()
        expressions.append(backend.Constant(scatter) @ problem.symbols[name])

    if additional_cells:
        rows = np.asarray(
            [row + len(reactions) * column for row, column in additional_cells], dtype=int
        )
        columns = np.arange(len(additional_cells), dtype=int)
        scatter = coo_matrix(
            (np.ones(len(additional_cells)), (rows, columns)),
            shape=(size, len(additional_cells)),
        ).tocsr()
        expressions.append(backend.Constant(scatter) @ additional_variable)

    if not expressions:
        return backend.Constant(np.zeros(size))
    result = expressions[0]
    for expression in expressions[1:]:
        result = result + expression
    return result


def _select_cells(flow, cells, n_conditions):
    rows = np.asarray([row for row, _ in cells], dtype=int)
    columns = np.asarray([column for _, column in cells], dtype=int)
    if len(flow.shape) == 1:
        return flow[rows]
    return flow[(rows, columns)]


def _bound_matrices(lower, upper, n_reactions, n_conditions):
    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    if n_conditions == 1 and lower.shape == (n_reactions,):
        lower = lower[:, None]
        upper = upper[:, None]
    expected = (n_reactions, n_conditions)
    if lower.shape != expected or upper.shape != expected:
        raise RuntimeError(f"Flow bounds have shapes {lower.shape}/{upper.shape}; expected {expected}.")
    return lower, upper


def _assign_binary_start(problem, name, values):
    symbol = problem.symbols.get(name)
    if symbol is None or getattr(symbol, "_vartype", None) != VarType.BINARY:
        raise RuntimeError(f"Cannot seed missing activity binary {name!r}.")
    vector = np.asarray(values, dtype=float)
    if vector.size != int(np.prod(symbol.shape, dtype=int)):
        raise RuntimeError(f"Activity warm start for {name!r} has the wrong number of entries.")
    symbol.value = np.rint(vector).astype(bool).reshape(symbol.shape)


def _finite_array(value, label):
    if value is None:
        raise RuntimeError(f"Solved {label} has no numeric value.")
    result = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(result)):
        raise RuntimeError(f"Solved {label} contains non-finite values.")
    return result


def _integrality_tolerance():
    return 1e-5


__all__ = ["ActivityClassifier", "build_activity_classifier"]
