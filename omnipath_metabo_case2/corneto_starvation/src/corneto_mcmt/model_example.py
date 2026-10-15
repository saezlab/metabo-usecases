"""Small plotting helpers for the beginner-facing model example notebook.

These helpers only format simulated observations and solved CORNETO results.
The model, score preprocessing, and temporal iMAT fit remain in their own
modules and are called visibly from the notebook.
"""

from __future__ import annotations

from io import BytesIO
import shutil
import subprocess
from typing import Any

from IPython.display import Image
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage


_TIME_BOUNDARY = 3.5
_ACTIVE_POSITIVE = "#176b87"
_ACTIVE_NEGATIVE = "#b5443f"
_INACTIVE = "#d4d9de"


def _figure_png(fig: plt.Figure) -> Image:
    """Return a native notebook PNG that GitHub can render."""
    buffer = BytesIO()
    fig.savefig(buffer, format="png", bbox_inches="tight", dpi=140)
    plt.close(fig)
    return Image(data=buffer.getvalue(), format="png")


def _graph_html(graph: Any, **plot_options: Any) -> str:
    """Render a CORNETO graph through its Graph.plot Graphviz backend."""
    dot = graph.plot(renderer="graphviz", **plot_options)
    # Some conda Graphviz installs need their plugin registry generated once.
    # CORNETO returns its HTML fallback when that registry is missing.
    if not callable(getattr(dot, "pipe", None)):
        executable = shutil.which("dot")
        if executable is not None:
            subprocess.run([executable, "-c"], check=True, capture_output=True, text=True)
        dot = graph.plot(renderer="graphviz", **plot_options)
    if not callable(getattr(dot, "pipe", None)):
        raise RuntimeError("CORNETO Graph.plot did not return a Graphviz graph")
    svg = dot.pipe(format="svg").decode("utf-8")
    svg = svg.replace("<svg ", '<svg style="max-width:100%;height:auto;display:block" ', 1)
    return f'<div style="width:100%; overflow-x:auto">{svg}</div>'


def _reaction_label_attributes(graph: Any) -> tuple[dict[str, dict[str, str]], dict[int, str]]:
    """Label CORNETO's generated reaction nodes or direct reaction edges."""
    reaction_nodes: dict[str, dict[str, str]] = {}
    direct_edge_labels: dict[int, str] = {}
    for edge_index in range(graph.num_edges):
        sources, targets = graph.get_edge(edge_index)
        reaction_id = str(graph.get_attr_edge(edge_index).get("id", f"R_{edge_index}"))
        if len(sources) > 1 or len(targets) > 1:
            reaction_nodes[f"e_{edge_index}_center"] = {
                "label": reaction_id,
                "shape": "box",
                "fixedsize": "false",
                "fontsize": "9",
                "fontname": "Helvetica",
                "color": "#89939b",
                "fontcolor": "#44515b",
                "margin": "0.025,0.01",
            }
        else:
            direct_edge_labels[edge_index] = reaction_id
    return reaction_nodes, direct_edge_labels


def plot_full_model(graph: Any, layout: str = "dot") -> str:
    """Show the full CORNETO metabolic graph as an SVG overview.

    layout="dot" (Graphviz's hierarchical layout, the original default) does
    not scale -- confirmed hung >15h on a 4559-reaction core network.
    layout="neato" (force-directed) renders the same network in ~9 minutes.
    """
    return _graph_html(
        graph,
        layout=layout,
        graph_attr={"rankdir": "LR", "pad": "0.15", "nodesep": "0.18", "ranksep": "0.35"},
        node_attr={"fontname": "Helvetica", "fontsize": "8", "margin": "0.04,0.02"},
        edge_attr={"fontname": "Helvetica", "fontsize": "6", "arrowsize": "0.55"},
    )


def _cluster_observed_rows(values: np.ndarray) -> np.ndarray:
    """Cluster rows from observed values only, retaining sample order."""
    values = np.asarray(values, dtype=float)
    fallback = float(np.nanmedian(values)) if np.isfinite(values).any() else 0.0
    row_means = np.array(
        [row[np.isfinite(row)].mean() if np.isfinite(row).any() else fallback for row in values]
    )
    row_means = np.where(np.isfinite(row_means), row_means, fallback)
    filled = np.where(np.isfinite(values), values, row_means[:, None])
    if filled.shape[0] < 2:
        return np.arange(filled.shape[0])
    return leaves_list(linkage(filled, method="average", metric="euclidean"))


def plot_input_heatmaps(dataset: Any) -> tuple[Image, Image]:
    """Show proteomics and GPR-derived scores, clustering observed rows."""
    labels = list(dataset.conditions)
    measurements = np.asarray(dataset.proteomics.state_expression, dtype=float)
    scores = np.asarray(dataset.proteomics.reaction_score_matrix, dtype=float)
    gene_order = _cluster_observed_rows(measurements)
    reaction_order = _cluster_observed_rows(scores)

    measurement_cmap = plt.get_cmap("viridis").copy()
    measurement_cmap.set_bad("#a6a6a6")
    score_cmap = plt.get_cmap("RdBu_r").copy()
    score_cmap.set_bad("#a6a6a6")

    figures: list[Image] = []
    fig, ax = plt.subplots(figsize=(10.5, max(5.5, 0.19 * len(dataset.genes))), dpi=130)
    image = ax.imshow(
        np.ma.masked_invalid(measurements[gene_order]),
        aspect="auto",
        interpolation="nearest",
        cmap=measurement_cmap,
    )
    ax.set_title("Proteomics measurements", loc="left", fontsize=13, weight="bold")
    ax.set_yticks(np.arange(len(gene_order)))
    ax.set_yticklabels([dataset.genes[index] for index in gene_order], fontsize=7)
    ax.set_xticks(np.arange(len(labels)), labels)
    ax.set_xlabel("Timepoint (in order)")
    ax.set_ylabel("Genes, clustered by measured values")
    ax.axvline(_TIME_BOUNDARY, color="#222222", linestyle="--", linewidth=1.2)
    fig.colorbar(image, ax=ax, pad=0.02, label="Mean protein measurement; gray = missing")
    fig.text(0.5, 0.005, "Dashed line marks the t3 → t4 boundary.", ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    figures.append(_figure_png(fig))

    fig, ax = plt.subplots(figsize=(10.5, max(5.5, 0.19 * len(dataset.reactions))), dpi=130)
    image = ax.imshow(
        scores[reaction_order],
        aspect="auto",
        interpolation="nearest",
        cmap=score_cmap,
        vmin=-1,
        vmax=1,
    )
    ax.set_title("GPR-derived reaction scores", loc="left", fontsize=13, weight="bold")
    ax.set_yticks(np.arange(len(reaction_order)))
    ax.set_yticklabels([dataset.reactions[index] for index in reaction_order], fontsize=7)
    ax.set_xticks(np.arange(len(labels)), labels)
    ax.set_xlabel("Timepoint (in order)")
    ax.set_ylabel("Reactions, clustered by observed scores")
    ax.axvline(_TIME_BOUNDARY, color="#222222", linestyle="--", linewidth=1.2)
    fig.colorbar(image, ax=ax, pad=0.02, ticks=[-1, 0, 1], label="Reaction score (−1, 0, +1)")
    fig.text(0.5, 0.005, "Dashed line marks the t3 → t4 boundary.", ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    figures.append(_figure_png(fig))

    return figures[0], figures[1]


def score_count_table(dataset: Any) -> pd.DataFrame:
    """Count high, neutral/unavailable, and low reaction scores per timepoint."""
    matrix = np.asarray(dataset.proteomics.reaction_score_matrix)
    return pd.DataFrame(
        {
            "high (+1)": np.count_nonzero(matrix == 1, axis=0),
            "neutral or unavailable (0)": np.count_nonzero(matrix == 0, axis=0),
            "low (−1)": np.count_nonzero(matrix == -1, axis=0),
        },
        index=pd.Index(dataset.conditions, name="timepoint"),
    )


def plot_inference_gallery(graph: Any, result: Any, layout: str = "dot") -> str:
    """Show eight solved CORNETO views on one fixed inferred-active subgraph.

    The edge set is the union of reactions active at any timepoint. Each
    panel reuses that graph and edge ordering. Color uses ``result.activity``;
    ``result.selected`` is intentionally not used because it records a
    structural selector and may include zero-flux reactions.

    layout="dot" does not scale to large networks (see plot_full_model).
    Note this renders the SVG layout once per condition (8x), so the time
    cost of a slow layout multiplies -- a 9-minute single-render layout
    could mean ~70+ minutes here if the active-reaction subset stays large.
    """
    reactions = tuple(
        str(graph.get_attr_edge(index).get("id", f"R_{index}"))
        for index in range(graph.num_edges)
    )
    if reactions != tuple(result.reactions):
        raise ValueError("The result reaction order does not match the displayed CORNETO graph")
    if result.activity is None:
        raise ValueError("The result must include threshold-based reaction activity")

    activity = np.asarray(result.activity, dtype=bool)
    flux = np.asarray(result.flux, dtype=float)
    if activity.shape != flux.shape or activity.shape != (len(reactions), len(result.conditions)):
        raise ValueError("Flux and activity must be reaction-by-timepoint arrays")

    selected_indices = np.flatnonzero(activity.any(axis=1))
    if selected_indices.size == 0:
        raise ValueError("No inferred active reactions are available to plot")
    context_graph = graph.edge_subgraph(selected_indices)
    maximum_flux = max(float(np.max(np.abs(flux[selected_indices]))), 1e-12)
    reaction_nodes, direct_edge_labels = _reaction_label_attributes(context_graph)

    graph_attrs = {"rankdir": "LR", "pad": "0.08", "nodesep": "0.14", "ranksep": "0.28"}
    node_attrs = {"fontname": "Helvetica", "fontsize": "8", "margin": "0.03,0.015"}
    edge_attrs = {"fontname": "Helvetica", "fontsize": "8", "arrowsize": "0.5"}
    panels: list[str] = []
    for time_index, condition in enumerate(result.conditions):
        custom_edge_attr: dict[int, dict[str, str]] = {}
        for local_index, reaction_index in enumerate(selected_indices):
            value = float(flux[reaction_index, time_index])
            is_active = bool(activity[reaction_index, time_index])
            width = 0.75 + 3.0 * np.sqrt(abs(value) / maximum_flux) if is_active else 0.65
            custom_edge_attr[int(local_index)] = {
                "color": (
                    _ACTIVE_POSITIVE if value > 0 else _ACTIVE_NEGATIVE
                ) if is_active else _INACTIVE,
                "fontcolor": "#52606a" if is_active else "#a5abb1",
                "penwidth": f"{width:.2f}",
            }
            if local_index in direct_edge_labels:
                custom_edge_attr[int(local_index)]["label"] = direct_edge_labels[int(local_index)]
        svg = _graph_html(
            context_graph,
            layout=layout,
            graph_attr=graph_attrs,
            node_attr=node_attrs,
            edge_attr=edge_attrs,
            custom_edge_attr=custom_edge_attr,
            custom_vertex_attr=reaction_nodes,
        )
        panels.append(
            f'<section style="min-width:0;border:1px solid #d6dbe0;padding:6px">'
            f'<h4 style="margin:2px 0 4px;text-align:center">{condition}</h4>{svg}</section>'
        )

    legend = (
        '<div style="display:flex;gap:18px;flex-wrap:wrap;margin:8px 0;font-size:0.9em">'
        f'<span><b style="color:{_ACTIVE_POSITIVE}">━━</b> active, positive flux</span>'
        f'<span><b style="color:{_ACTIVE_NEGATIVE}">━━</b> active, reverse flux</span>'
        f'<span><b style="color:{_INACTIVE}">━━</b> inactive at this timepoint</span>'
        '<span>Active edge width follows absolute flux.</span></div>'
    )
    return (
        f'<p>All eight panels use the same {len(selected_indices)}-reaction CORNETO subgraph. '
        "Only the timepoint's flux and activity styles change.</p>"
        + legend
        + '<div style="display:grid;grid-template-columns:minmax(0,1fr);gap:12px">'
        + "".join(panels)
        + "</div>"
    )


def plot_activity_overview(result: Any) -> Image:
    """Show inferred threshold activity before any simulated truth is compared."""
    activity = np.asarray(result.activity, dtype=bool)
    fig, ax = plt.subplots(figsize=(10, 9), dpi=135)
    image = ax.imshow(activity, aspect="auto", interpolation="nearest", cmap="Blues", vmin=0, vmax=1)
    ax.set_title("Inferred reaction activity", loc="left", fontsize=13, weight="bold")
    ax.set_xticks(np.arange(len(result.conditions)), result.conditions)
    ax.set_yticks(np.arange(len(result.reactions)))
    ax.set_yticklabels(result.reactions, fontsize=7)
    ax.set_xlabel("Timepoint (in order)")
    ax.set_ylabel("Reaction")
    ax.axvline(_TIME_BOUNDARY, color="#a73935", linestyle="--", linewidth=1.6)
    ax.text(3.48, -1.5, "t3 → t4", ha="right", va="bottom", color="#8f302c", fontsize=9)
    colorbar = fig.colorbar(image, ax=ax, pad=0.02, ticks=[0, 1])
    colorbar.ax.set_yticklabels(["inactive", "active"])
    fig.tight_layout()
    return _figure_png(fig)


def evaluation_table(evaluation: Any) -> pd.DataFrame:
    """Summarize network accuracy and each transition class for the guide."""
    rows = [
        {
            "Measure": "Active-network F1",
            "Value": f"{evaluation.active_f1:.3f}",
            "Meaning": "Active/inactive agreement across all reactions and timepoints.",
        },
        {
            "Measure": "Switch detection F1",
            "Value": f"{evaluation.switch_f1:.3f}",
            "Meaning": "Any on/off change at the right reaction and timepoint, regardless of direction.",
        },
        {
            "Measure": "Switch precision / recall",
            "Value": f"{evaluation.switch_precision:.3f} / {evaluation.switch_recall:.3f}",
            "Meaning": "Fraction of predicted switches that are real / fraction of real switches found.",
        },
    ]
    transition_classes = (
        ("Off→off F1", "Inactive at both timepoints."),
        ("Off→on F1", "Becomes active."),
        ("On→off F1", "Becomes inactive."),
        ("On→on F1", "Active at both timepoints."),
    )
    confusion = np.asarray(evaluation.transition_confusion, dtype=int)
    defined = confusion.sum(axis=0) + confusion.sum(axis=1) > 0
    for index, (name, meaning) in enumerate(transition_classes):
        rows.append(
            {
                "Measure": name,
                "Value": f"{evaluation.transition_class_f1[index]:.3f}" if defined[index] else "—",
                "Meaning": meaning,
            }
        )
    rows.extend(
        [
            {
                "Measure": "Four-class transition macro F1",
                "Value": f"{evaluation.transition_macro_f1:.3f}",
                "Meaning": "Mean of the defined transition-class F1 values.",
            },
            {
                "Measure": "Stationary false switches",
                "Value": str(evaluation.stationary_false_switches),
                "Meaning": "Inferred switches where the simulated regime stays the same.",
            },
        ]
    )
    return pd.DataFrame(rows)
