from __future__ import annotations

import json
import math
import os
import textwrap
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", str(Path("/tmp/matplotlib").resolve()))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd


ROOT = Path("outputs")
CSV_DIR = ROOT / "csvs"
RICH_DIR = ROOT / "figures"

SOURCE_COLORS = {
    "chembl": "#31688e",
    "drugcentral": "#f28e2b",
    "foodb": "#59a14f",
    "hmdb": "#e15759",
    "macdb": "#8f63b8",
    "metatlas": "#9c755f",
    "pfocr": "#d66ab0",
    "rhea": "#4e79a7",
    "stitch": "#7f7f7f",
    "chebi": "#b07aa1",
    "swisslipids": "#76b7b2",
}

CATEGORY_COLORS = {
    "Center": "#cc3d3d",
    "ID": "#2f78b7",
    "Structure": "#21b7c6",
    "Food": "#44a35f",
    "Metabolic reaction": "#9a6cc8",
    "Class": "#f28e2b",
}

TISSUE_COLORS = {
    "Blood": "#c6403d",
    "Serum": "#f28e2b",
    "Plasma": "#8c564b",
    "Tissue": "#59a14f",
    "Stool": "#4e79a7",
    "Urine": "#9467bd",
    "Sweat": "#7f7f7f",
    "Cells": "#e377c2",
}

TYPE_COLORS = {
    "Chemical:OM:0037": "#2f78b7",
    "Protein:MI:0326": "#e15759",
    "Cv Term:OM:0012": "#8f63b8",
    "Food:OM:0020": "#44a35f",
    "Reaction:OM:0015": "#f28e2b",
    "Transport:OM:0035": "#9c755f",
}


def set_style() -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": "#fbfaf6",
            "axes.facecolor": "#fbfaf6",
            "axes.edgecolor": "#333333",
            "axes.labelcolor": "#111111",
            "axes.titleweight": "bold",
            "axes.titlesize": 15,
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "grid.color": "#d8d2c3",
            "grid.linewidth": 0.8,
            "legend.frameon": True,
            "savefig.facecolor": "#fbfaf6",
        }
    )


def split_pipe(value: Any) -> list[str]:
    if pd.isna(value):
        return []
    parts = [part.strip() for part in str(value).split("|")]
    return [part for part in parts if part and part.lower() != "nan"]


def shorten(value: Any, width: int = 42) -> str:
    text = "" if pd.isna(value) else str(value)
    return textwrap.shorten(text, width=width, placeholder="...")


def wrap(value: Any, width: int = 26) -> str:
    text = "" if pd.isna(value) else str(value)
    return "\n".join(textwrap.wrap(text, width=width, break_long_words=False)) or text


def normalize_label(value: Any) -> str:
    text = "" if pd.isna(value) else str(value)
    return " ".join(text.casefold().split())


def preferred_label_score(value: Any) -> tuple[int, int, str]:
    text = "" if pd.isna(value) else str(value)
    tail = text.split(":", 1)[-1].strip()
    all_upper_penalty = int(len(tail) > 3 and tail.upper() == tail)
    all_lower_penalty = int(len(tail) > 3 and tail.lower() == tail)
    return (all_upper_penalty + all_lower_penalty, len(text), text)


def source_color(source: str) -> str:
    return SOURCE_COLORS.get(str(source), "#6f6f6f")


def load_tables(csv_dir: Path = CSV_DIR) -> dict[str, pd.DataFrame]:
    tables: dict[str, pd.DataFrame] = {}
    for csv_path in sorted(csv_dir.glob("*.csv")):
        try:
            tables[csv_path.stem] = pd.read_csv(csv_path)
        except pd.errors.EmptyDataError:
            tables[csv_path.stem] = pd.DataFrame()
    return tables


def filter_resolved_for_plots(df: pd.DataFrame) -> pd.DataFrame:
    for col in ("matched_seed_resolution_status", "seed_resolution_status", "query_seed_resolution_status"):
        if col in df.columns:
            return df[df[col].fillna("").eq("resolved")].copy()
    return df


def resolved_tables_for_plots(tables: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    return {name: filter_resolved_for_plots(df) for name, df in tables.items()}


def savefig(fig: plt.Figure, out_path: Path) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=240, bbox_inches="tight")
    plt.close(fig)
    return out_path


def explode_relation_sources(direct_relations: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for _, row in direct_relations.iterrows():
        for source in split_pipe(row.get("sources")):
            rows.append(
                {
                    "relation_pk": row.get("relation_pk"),
                    "interaction_type": row.get("interaction_type") or row.get("predicate"),
                    "relation_category": row.get("relation_category"),
                    "source": source,
                    "partner_entity_type": row.get("partner_entity_type") or row.get("other_entity_types"),
                    "partner_label": row.get("partner_label") or row.get("other_entity_labels"),
                }
            )
    return pd.DataFrame(rows)


def make_extra_mindmap(tables: dict[str, pd.DataFrame], out_dir: Path) -> Path:
    data = tables["fig_mindmap_data"].copy()
    nodes = data[data["record_type"].eq("node")].rename(columns={"evidence": "evidence"}).copy()
    edges = data[data["record_type"].eq("edge")].rename(
        columns={"source_node": "source", "target_node": "target", "weight": "weight"}
    ).copy()
    if nodes.empty or edges.empty:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "No resolved mindmap data available", ha="center", va="center")
        ax.axis("off")
        return savefig(fig, out_dir / "extra_evidence_atlas_mindmap.png")
    nodes = nodes[~(nodes["category"].eq("ID") & nodes["label"].fillna("").str.startswith("Name:"))].copy()
    nodes["_norm_label"] = nodes["label"].map(normalize_label)
    nodes["_label_score"] = nodes["label"].map(preferred_label_score)
    nodes = (
        nodes.sort_values(["category", "_norm_label", "_label_score", "evidence"], ascending=[True, True, True, False])
        .drop_duplicates(["category", "_norm_label"], keep="first")
        .drop(columns=["_norm_label", "_label_score"])
    )
    kept_node_ids = set(nodes["node_id"])
    edges = edges[edges["source"].isin(kept_node_ids) & edges["target"].isin(kept_node_ids)].copy()
    center = nodes.loc[nodes["category"].eq("Center"), "node_id"].iloc[0]

    graph = nx.Graph()
    for _, row in nodes.iterrows():
        graph.add_node(
            row["node_id"],
            label=row["label"],
            category=row["category"],
            evidence=float(row.get("evidence", 1) or 1),
        )
    for _, row in edges.iterrows():
        graph.add_edge(row["source"], row["target"], category=row["category"], weight=float(row.get("weight", 1) or 1))

    pos = {center: (0.0, 0.0)}
    sector_specs = {
        "Structure": (65, 135, 4.0),
        "ID": (150, 235, 4.6),
        "Class": (230, 305, 3.9),
        "Metabolic reaction": (255, 340, 5.0),
        "Food": (-15, 45, 4.5),
    }
    for category, (start, end, radius) in sector_specs.items():
        category_nodes = nodes[nodes["category"].eq(category)].copy()
        if category_nodes.empty:
            continue
        category_nodes = category_nodes.sort_values(["evidence", "label"], ascending=[False, True])
        angles = [0.5 * (start + end)] if len(category_nodes) == 1 else np.linspace(start, end, len(category_nodes))
        for node_id, angle in zip(category_nodes["node_id"], angles):
            rad = math.radians(float(angle))
            jitter = 0.13 * math.sin(len(str(node_id)))
            pos[node_id] = ((radius + jitter) * math.cos(rad), (radius + jitter) * math.sin(rad))

    fig, ax = plt.subplots(figsize=(17, 12))
    ax.axis("off")

    for category in ["Center", "ID", "Structure", "Food", "Metabolic reaction", "Class"]:
        nodelist = [node for node, attrs in graph.nodes(data=True) if attrs["category"] == category]
        if not nodelist:
            continue
        sizes = [780 if category == "Center" else 240 + min(graph.nodes[node]["evidence"], 8) * 95 for node in nodelist]
        nx.draw_networkx_nodes(
            graph,
            pos,
            nodelist=nodelist,
            node_color=CATEGORY_COLORS.get(category, "#777777"),
            node_size=sizes,
            alpha=0.92,
            linewidths=1.1,
            edgecolors="#fbfaf6",
            ax=ax,
        )

    edge_colors = [CATEGORY_COLORS.get(graph.edges[edge].get("category"), "#999999") for edge in graph.edges]
    edge_widths = [0.7 + min(float(graph.edges[edge].get("weight") or 1), 8) * 0.25 for edge in graph.edges]
    nx.draw_networkx_edges(graph, pos, edge_color=edge_colors, width=edge_widths, alpha=0.35, ax=ax)

    labels = {
        node: wrap(attrs["label"], 42) if attrs["category"] == "Metabolic reaction" else shorten(attrs["label"], 48)
        for node, attrs in graph.nodes(data=True)
    }
    nx.draw_networkx_labels(graph, pos, labels=labels, font_size=7.3, font_color="#151515", ax=ax)

    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=CATEGORY_COLORS[key], label=key, markersize=9)
        for key in ["ID", "Structure", "Food", "Metabolic reaction", "Class"]
        if key in set(nodes["category"])
    ]
    ax.legend(handles=handles, title="Node type", loc="upper left", bbox_to_anchor=(1.01, 1.0))
    return savefig(fig, out_dir / "extra_evidence_atlas_mindmap.png")


def make_fig3_panel_D_cancer_stacked(tables: dict[str, pd.DataFrame], out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = tables["fig3_panel_D_data"].copy()
    if summary.empty:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "No resolved cancer association data available", ha="center", va="center")
        ax.axis("off")
        return savefig(fig, out_dir / "fig3_panel_D_cancer_associations_sample_type.png")
    order = summary.groupby("disease_type")["evidence_count"].sum().sort_values(ascending=False).index
    pivot = summary.pivot_table(index="disease_type", columns="sample_type", values="evidence_count", aggfunc="sum", fill_value=0)
    tissue_order = [col for col in pivot.sum(axis=0).sort_values(ascending=False).index]
    colors = [TISSUE_COLORS.get(tissue, "#6f6f6f") for tissue in tissue_order]
    pivot = pivot.reindex(index=order, columns=tissue_order)

    fig, ax = plt.subplots(figsize=(13.5, 6.8))
    pivot.plot(kind="bar", stacked=True, color=colors, ax=ax, width=0.72)
    ax.set_xlabel("Disease type")
    ax.set_ylabel("Evidence count (unique disease/tissue/study units)")
    ax.grid(axis="y", alpha=0.55)
    ax.legend(title="Sample Type", bbox_to_anchor=(1.02, 1), loc="upper left")
    ax.set_xticklabels([wrap(label, 18) for label in pivot.index], rotation=35, ha="right")
    for i, total in enumerate(pivot.sum(axis=1)):
        ax.text(i, total + 0.35, str(int(total)), ha="center", va="bottom", fontsize=9, color="#333333")
    return savefig(fig, out_dir / "fig3_panel_D_cancer_associations_sample_type.png")


def make_fig3_panel_C_interaction_type_barplot(tables: dict[str, pd.DataFrame], out_dir: Path) -> Path:
    counts = tables["fig3_panel_C_data"].rename(columns={"relation_count": "count"}).copy()
    if counts.empty:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "No resolved interaction data available", ha="center", va="center")
        ax.axis("off")
        return savefig(fig, out_dir / "fig3_panel_C_interaction_type_barplot.png")
    order = counts.groupby("interaction_type")["count"].sum().sort_values(ascending=False).index
    source_order = counts.groupby("source")["count"].sum().sort_values(ascending=False).index
    pivot = counts.pivot_table(index="interaction_type", columns="source", values="count", aggfunc="sum", fill_value=0)
    pivot = pivot.reindex(index=order, columns=source_order)

    fig, ax = plt.subplots(figsize=(12.5, 7))
    colors = [source_color(source) for source in pivot.columns]
    pivot.plot(kind="bar", stacked=True, color=colors, width=0.68, ax=ax)
    ax.set_xlabel("Interaction type")
    ax.set_ylabel("Relation count")
    ax.grid(axis="y", alpha=0.45)
    ax.set_xticklabels([label.replace("_", "\n") for label in pivot.index], rotation=0)
    ax.legend(title="Source", bbox_to_anchor=(1.02, 1), loc="upper left")
    for i, total in enumerate(pivot.sum(axis=1)):
        ax.text(i, total + 0.35, str(int(total)), ha="center", va="bottom", fontsize=9)
    return savefig(fig, out_dir / "fig3_panel_C_interaction_type_barplot.png")


def make_extra_target_type_by_source(tables: dict[str, pd.DataFrame], out_dir: Path) -> Path:
    rel_sources = explode_relation_sources(tables["relations"])
    rows: list[dict[str, Any]] = []
    for _, row in rel_sources.iterrows():
        for entity_type in split_pipe(row.get("partner_entity_type")):
            rows.append({"source": row["source"], "target_type": entity_type, "relation_pk": row["relation_pk"]})
    df = pd.DataFrame(rows)
    counts = df.groupby(["source", "target_type"])["relation_pk"].nunique().reset_index(name="count")
    source_order = counts.groupby("source")["count"].sum().sort_values(ascending=False).index
    type_order = counts.groupby("target_type")["count"].sum().sort_values(ascending=False).index
    pivot = counts.pivot_table(index="source", columns="target_type", values="count", aggfunc="sum", fill_value=0)
    pivot = pivot.reindex(index=source_order, columns=type_order)

    fig, ax = plt.subplots(figsize=(11, 6.5))
    colors = [TYPE_COLORS.get(col, "#777777") for col in pivot.columns]
    pivot.plot(kind="barh", stacked=True, color=colors, ax=ax, width=0.72)
    ax.set_xlabel("Direct relation count")
    ax.set_ylabel("Source")
    ax.grid(axis="x", alpha=0.42)
    ax.legend(title="Target type", bbox_to_anchor=(1.02, 1), loc="upper left")
    return savefig(fig, out_dir / "extra_target_type_composition_by_source.png")


def make_all_visualizations(root: Path = ROOT, out_dir: Path | None = None) -> pd.DataFrame:
    set_style()
    csv_dir = root / "csvs"
    out_dir = out_dir or (root / "figures")
    out_dir.mkdir(parents=True, exist_ok=True)
    tables = resolved_tables_for_plots(load_tables(csv_dir))

    figure_funcs = [
        ("main_fig3_panel", "C", "Interaction type by source", make_fig3_panel_C_interaction_type_barplot),
        ("main_fig3_panel", "D", "Cancer associations by sample type", make_fig3_panel_D_cancer_stacked),
        ("extra", "", "Evidence atlas mindmap", make_extra_mindmap),
        ("extra", "", "Target type by source", make_extra_target_type_by_source),
    ]

    rows: list[dict[str, Any]] = []
    for figure_type, panel, title, func in figure_funcs:
        path = func(tables, out_dir)
        rows.append({"figure_type": figure_type, "panel": panel, "title": title, "path": str(path)})

    index = pd.DataFrame(rows)
    index.to_csv(out_dir / "figure_index.csv", index=False)
    (out_dir / "figure_index.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
    return index


if __name__ == "__main__":
    index = make_all_visualizations()
    print(index.to_string(index=False))
