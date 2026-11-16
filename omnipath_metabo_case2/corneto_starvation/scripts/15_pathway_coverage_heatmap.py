"""Pathway-level summary of Script 14's metabolomics-anchored network --
the full-graph SVG from that script (7,380 nodes, 23,602 edges) is not
practically viewable, so this aggregates to Mouse-GEM's own 151 subsystem
annotations (SBML groups extension; every one of the 12,987 reactions
belongs to exactly one subsystem -- a clean partition, not an overlapping
tagging) instead.

For each subsystem: what fraction of its reactions were selected, and among
those selected, how many are explained by which evidence (touching a
variable plasma/tissue metabolite vs. proteomics support) -- so pathways
that are well-covered by real evidence are visually distinct from ones that
got swept in mostly by proteomics alone or barely survived at all.

Outputs (corneto_starvation/result/metabolomics_anchored_network/):
  pathway_coverage.csv            per-subsystem reaction/evidence counts
  pathway_coverage_heatmap.png    all 151 subsystems, sorted by % selected
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd

PKG_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[2]
MOUSE_GEM_XML = PKG_ROOT / "assets" / "models" / "Mouse-GEM-v1.8.0.xml"
OUTDIR = PKG_ROOT / "result" / "metabolomics_anchored_network"
SELECTED_REACTIONS_PATH = OUTDIR / "selected_reactions.csv"

# dataviz skill reference palette, sequential blue ramp (steps 100/300/500/700).
BLUE_RAMP = LinearSegmentedColormap.from_list(
    "blue_sequential", ["#cde2fb", "#6da7ec", "#256abf", "#0d366b"],
)
MIN_SUBSYSTEM_SIZE = 3  # drop subsystems with fewer reactions than this (too noisy)


def parse_subsystems(xml_path: Path) -> dict[str, list[str]]:
    content = xml_path.read_text()
    groups = re.findall(
        r'<groups:group[^>]*groups:name="([^"]*)"[^>]*>(.*?)</groups:group>', content, re.S,
    )
    subsystems: dict[str, list[str]] = {}
    for name, body in groups:
        subsystems[name] = re.findall(r'groups:idRef="(MAR\d+)"', body)
    return subsystems


def main() -> None:
    subsystems = parse_subsystems(MOUSE_GEM_XML)
    total_reactions = sum(len(v) for v in subsystems.values())
    print(f"Mouse-GEM subsystems: {len(subsystems)}, covering {total_reactions} reactions")

    selected = pd.read_csv(SELECTED_REACTIONS_PATH, index_col=0)
    selected_set = set(selected.index)
    touching = set(selected.index[selected["touches_variable_metabolite"].fillna("") != ""])
    proteomics_supported = set(selected.index[selected["proteomics_score"] > 0])

    rows = []
    for name, reactions in subsystems.items():
        if len(reactions) < MIN_SUBSYSTEM_SIZE:
            continue
        reaction_set = set(reactions)
        n_total = len(reaction_set)
        n_selected = len(reaction_set & selected_set)
        n_touching = len(reaction_set & touching)
        n_proteomics = len(reaction_set & proteomics_supported)
        rows.append({
            "subsystem": name,
            "n_total": n_total,
            "n_selected": n_selected,
            "pct_selected": 100 * n_selected / n_total,
            "pct_touching_of_selected": 100 * n_touching / n_selected if n_selected else 0.0,
            "pct_proteomics_of_selected": 100 * n_proteomics / n_selected if n_selected else 0.0,
        })
    coverage = pd.DataFrame(rows).sort_values("pct_selected", ascending=False).reset_index(drop=True)
    coverage.to_csv(OUTDIR / "pathway_coverage.csv", index=False)
    print(f"Subsystems with >= {MIN_SUBSYSTEM_SIZE} reactions: {len(coverage)}")

    columns = [
        ("pct_selected", "% of pathway\nselected"),
        ("pct_touching_of_selected", "% selected touching\na variable metabolite"),
        ("pct_proteomics_of_selected", "% selected\nproteomics-supported"),
    ]
    matrix = coverage[[c for c, _ in columns]].to_numpy()

    n_rows = len(coverage)
    fig, ax = plt.subplots(figsize=(7.5, max(6, 0.16 * n_rows)))
    im = ax.imshow(matrix, aspect="auto", cmap=BLUE_RAMP, vmin=0, vmax=100)

    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(
        [f"{name} ({n})" for name, n in zip(coverage["subsystem"], coverage["n_total"])], fontsize=6,
    )
    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels([label for _, label in columns], fontsize=8)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("%", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    ax.set_title(
        f"Script 14 network: pathway-level coverage ({len(selected)}/{12987} reactions selected)\n"
        f"subsystems sorted by % selected, label shows (reactions in subsystem)",
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(OUTDIR / "pathway_coverage_heatmap.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUTDIR / 'pathway_coverage_heatmap.png'}")


if __name__ == "__main__":
    main()
