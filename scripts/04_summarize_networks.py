"""Summarize the 16 reconstructed state networks as pathway-level figures, so
the rewiring can be read without opening every per-state network.

From 03's reaction x state membership matrix, count active reactions per
pathway at each state. Produces, per genotype, a pathway x timepoint activity
heatmap (how many reactions each pathway uses over the time course) and a
WT - ob/ob difference heatmap (which pathways each genotype runs more).

Outputs (result/state_network_reconstruction/):
  pathway_activity.csv              pathway x state active-reaction counts
  pathway_activity_heatmap.png      WT and ob/ob activity, pathway x timepoint
  pathway_difference_heatmap.png    WT - ob/ob, pathway x timepoint
  pathway_fraction_heatmap.png      activity as a fraction of each pathway's
                                    reactions (0-100%), so pathways of different
                                    sizes are comparable
  network_size_trajectory.png       total reactions per state over time, WT vs ob/ob
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

PKG_SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(PKG_SRC))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

network = importlib.import_module("01_metabolic_task_network")  # noqa: E402

PKG_ROOT = Path(__file__).resolve().parents[1]
OUTDIR = PKG_ROOT / "result" / "state_network_reconstruction"
MEMBERSHIP = OUTDIR / "selected_reactions_by_state.csv"

GENOTYPES = ("WT", "ob_ob")
STATE_SUFFIXES = ("0h", "2h", "4h", "6h", "8h", "12h", "16h", "24h")


def _annotate(ax, values, fmt="{:.0f}", dark_above=None):
    """Write each cell's value; flip text color on dark cells for contrast."""
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            v = values[i, j]
            color = "white" if (dark_above is not None and abs(v) >= dark_above) else "#1b1f23"
            ax.text(j, i, fmt.format(v), ha="center", va="center", fontsize=7, color=color)


def main() -> None:
    if not MEMBERSHIP.exists():
        raise RuntimeError(
            f"{MEMBERSHIP} not found -- run 03_state_network_reconstruction.py first."
        )

    graph = network.build_full_metabolic_network()
    pathway = {
        graph.get_attr_edge(i)["id"]: (graph.get_attr_edge(i).get("pathway") or "Unassigned")
        for i in range(graph.num_edges)
    }

    membership = pd.read_csv(MEMBERSHIP, index_col=0).astype(bool)
    membership["pathway"] = membership.index.map(pathway)
    states = [f"{g}_{t}" for g in GENOTYPES for t in STATE_SUFFIXES]
    activity = membership.groupby("pathway")[states].sum()

    # Order pathways by how WT-biased they are, so the difference reads as a
    # gradient (WT-heavy pathways on top, ob/ob-heavy at the bottom).
    wt_cols = [f"WT_{t}" for t in STATE_SUFFIXES]
    ob_cols = [f"ob_ob_{t}" for t in STATE_SUFFIXES]
    bias = activity[wt_cols].mean(axis=1) - activity[ob_cols].mean(axis=1)
    activity = activity.loc[bias.sort_values(ascending=False).index]
    activity.to_csv(OUTDIR / "pathway_activity.csv")

    pathways = list(activity.index)
    hours = list(STATE_SUFFIXES)
    wt = activity[wt_cols].to_numpy()
    ob = activity[ob_cols].to_numpy()

    # Fraction of each pathway's reactions that are active (fixed denominator =
    # the pathway's size in the 117-reaction network), so pathways of different
    # sizes and states of different total size are comparable.
    pathway_total = membership.groupby("pathway").size().reindex(pathways).to_numpy()
    wt_frac = 100.0 * wt / pathway_total[:, None]
    ob_frac = 100.0 * ob / pathway_total[:, None]

    # --- Figure 1: per-genotype activity (sequential: magnitude, one hue) ---
    vmax = int(max(wt.max(), ob.max()))
    fig, axes = plt.subplots(1, 2, figsize=(12, 7), sharey=True)
    for ax, data, name in zip(axes, (wt, ob), ("WT", "ob/ob")):
        im = ax.imshow(data, aspect="auto", cmap="Blues", vmin=0, vmax=vmax)
        ax.set_xticks(range(len(hours)), hours, fontsize=8)
        ax.set_yticks(range(len(pathways)), pathways, fontsize=8)
        ax.set_title(name, fontsize=11, fontweight="bold")
        ax.set_xlabel("timepoint", fontsize=9)
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        _annotate(ax, data, dark_above=vmax * 0.6)
    fig.colorbar(im, ax=axes, label="active reactions", fraction=0.025, pad=0.02)
    fig.suptitle("Pathway activity in the reconstructed state networks", fontsize=12)
    fig.savefig(OUTDIR / "pathway_activity_heatmap.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    # --- Figure 2: WT - ob/ob (diverging: polarity, two hues + gray midpoint) ---
    diff = wt - ob
    bound = int(max(1, np.abs(diff).max()))
    fig, ax = plt.subplots(figsize=(7.5, 7))
    im = ax.imshow(diff, aspect="auto", cmap="RdBu_r", vmin=-bound, vmax=bound)
    ax.set_xticks(range(len(hours)), hours, fontsize=8)
    ax.set_yticks(range(len(pathways)), pathways, fontsize=8)
    ax.set_xlabel("timepoint", fontsize=9)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    _annotate(ax, diff, fmt="{:+.0f}", dark_above=bound * 0.6)
    cbar = fig.colorbar(im, ax=ax, label="WT - ob/ob active reactions", fraction=0.045, pad=0.03)
    cbar.ax.text(0.5, 1.02, "WT more", ha="center", va="bottom", transform=cbar.ax.transAxes, fontsize=7)
    cbar.ax.text(0.5, -0.02, "ob/ob more", ha="center", va="top", transform=cbar.ax.transAxes, fontsize=7)
    ax.set_title("Genotype difference in pathway activity", fontsize=12)
    fig.savefig(OUTDIR / "pathway_difference_heatmap.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    # --- Figure 3: activity as a fraction of pathway size (0-100%) ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 7), sharey=True)
    for ax, data, name in zip(axes, (wt_frac, ob_frac), ("WT", "ob/ob")):
        im = ax.imshow(data, aspect="auto", cmap="Blues", vmin=0, vmax=100)
        ax.set_xticks(range(len(hours)), hours, fontsize=8)
        ax.set_yticks(range(len(pathways)), pathways, fontsize=8)
        ax.set_title(name, fontsize=11, fontweight="bold")
        ax.set_xlabel("timepoint", fontsize=9)
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        _annotate(ax, data, fmt="{:.0f}%", dark_above=60)
    fig.colorbar(im, ax=axes, label="% of pathway's reactions active", fraction=0.025, pad=0.02)
    fig.suptitle("Pathway activity as a fraction of pathway size", fontsize=12)
    fig.savefig(OUTDIR / "pathway_fraction_heatmap.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    # --- Figure 4: total network size over the time course (line) ---
    totals = membership[states].sum()
    wt_total = [int(totals[f"WT_{t}"]) for t in STATE_SUFFIXES]
    ob_total = [int(totals[f"ob_ob_{t}"]) for t in STATE_SUFFIXES]
    x = range(len(STATE_SUFFIXES))
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot(x, wt_total, marker="o", markersize=7, linewidth=2, color="#2a6f97")
    ax.plot(x, ob_total, marker="o", markersize=7, linewidth=2, color="#e8833a")
    ax.text(x[-1] + 0.12, wt_total[-1], "WT", color="#2a6f97", va="center", fontweight="bold", fontsize=10)
    ax.text(x[-1] + 0.12, ob_total[-1], "ob/ob", color="#e8833a", va="center", fontweight="bold", fontsize=10)
    ax.set_xticks(list(x), hours)
    ax.set_xlabel("timepoint", fontsize=9)
    ax.set_ylabel("reactions in reconstructed network", fontsize=9)
    ax.set_title("Reconstructed network size over the starvation time course", fontsize=11)
    ax.set_ylim(0, graph.num_edges)
    ax.grid(axis="y", color="#edf0f2")
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.savefig(OUTDIR / "network_size_trajectory.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote pathway_activity.csv and 4 figures to {OUTDIR}")
    print("\nMost WT-biased pathways (mean WT - ob/ob over the time course):")
    print(bias.sort_values(ascending=False).head(6).to_string())


if __name__ == "__main__":
    main()
