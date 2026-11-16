"""Per-gene, per-condition proteomics-vs-transcriptomics correlation, per
Pablo's suggestion: Script 10 pooled all (gene, state) pairs into one
overall R^2, which hides how individual genes behave. Here, for each
shared gene, compute Pearson r/R^2 twice -- once across the 8 WT
timepoints, once across the 8 ob/ob timepoints -- so we can see whether a
gene's protein and transcript levels track together in WT, in ob/ob,
both, or neither.

Pearson correlation is invariant to per-gene z-scoring (an affine
transform), so raw values are used directly -- no z-scoring needed here,
unlike Script 10 where z-scoring mattered because *different genes'*
values were pooled together.

Outputs (corneto_starvation/result/QC/):
  proteomics_vs_transcriptomics_per_gene_r2.csv   r/R^2/p per gene per condition, classification
  proteomics_vs_transcriptomics_per_gene_scatter.png   R^2_WT vs R^2_ob/ob, one point per gene
  proteomics_vs_transcriptomics_per_gene_top20.png   PNG table, genes ranked by max(R^2_WT, R^2_ob/ob)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

REPO_ROOT = Path(__file__).resolve().parents[1]
PROTEIN_PATH = REPO_ROOT / "result" / "gene_scores" / "all_states_gene_scores.csv"
TRANSCRIPT_PATH = REPO_ROOT / "result" / "transcript_scores" / "all_states_scores.csv"
OUTDIR = REPO_ROOT / "result" / "QC"

TIMEPOINT_ORDER = (0, 2, 4, 6, 8, 12, 16, 24)
WT_COLS = [f"WT_{t}h" for t in TIMEPOINT_ORDER]
OBOB_COLS = [f"ob_ob_{t}h" for t in TIMEPOINT_ORDER]

# Single threshold: a condition "counts" as correlated if its R^2 clears this.
R2_THRESHOLD = 0.5


def per_condition_corr(protein: pd.DataFrame, transcript: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    rows = []
    for gene in protein.index:
        x = transcript.loc[gene, cols].to_numpy(dtype=float)
        y = protein.loc[gene, cols].to_numpy(dtype=float)
        if np.std(x) == 0 or np.std(y) == 0:
            rows.append({"gene": gene, "r": np.nan, "r2": np.nan, "p": np.nan})
            continue
        r, p = pearsonr(x, y)
        rows.append({"gene": gene, "r": r, "r2": r ** 2, "p": p})
    return pd.DataFrame(rows).set_index("gene")


def classify(row: pd.Series) -> str:
    r2_wt, r2_ob = row["r2_WT"], row["r2_ob_ob"]
    if np.isnan(r2_wt) or np.isnan(r2_ob):
        return "undefined (zero variance in one condition)"
    wt_hit = r2_wt >= R2_THRESHOLD
    ob_hit = r2_ob >= R2_THRESHOLD
    if wt_hit and ob_hit:
        return "both"
    if wt_hit:
        return "WT only"
    if ob_hit:
        return "ob_ob only"
    return "neither"


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    protein = pd.read_csv(PROTEIN_PATH, index_col=0)
    transcript = pd.read_csv(TRANSCRIPT_PATH, index_col=0)

    shared = sorted(set(protein.index) & set(transcript.index))
    protein = protein.loc[shared]
    transcript = transcript.loc[shared]
    print(f"{len(shared)} shared genes")

    wt = per_condition_corr(protein, transcript, WT_COLS).add_suffix("_WT")
    ob = per_condition_corr(protein, transcript, OBOB_COLS).add_suffix("_ob_ob")
    merged = wt.join(ob)
    merged["classification"] = merged.apply(classify, axis=1)
    merged["r2_max"] = merged[["r2_WT", "r2_ob_ob"]].max(axis=1)
    merged = merged.sort_values("r2_max", ascending=False)

    n_undefined = (merged["classification"] == "undefined (zero variance in one condition)").sum()
    print(f"{n_undefined} genes dropped (zero variance within WT or ob/ob)")
    counts = merged["classification"].value_counts()
    print(counts.to_string())

    merged.to_csv(OUTDIR / "proteomics_vs_transcriptomics_per_gene_r2.csv")

    plot_df = merged.dropna(subset=["r2_WT", "r2_ob_ob"])

    # Scatter: R^2_ob/ob vs R^2_WT, one point per gene.
    color_map = {
        "both": "#0d6e64",
        "WT only": "#2b6cb0",
        "ob_ob only": "#a73935",
        "neither": "#999999",
    }
    fig, ax = plt.subplots(figsize=(7, 7))
    for label in ["both", "WT only", "ob_ob only", "neither"]:
        sub = plot_df[plot_df["classification"] == label]
        ax.scatter(
            sub["r2_ob_ob"], sub["r2_WT"], s=10, alpha=0.6,
            color=color_map[label], edgecolors="none", label=f"{label} (n={len(sub)})",
        )
    r_half_r2 = 0.5 ** 2  # R^2 corresponding to r = 0.5
    for value, style, color in [(R2_THRESHOLD, "--", "#999999"), (r_half_r2, ":", "#555555")]:
        ax.axhline(value, color=color, linewidth=1, zorder=0, linestyle=style)
        ax.axvline(value, color=color, linewidth=1, zorder=0, linestyle=style)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("R² (ob/ob, 8 timepoints)")
    ax.set_ylabel("R² (WT, 8 timepoints)")
    ax.set_title(f"Protein vs transcript agreement, per gene per condition\n{len(plot_df)} genes with defined R² in both conditions, classification threshold R²={R2_THRESHOLD}", fontsize=11)
    from matplotlib.lines import Line2D
    line_handles = [
        Line2D([0], [0], color="#999999", linestyle="--", label=f"R²={R2_THRESHOLD} (r≈{R2_THRESHOLD ** 0.5:.2f}) — classification threshold"),
        Line2D([0], [0], color="#555555", linestyle=":", label=f"R²={r_half_r2:.2f} (r=0.5) — reference only"),
    ]
    ax.legend(handles=ax.get_legend_handles_labels()[0] + line_handles, loc="upper left", fontsize=8, markerscale=1.5)
    fig.tight_layout()
    fig.savefig(OUTDIR / "proteomics_vs_transcriptomics_per_gene_scatter.png", dpi=180)
    plt.close(fig)
    print(f"wrote {OUTDIR / 'proteomics_vs_transcriptomics_per_gene_scatter.png'}")

    # Top-20 table (by max R^2 across the two conditions).
    top20 = merged.dropna(subset=["r2_WT", "r2_ob_ob"]).head(20).reset_index()
    table_data = top20[["gene", "r2_WT", "r2_ob_ob", "classification"]].copy()
    table_data["r2_WT"] = table_data["r2_WT"].round(3)
    table_data["r2_ob_ob"] = table_data["r2_ob_ob"].round(3)

    fig, ax = plt.subplots(figsize=(8, 0.33 * len(table_data) + 0.6))
    ax.axis("off")
    tbl = ax.table(
        cellText=table_data.values,
        colLabels=["gene", "R² (WT)", "R² (ob/ob)", "classification"],
        loc="center", cellLoc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 1.4)
    ax.set_title("Top 20 genes by max(R²_WT, R²_ob/ob)", pad=10)
    fig.savefig(OUTDIR / "proteomics_vs_transcriptomics_per_gene_top20.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUTDIR / 'proteomics_vs_transcriptomics_per_gene_top20.png'}")


if __name__ == "__main__":
    main()
