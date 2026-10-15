"""Build presentation-ready heatmaps summarizing the input data: how each
modality was selected/QC'd, and what it looks like across the 16
(genotype, timepoint) states. One PNG per modality, row-standardized
(z-score per gene/metabolite) so each molecule's own relative pattern over
time is visible regardless of its absolute abundance scale, with rows
hierarchically clustered for readability and a visible WT | ob_ob divider.

Outputs (corneto_starvation/result/QC/):
  proteomics_top_variable.png   top-120-variance genes (3452 total -- too
                                 many to show individually, so we show the
                                 genes whose relative abundance moves the
                                 most across the time course)
  metabolite_tissue.png         all 171 tissue metabolites
  metabolite_plasma.png         all 98 plasma metabolites
  lipid.png                     all 16 lipid species
  data_overview_summary.csv     one row per modality: n measured, QC rule,
                                 replicate design, coverage
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.cluster.hierarchy import leaves_list, linkage

REPO_ROOT = Path(__file__).resolve().parents[1]
OUTDIR = REPO_ROOT / "result" / "QC"

GENE_SCORES = REPO_ROOT / "result" / "gene_scores" / "all_states_gene_scores.csv"
MET_DIR = REPO_ROOT / "result" / "metabolite_scores"

TIMEPOINT_ORDER = (0, 2, 4, 6, 8, 12, 16, 24)
STATE_ORDER = [f"{g}_{t}h" for g in ("WT", "ob_ob") for t in TIMEPOINT_ORDER]


def cluster_order(z: np.ndarray) -> np.ndarray:
    if z.shape[0] < 3:
        return np.arange(z.shape[0])
    finite = np.nan_to_num(z, nan=0.0)
    link = linkage(finite, method="average", metric="euclidean")
    return leaves_list(link)


def zscore_rows(wide: pd.DataFrame) -> pd.DataFrame:
    mean = wide.mean(axis=1, skipna=True)
    std = wide.std(axis=1, skipna=True).replace(0, np.nan)
    return wide.sub(mean, axis=0).div(std, axis=0)


def plot_heatmap(
    wide: pd.DataFrame,
    title: str,
    outpath: Path,
    *,
    row_label_every: int = 1,
    figsize_height_per_row: float = 0.16,
) -> None:
    wide = wide[STATE_ORDER]
    z = zscore_rows(wide)
    order = cluster_order(z.to_numpy())
    z = z.iloc[order]

    n_rows = z.shape[0]
    fig_h = max(3.0, n_rows * figsize_height_per_row)
    fig, ax = plt.subplots(figsize=(10, fig_h))

    vmax = np.nanpercentile(np.abs(z.to_numpy()), 95)
    vmax = max(vmax, 0.5)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    im = ax.imshow(z.to_numpy(), aspect="auto", cmap="RdBu_r", norm=norm, interpolation="nearest")

    ax.set_xticks(range(len(STATE_ORDER)))
    ax.set_xticklabels(STATE_ORDER, rotation=90, fontsize=7)
    ax.axvline(7.5, color="black", linewidth=1.5)

    if row_label_every and n_rows <= 180:
        ax.set_yticks(range(0, n_rows, row_label_every))
        ax.set_yticklabels(z.index[::row_label_every], fontsize=5.5)
    else:
        ax.set_yticks([])

    ax.set_title(title, fontsize=12, fontweight="bold")
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("z-score (per row, relative to its own 16-state mean)", fontsize=8)
    fig.tight_layout()
    fig.savefig(outpath, dpi=180)
    plt.close(fig)
    print(f"wrote {outpath} ({n_rows} rows)")


def build_proteomics_heatmap(top_n: int = 120) -> None:
    wide = pd.read_csv(GENE_SCORES, index_col=0)
    variance = wide.var(axis=1, skipna=True)
    top_genes = variance.sort_values(ascending=False).head(top_n).index
    plot_heatmap(
        wide.loc[top_genes],
        f"Proteomics: top {top_n} most-variable genes (of {wide.shape[0]} measured)",
        OUTDIR / "proteomics_top_variable.png",
    )


def build_metabolomics_heatmap(modality: str, label: str) -> pd.DataFrame:
    wide = pd.read_csv(MET_DIR / modality / "all_states_scores.csv", index_col=0)
    plot_heatmap(wide, f"{label} ({wide.shape[0]} measured, full set)", OUTDIR / f"{modality}.png")
    return wide


def write_summary() -> None:
    rows = []

    prot = pd.read_csv(GENE_SCORES, index_col=0)
    rows.append({
        "modality": "Proteomics",
        "n_measured": prot.shape[0],
        "source": "data/ads2547_data_file_s1.h5 /protein",
        "replicate_design": "16 states x 5 replicates",
        "qc_rule": "log2(x+1); keep if >=3/5 replicates non-missing; median-aggregate",
    })
    for modality, label, key in [
        ("metabolite", "Tissue metabolites (intracellular)", "/metabolite"),
        ("plasma_metabolite", "Plasma metabolites (extracellular)", "/plasma_metabolite"),
        ("lipid", "Lipid species", "/lipid"),
    ]:
        wide = pd.read_csv(MET_DIR / modality / "all_states_scores.csv", index_col=0)
        rows.append({
            "modality": label,
            "n_measured": wide.shape[0],
            "source": f"data/ads2547_data_file_s1.h5 {key}",
            "replicate_design": "16 states x 5 replicates",
            "qc_rule": "log2(x+1); keep if >=3/5 replicates non-missing; median-aggregate",
        })
    pd.DataFrame(rows).to_csv(OUTDIR / "data_overview_summary.csv", index=False)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    build_proteomics_heatmap()
    build_metabolomics_heatmap("metabolite", "Tissue metabolites (intracellular)")
    build_metabolomics_heatmap("plasma_metabolite", "Plasma metabolites (extracellular)")
    build_metabolomics_heatmap("lipid", "Lipid species")
    write_summary()


if __name__ == "__main__":
    main()
