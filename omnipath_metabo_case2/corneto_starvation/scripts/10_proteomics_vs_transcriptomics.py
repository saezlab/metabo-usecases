"""Compare proteomics and transcriptomics for the genes measured in both,
per Christina's request: how well does enzyme change agree between the two
modalities, to justify using proteomics as the primary evidence (higher
weight) and transcriptomics as a fallback only for enzymes without a
protein measurement.

Method: restrict to genes present in both Script 01's protein scores and
Script 09's transcript scores. Each gene's 16-state trajectory is
z-scored within its own modality (so absolute units -- protein intensity
vs TPM -- don't matter, only each gene's *relative pattern* across states
is compared, which is what "enzyme change" means here). Every
(gene, state) pair becomes one point; R^2 is the squared Pearson
correlation across all of them, plus the same broken out per state (to
see whether agreement is roughly constant across the time course or
concentrated in some states).

Outputs (corneto_starvation/result/QC/):
  proteomics_vs_transcriptomics_scatter.png   all (gene, state) points, R^2 annotated
  proteomics_vs_transcriptomics_per_state_r2.csv
  proteomics_vs_transcriptomics_gene_coverage.csv   how many genes overlap vs protein-only vs transcript-only
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[1]
PROTEIN_PATH = REPO_ROOT / "result" / "gene_scores" / "all_states_gene_scores.csv"
TRANSCRIPT_PATH = REPO_ROOT / "result" / "transcript_scores" / "all_states_scores.csv"
OUTDIR = REPO_ROOT / "result" / "QC"

TIMEPOINT_ORDER = (0, 2, 4, 6, 8, 12, 16, 24)
STATE_ORDER = [f"{g}_{t}h" for g in ("WT", "ob_ob") for t in TIMEPOINT_ORDER]


def zscore_rows(wide: pd.DataFrame) -> pd.DataFrame:
    mean = wide.mean(axis=1, skipna=True)
    std = wide.std(axis=1, skipna=True).replace(0, np.nan)
    return wide.sub(mean, axis=0).div(std, axis=0)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    protein = pd.read_csv(PROTEIN_PATH, index_col=0)[STATE_ORDER]
    transcript = pd.read_csv(TRANSCRIPT_PATH, index_col=0)[STATE_ORDER]

    protein_genes = set(protein.index)
    transcript_genes = set(transcript.index)
    shared = sorted(protein_genes & transcript_genes)

    coverage = pd.DataFrame([
        {"set": "protein only (no transcript measurement)", "n_genes": len(protein_genes - transcript_genes)},
        {"set": "transcript only (no protein measurement)", "n_genes": len(transcript_genes - protein_genes)},
        {"set": "measured in both", "n_genes": len(shared)},
    ])
    coverage.to_csv(OUTDIR / "proteomics_vs_transcriptomics_gene_coverage.csv", index=False)
    print(coverage.to_string(index=False))

    prot_z = zscore_rows(protein.loc[shared])
    tx_z = zscore_rows(transcript.loc[shared])

    # Drop genes with zero variance in either modality (std=NaN after zscore_rows).
    valid_genes = prot_z.dropna(how="any").index.intersection(tx_z.dropna(how="any").index)
    prot_z = prot_z.loc[valid_genes]
    tx_z = tx_z.loc[valid_genes]
    print(f"\n{len(shared)} shared genes; {len(valid_genes)} with non-zero variance in both modalities, used below")

    # Per-state R^2
    per_state_rows = []
    for state in STATE_ORDER:
        x = tx_z[state].to_numpy()
        y = prot_z[state].to_numpy()
        r = np.corrcoef(x, y)[0, 1]
        per_state_rows.append({"state": state, "n_genes": len(x), "r": r, "r2": r ** 2})
    per_state_df = pd.DataFrame(per_state_rows)
    per_state_df.to_csv(OUTDIR / "proteomics_vs_transcriptomics_per_state_r2.csv", index=False)
    print()
    print(per_state_df.to_string(index=False))

    # Overall scatter: every (gene, state) pair pooled.
    x_all = tx_z.to_numpy().ravel()
    y_all = prot_z.to_numpy().ravel()
    mask = np.isfinite(x_all) & np.isfinite(y_all)
    x_all, y_all = x_all[mask], y_all[mask]
    r_overall = np.corrcoef(x_all, y_all)[0, 1]
    r2_overall = r_overall ** 2

    slope, intercept = np.polyfit(x_all, y_all, 1)

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(x_all, y_all, s=4, alpha=0.08, color="#0d6e64", edgecolors="none", rasterized=True)
    xs = np.linspace(x_all.min(), x_all.max(), 100)
    ax.plot(xs, slope * xs + intercept, color="#a73935", linewidth=1.5, label=f"fit: y = {slope:.2f}x + {intercept:.2f}")
    ax.axhline(0, color="#cccccc", linewidth=0.8, zorder=0)
    ax.axvline(0, color="#cccccc", linewidth=0.8, zorder=0)
    ax.set_xlabel("Transcript z-score (per gene, across 16 states)")
    ax.set_ylabel("Protein z-score (per gene, across 16 states)")
    ax.set_title(f"Proteomics vs transcriptomics\n{len(valid_genes)} shared genes x 16 states = {len(x_all)} points, R² = {r2_overall:.3f}")
    ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUTDIR / "proteomics_vs_transcriptomics_scatter.png", dpi=180)
    plt.close(fig)
    print(f"\nOverall R² = {r2_overall:.4f} (r = {r_overall:.4f}) across {len(x_all)} (gene, state) points")
    print(f"wrote {OUTDIR / 'proteomics_vs_transcriptomics_scatter.png'}")


if __name__ == "__main__":
    main()
