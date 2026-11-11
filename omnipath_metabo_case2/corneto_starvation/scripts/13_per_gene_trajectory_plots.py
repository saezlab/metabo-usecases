"""Per-colleague follow-up on Script 12's per-gene R^2 table: pick one gene
with large R^2 in WT only, one in ob/ob only, and one in both, and plot the
protein/transcript trajectories across the 8 timepoints so we can eyeball
whether the correlation (and its genotype-specificity) makes biological
sense.

Gene choices (high R^2 within each classification, filtered for a
recognizable metabolic role rather than a structural/unannotated gene;
"both" picked by highest min(R2_WT, R2_ob/ob), not max, so it's genuinely
strong in both genotypes rather than just clearing the 0.5 threshold in
one of them):
  both:        Scd1 -- stearoyl-CoA desaturase 1, canonical lipogenic
               enzyme, fasting-suppressed in both genotypes
  WT only:     Acadm -- medium-chain acyl-CoA dehydrogenase (fatty-acid
               beta-oxidation), induced by fasting
  ob/ob only:  Fasn -- fatty acid synthase, canonical lipogenic enzyme

Values are z-scored per gene per modality across all 16 states (same
convention as Script 10), purely so protein and transcript -- different
units -- sit on one axis; it does not change the correlation itself.

Output (corneto_starvation/result/QC/):
  per_gene_trajectory_examples.png
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
R2_PATH = REPO_ROOT / "result" / "QC" / "proteomics_vs_transcriptomics_per_gene_r2.csv"
OUTDIR = REPO_ROOT / "result" / "QC"

TIMEPOINT_ORDER = (0, 2, 4, 6, 8, 12, 16, 24)
GENOTYPES = ("WT", "ob_ob")
GENES = ("Scd1", "Acadm", "Fasn")  # (both, WT only, ob/ob only)


def zscore_row(row: pd.Series) -> pd.Series:
    return (row - row.mean()) / row.std()


def main() -> None:
    protein = pd.read_csv(PROTEIN_PATH, index_col=0)
    transcript = pd.read_csv(TRANSCRIPT_PATH, index_col=0)
    r2 = pd.read_csv(R2_PATH, index_col=0)

    all_cols = [f"{g}_{t}h" for g in GENOTYPES for t in TIMEPOINT_ORDER]
    prot_z = protein.loc[list(GENES), all_cols].apply(zscore_row, axis=1)
    tx_z = transcript.loc[list(GENES), all_cols].apply(zscore_row, axis=1)

    fig, axes = plt.subplots(len(GENES), len(GENOTYPES), figsize=(9, 3 * len(GENES)), sharex=True)

    for i, gene in enumerate(GENES):
        row = r2.loc[gene]
        for j, geno in enumerate(GENOTYPES):
            ax = axes[i, j]
            cols = [f"{geno}_{t}h" for t in TIMEPOINT_ORDER]
            ax.plot(TIMEPOINT_ORDER, prot_z.loc[gene, cols], marker="o", color="#0d6e64", label="protein")
            ax.plot(TIMEPOINT_ORDER, tx_z.loc[gene, cols], marker="s", linestyle="--", color="#a73935", label="transcript")
            r2_val = row["r2_WT"] if geno == "WT" else row["r2_ob_ob"]
            ax.set_title(f"{gene} -- {geno}  (R²={r2_val:.2f})", fontsize=10)
            ax.axhline(0, color="#cccccc", linewidth=0.7, zorder=0)
            if i == len(GENES) - 1:
                ax.set_xlabel("hours fasted")
            if j == 0:
                ax.set_ylabel("z-score")
            if i == 0 and j == 0:
                ax.legend(fontsize=8, loc="best")

    fig.suptitle(f"Trajectories for one gene per classification\n({', '.join(GENES)} = both / WT only / ob-ob only)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUTDIR / "per_gene_trajectory_examples.png", dpi=180)
    plt.close(fig)
    print(f"wrote {OUTDIR / 'per_gene_trajectory_examples.png'}")

    print(r2.loc[list(GENES), ["r2_WT", "r2_ob_ob", "classification"]])


if __name__ == "__main__":
    main()
