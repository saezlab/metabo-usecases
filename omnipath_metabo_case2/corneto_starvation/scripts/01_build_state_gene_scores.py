"""Build per-(genotype,timepoint)-state gene abundance scores from the
Morita et al. liver proteomics time course (data/ads2547_data_file_s1.h5).

16 states = 2 genotypes (WT, ob/ob) x 8 timepoints (0/2/4/6/8/12/16/24h),
5 replicates each. For each state: log2(x+1) transform, drop genes with
fewer than MIN_REPLICATES non-missing replicates, median-aggregate the rest.

Outputs (corneto_starvation/result/gene_scores/):
  <state>_gene_scores.csv   one per state, columns gene_symbol,score
  all_states_gene_scores.csv  wide table, genes x 16 states
  qc_summary.csv            per-state kept/dropped gene counts and quantiles
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_H5 = REPO_ROOT / "data" / "ads2547_data_file_s1.h5"
OUTDIR = Path(__file__).resolve().parents[1] / "result" / "gene_scores"

MIN_REPLICATES = 3


def parse_row_gene_symbol(row_id: str) -> str:
    """"protein;Pck1;ENSMUSG...;mmu:..." -> "Pck1"."""
    fields = row_id.split(";")
    if len(fields) != 4 or fields[0] != "protein" or not fields[1]:
        raise ValueError(f"Unexpected /protein Row format: {row_id!r}")
    return fields[1]


def sanitize_genotype(genotype: str) -> str:
    return {"WT": "WT", "ob/ob": "ob_ob"}[genotype]


def load_long_protein_table() -> pd.DataFrame:
    with pd.HDFStore(DATA_H5, mode="r") as store:
        protein = store["/protein"]
        sample_info = store["/Sample_information"]

    gene_symbols = [parse_row_gene_symbol(row) for row in protein["Row"]]
    if len(set(gene_symbols)) != len(gene_symbols):
        duplicates = pd.Series(gene_symbols)[pd.Series(gene_symbols).duplicated()].unique()
        raise ValueError(f"Duplicate gene symbols in /protein Row column: {duplicates[:10]!r}")

    wide = protein.drop(columns="Row").copy()
    wide.index = gene_symbols
    wide.index.name = "gene_symbol"

    long = wide.reset_index().melt(
        id_vars="gene_symbol", var_name="MouseID", value_name="intensity"
    )
    long = long.merge(sample_info, on="MouseID", how="left")
    if long["genotype"].isna().any():
        missing = long.loc[long["genotype"].isna(), "MouseID"].unique()
        raise ValueError(f"Sample IDs in /protein not found in /Sample_information: {missing!r}")
    return long


def aggregate_state(
    long_df: pd.DataFrame,
    genotype: str,
    time: int,
    *,
    min_replicates: int = MIN_REPLICATES,
) -> tuple[pd.Series, dict]:
    state_df = long_df[(long_df["genotype"] == genotype) & (long_df["time"] == time)]
    n_replicates_total = state_df["MouseID"].nunique()

    log_df = state_df.copy()
    log_df["log_intensity"] = np.log2(log_df["intensity"] + 1.0)

    counts = log_df.groupby("gene_symbol")["log_intensity"].count()
    kept_genes = counts[counts >= min_replicates].index

    kept = log_df[log_df["gene_symbol"].isin(kept_genes)]
    scores = kept.groupby("gene_symbol")["log_intensity"].median()

    qc = {
        "n_replicates_total": n_replicates_total,
        "n_genes_total": long_df["gene_symbol"].nunique(),
        "n_genes_kept": len(scores),
        "n_genes_dropped_low_replicate": long_df["gene_symbol"].nunique() - len(scores),
        "score_min": float(scores.min()) if len(scores) else float("nan"),
        "score_median": float(scores.median()) if len(scores) else float("nan"),
        "score_max": float(scores.max()) if len(scores) else float("nan"),
    }
    return scores, qc


def build_all_states(long_df: pd.DataFrame) -> tuple[dict[str, pd.Series], pd.DataFrame]:
    states: dict[str, pd.Series] = {}
    qc_rows = []
    combos = (
        long_df[["genotype", "time"]].drop_duplicates().sort_values(["genotype", "time"]).itertuples(index=False)
    )
    for genotype, time in combos:
        state_key = f"{sanitize_genotype(genotype)}_{time}h"
        scores, qc = aggregate_state(long_df, genotype, time)
        states[state_key] = scores
        qc_rows.append({"state": state_key, "genotype": genotype, "time": time, **qc})
    qc_df = pd.DataFrame(qc_rows)
    return states, qc_df


def write_outputs(states: dict[str, pd.Series], qc_df: pd.DataFrame, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    for state_key, scores in states.items():
        scores.rename("score").rename_axis("gene_symbol").reset_index().to_csv(
            outdir / f"{state_key}_gene_scores.csv", index=False
        )

    wide = pd.DataFrame(states)
    wide.index.name = "gene_symbol"
    wide = wide[sorted(wide.columns)]
    wide.to_csv(outdir / "all_states_gene_scores.csv")

    qc_df.to_csv(outdir / "qc_summary.csv", index=False)


def main() -> None:
    long_df = load_long_protein_table()
    states, qc_df = build_all_states(long_df)
    write_outputs(states, qc_df, OUTDIR)
    print(f"Wrote {len(states)} state gene-score files to {OUTDIR}")
    print(qc_df.to_string(index=False))


if __name__ == "__main__":
    main()
