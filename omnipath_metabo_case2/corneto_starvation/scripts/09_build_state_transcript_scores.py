"""Build per-(genotype,timepoint)-state transcript (TPM) scores from the
Morita et al. liver transcriptomics time course. Mirrors Script 01's
protein pipeline exactly, with one addition: 8 of 14,025 transcript rows
share a gene symbol across two Ensembl IDs (old/new annotation of the same
gene) -- their TPM values are summed per sample before aggregation, since
TPM is additive across transcripts/annotations of the same gene.

Source: data/ads2547_data_file_s1.xlsx, sheet "transcript" -- this sheet
was NOT carried into ads2547_data_file_s1.h5 (only Sample_information,
metabolite, lipid, FFAandAcyls, plasma metabolite, protein were), so this
script reads the .xlsx directly instead of the .h5.

16 states = 2 genotypes (WT, ob/ob) x 8 timepoints (0/2/4/6/8/12/16/24h),
5 replicates each. For each state: log2(TPM+1) transform, drop genes with
fewer than MIN_REPLICATES non-missing replicates, median-aggregate the rest
-- identical QC rule to Script 01's proteomics pipeline, so the two
modalities are directly comparable gene-by-gene.

Outputs (corneto_starvation/result/transcript_scores/):
  <state>_scores.csv        one per state, columns gene_symbol,score
  all_states_scores.csv     wide table, genes x 16 states
  qc_summary.csv            per-state kept/dropped gene counts and quantiles
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_XLSX = REPO_ROOT / "data" / "ads2547_data_file_s1.xlsx"
OUTDIR = Path(__file__).resolve().parents[1] / "result" / "transcript_scores"

MIN_REPLICATES = 3


def parse_row_gene_symbol(row_id: str) -> str:
    """"transcript;Serpina1e;ENSMUSG00000072849;mmu:20704" -> "Serpina1e"."""
    fields = row_id.split(";")
    if len(fields) != 4 or fields[0] != "transcript" or not fields[1]:
        raise ValueError(f"Unexpected transcript Row format: {row_id!r}")
    return fields[1]


def sanitize_genotype(genotype: str) -> str:
    return {"WT": "WT", "ob/ob": "ob_ob"}[genotype]


def load_long_transcript_table() -> pd.DataFrame:
    with pd.HDFStore(REPO_ROOT / "data" / "ads2547_data_file_s1.h5", mode="r") as store:
        sample_info = store["/Sample_information"]

    raw = pd.read_excel(DATA_XLSX, sheet_name="transcript")
    gene_symbols = [parse_row_gene_symbol(row) for row in raw["Row"]]

    wide = raw.drop(columns="Row").copy()
    wide.index = gene_symbols
    wide.index.name = "gene_symbol"
    # Sum TPM across rows sharing a gene symbol (old/new Ensembl annotation
    # of the same gene) before anything else touches the data.
    wide = wide.groupby("gene_symbol").sum()

    long = wide.reset_index().melt(
        id_vars="gene_symbol", var_name="MouseID", value_name="intensity"
    )
    long = long.merge(sample_info, on="MouseID", how="left")
    if long["genotype"].isna().any():
        missing = long.loc[long["genotype"].isna(), "MouseID"].unique()
        raise ValueError(f"Sample IDs in transcript sheet not found in /Sample_information: {missing!r}")
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
            outdir / f"{state_key}_scores.csv", index=False
        )

    wide = pd.DataFrame(states)
    wide.index.name = "gene_symbol"
    wide = wide[sorted(wide.columns)]
    wide.to_csv(outdir / "all_states_scores.csv")

    qc_df.to_csv(outdir / "qc_summary.csv", index=False)


def main() -> None:
    long_df = load_long_transcript_table()
    states, qc_df = build_all_states(long_df)
    write_outputs(states, qc_df, OUTDIR)
    print(f"Wrote {len(states)} state transcript-score files to {OUTDIR}")
    print(qc_df.to_string(index=False))


if __name__ == "__main__":
    main()
