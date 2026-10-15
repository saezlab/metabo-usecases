"""Build per-(genotype,timepoint)-state metabolomics abundance scores from
the Morita et al. liver metabolomics time course (data/ads2547_data_file_s1.h5),
mirroring Script 01's protein pipeline exactly so the two modalities are
directly comparable.

Covers three modalities stored in s1.h5 (raw per-replicate values, same
16-state x 5-replicate design as /protein):
  /metabolite          171 tissue (intracellular) metabolites
  /plasma_metabolite    98 plasma (extracellular) metabolites
  /lipid                16 lipid species

16 states = 2 genotypes (WT, ob/ob) x 8 timepoints (0/2/4/6/8/12/16/24h),
5 replicates each. For each state: log2(x+1) transform, drop
metabolites with fewer than MIN_REPLICATES non-missing replicates,
median-aggregate the rest -- identical QC rule to Script 01.

Also computes, per modality, adjacent-timepoint fold change (Tn vs Tn-1,
within each genotype's own trajectory) on the *linear* (pre-log) median
intensity, for the qualitative up/down direction Pablo's metabolomics-anchored
network method needs. T0 has no Tn-1 and gets FC=NaN / direction="NA".

Outputs (corneto_starvation/result/metabolite_scores/<modality>/):
  <state>_scores.csv         one per state, columns metabolite_name,score
  all_states_scores.csv      wide table, metabolites x 16 states (log2 scores)
  qc_summary.csv             per-state kept/dropped counts and quantiles
  adjacent_fold_change.csv   metabolites x 16 states, FC vs previous timepoint
  adjacent_direction.csv     metabolites x 16 states, "Up"/"Down"/"NS"-like call
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_H5 = REPO_ROOT / "data" / "ads2547_data_file_s1.h5"
OUTDIR = Path(__file__).resolve().parents[1] / "result" / "metabolite_scores"

MIN_REPLICATES = 3
# A metabolite direction is called "Up"/"Down" only if |log2 FC| exceeds this;
# otherwise "NS" (not meaningfully changed from the previous timepoint).
FC_LOG2_THRESHOLD = 0.3  # ~1.23-fold

MODALITIES = {
    "metabolite": "/metabolite",
    "plasma_metabolite": "/plasma_metabolite",
    "lipid": "/lipid",
}

TIMEPOINT_ORDER = (0, 2, 4, 6, 8, 12, 16, 24)


def parse_row_name(row_id: str) -> str:
    """"metabolite;F6P;C00085" / "plasma metabolite;Lactate;C00186" -> "F6P" / "Lactate"."""
    fields = row_id.split(";")
    if len(fields) != 3 or not fields[1]:
        raise ValueError(f"Unexpected metabolomics Row format: {row_id!r}")
    return fields[1]


def sanitize_genotype(genotype: str) -> str:
    return {"WT": "WT", "ob/ob": "ob_ob"}[genotype]


def load_long_table(hdf_key: str) -> pd.DataFrame:
    with pd.HDFStore(DATA_H5, mode="r") as store:
        raw = store[hdf_key]
        sample_info = store["/Sample_information"]

    names = [parse_row_name(row) for row in raw["Row"]]
    if len(set(names)) != len(names):
        dup = pd.Series(names)[pd.Series(names).duplicated()].unique()
        raise ValueError(f"Duplicate metabolite names in {hdf_key} Row column: {dup[:10]!r}")

    wide = raw.drop(columns="Row").copy()
    wide.index = names
    wide.index.name = "metabolite_name"

    long = wide.reset_index().melt(
        id_vars="metabolite_name", var_name="MouseID", value_name="intensity"
    )
    long = long.merge(sample_info, on="MouseID", how="left")
    if long["genotype"].isna().any():
        missing = long.loc[long["genotype"].isna(), "MouseID"].unique()
        raise ValueError(f"Sample IDs in {hdf_key} not found in /Sample_information: {missing!r}")
    return long


def aggregate_state(
    long_df: pd.DataFrame,
    genotype: str,
    time: int,
    *,
    min_replicates: int = MIN_REPLICATES,
) -> tuple[pd.Series, pd.Series, dict]:
    """Returns (log2(x+1) median score, linear median intensity, qc dict)."""
    state_df = long_df[(long_df["genotype"] == genotype) & (long_df["time"] == time)]
    n_replicates_total = state_df["MouseID"].nunique()

    work = state_df.copy()
    work["log_intensity"] = np.log2(work["intensity"] + 1.0)

    counts = work.groupby("metabolite_name")["log_intensity"].count()
    kept = counts[counts >= min_replicates].index

    kept_df = work[work["metabolite_name"].isin(kept)]
    log_scores = kept_df.groupby("metabolite_name")["log_intensity"].median()
    linear_scores = kept_df.groupby("metabolite_name")["intensity"].median()

    qc = {
        "n_replicates_total": n_replicates_total,
        "n_metabolites_total": long_df["metabolite_name"].nunique(),
        "n_metabolites_kept": len(log_scores),
        "n_metabolites_dropped_low_replicate": long_df["metabolite_name"].nunique() - len(log_scores),
        "score_min": float(log_scores.min()) if len(log_scores) else float("nan"),
        "score_median": float(log_scores.median()) if len(log_scores) else float("nan"),
        "score_max": float(log_scores.max()) if len(log_scores) else float("nan"),
    }
    return log_scores, linear_scores, qc


def build_all_states(long_df: pd.DataFrame) -> tuple[dict[str, pd.Series], dict[str, pd.Series], pd.DataFrame]:
    log_states: dict[str, pd.Series] = {}
    linear_states: dict[str, pd.Series] = {}
    qc_rows = []
    combos = (
        long_df[["genotype", "time"]].drop_duplicates().sort_values(["genotype", "time"]).itertuples(index=False)
    )
    for genotype, time in combos:
        state_key = f"{sanitize_genotype(genotype)}_{time}h"
        log_scores, linear_scores, qc = aggregate_state(long_df, genotype, time)
        log_states[state_key] = log_scores
        linear_states[state_key] = linear_scores
        qc_rows.append({"state": state_key, "genotype": genotype, "time": time, **qc})
    qc_df = pd.DataFrame(qc_rows)
    return log_states, linear_states, qc_df


def compute_adjacent_fold_change(linear_wide: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """log2 fold change of each state's linear median intensity vs the SAME
    genotype's previous timepoint (T0 has no predecessor -> NaN/"NA")."""
    fc = pd.DataFrame(index=linear_wide.index, columns=linear_wide.columns, dtype=float)
    direction = pd.DataFrame(index=linear_wide.index, columns=linear_wide.columns, dtype=object)

    for genotype in ("WT", "ob_ob"):
        ordered_states = [f"{genotype}_{t}h" for t in TIMEPOINT_ORDER if f"{genotype}_{t}h" in linear_wide.columns]
        prev_state = None
        for state in ordered_states:
            if prev_state is None:
                fc[state] = np.nan
                direction[state] = "NA"
            else:
                ratio = (linear_wide[state] + 1.0) / (linear_wide[prev_state] + 1.0)
                log2fc = np.log2(ratio)
                fc[state] = log2fc
                direction[state] = np.select(
                    [log2fc >= FC_LOG2_THRESHOLD, log2fc <= -FC_LOG2_THRESHOLD],
                    ["Up", "Down"],
                    default="NS",
                )
                direction.loc[log2fc.isna(), state] = "NA"
            prev_state = state
    return fc, direction


def write_outputs(
    modality: str,
    log_states: dict[str, pd.Series],
    linear_states: dict[str, pd.Series],
    qc_df: pd.DataFrame,
    outdir: Path,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    for state_key, scores in log_states.items():
        scores.rename("score").rename_axis("metabolite_name").reset_index().to_csv(
            outdir / f"{state_key}_scores.csv", index=False
        )

    log_wide = pd.DataFrame(log_states)
    log_wide.index.name = "metabolite_name"
    log_wide = log_wide[sorted(log_wide.columns)]
    log_wide.to_csv(outdir / "all_states_scores.csv")

    linear_wide = pd.DataFrame(linear_states)
    linear_wide.index.name = "metabolite_name"
    linear_wide = linear_wide[sorted(linear_wide.columns)]

    fc, direction = compute_adjacent_fold_change(linear_wide)
    fc.rename_axis("metabolite_name").to_csv(outdir / "adjacent_fold_change.csv")
    direction.rename_axis("metabolite_name").to_csv(outdir / "adjacent_direction.csv")

    qc_df.to_csv(outdir / "qc_summary.csv", index=False)
    print(f"[{modality}] wrote {len(log_states)} state files to {outdir}")


def main() -> None:
    for modality, hdf_key in MODALITIES.items():
        long_df = load_long_table(hdf_key)
        log_states, linear_states, qc_df = build_all_states(long_df)
        write_outputs(modality, log_states, linear_states, qc_df, OUTDIR / modality)
        print(qc_df[["state", "n_metabolites_kept", "n_metabolites_dropped_low_replicate"]].to_string(index=False))
        print()


if __name__ == "__main__":
    main()
