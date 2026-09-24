from pathlib import Path

import pandas as pd
from omnipath_metabo.datasets import cosmos

# Paths relative to the repository root (this file: analyses/<name>/).
REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_FILE = REPO_ROOT / "data/derived/cancer-cell-lines/ExtendedDataTable_DifferentialAnalysis_Shorthouse_LungMutations.xlsx"
OUT_DIR = REPO_ROOT / "data/derived/cancer-cell-lines-cosmos"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Load and filter DEMs ──────────────────────────────────────────────────────

sheets = {"EGFR": "EGFR_filt_limma", "KRAS": "KRAS_filt_limma"}
dems = {}

for label, sheet in sheets.items():
    df = pd.read_excel(DATA_FILE, sheet_name=sheet)
    filtered = df[(df["P.Value"] < 0.05) & (df["logFC"].abs() > 0.5)]
    dems[label] = filtered
    print(f"\n{label}: {len(filtered)} DEMs (p<0.05, |logFC|>0.5)")
    print(filtered[["metabolite_name", "logFC", "P.Value", "chebi"]].to_string(index=False))

# ── Build COSMOS PKN ──────────────────────────────────────────────────────────

print("\n\n=== Building allosteric network ===")
allosteric = cosmos.build_allosteric()
df_a = pd.DataFrame(allosteric.network)
print(df_a.groupby(["resource", "interaction_type"]).size())
df_a.to_csv(f"{OUT_DIR}/pkn_allosteric.csv", index=False)
print(f"Saved: {OUT_DIR}/pkn_allosteric.csv  ({len(df_a)} edges)")

print("\n=== Building enzyme-metabolite network ===")
enzyme_met = cosmos.build_enzyme_metabolite(gem={"include_orphans": False})
df_e = pd.DataFrame(enzyme_met.network)
print(df_e.groupby(["resource", "interaction_type"]).size())
df_e.to_csv(f"{OUT_DIR}/pkn_enzyme_metabolite.csv", index=False)
print(f"Saved: {OUT_DIR}/pkn_enzyme_metabolite.csv  ({len(df_e)} edges)")

