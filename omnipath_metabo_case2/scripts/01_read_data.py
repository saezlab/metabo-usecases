"""
Read data from ads2547_data_file_s1.h5 (measurements) and
ads2547_data_file_s2.h5 (DEA results).
"""

import tables
import pandas as pd
from pathlib import Path

DATA_DIR = Path("omnipath_metabo_case2/data")
S1_H5 = DATA_DIR / "ads2547_data_file_s1.h5"
S2_H5 = DATA_DIR / "ads2547_data_file_s2.h5"

# ── S1: time-course measurements ───────────────────────────────────────────────
with pd.HDFStore(S1_H5, mode="r") as store:
    sample_info   = store["/Sample_information"]
    metabolite    = store["/metabolite"]
    lipid         = store["/lipid"]
    ffa_acyls     = store["/FFAandAcyls"]
    plasma_metab  = store["/plasma_metabolite"]
    protein       = store["/protein"]
    lipid_all     = store["/Lipid_all"]
    acylcarnitine = store["/Acylcarnitine_AcylCoA_all"]

# ── S2: DEA results (AUC, t-half, p/q-values, WT/ob/ob change) ────────────────
with pd.HDFStore(S2_H5, mode="r") as store:
    dea_metabolite   = store["/metabolite"]
    dea_lipid        = store["/lipid"]
    dea_ffa_acyls    = store["/FFAandAcyls"]
    dea_protein      = store["/protein"]
    dea_plasma_metab = store["/plasma_metabolite"]

# ── Quick overview ─────────────────────────────────────────────────────────────
s1_datasets = {
    "sample_info":    sample_info,
    "metabolite":     metabolite,
    "lipid":          lipid,
    "ffa_acyls":      ffa_acyls,
    "plasma_metab":   plasma_metab,
    "protein":        protein,
    "lipid_all":      lipid_all,
    "acylcarnitine":  acylcarnitine,
}

s2_dea = {
    "dea_metabolite":   dea_metabolite,
    "dea_lipid":        dea_lipid,
    "dea_ffa_acyls":    dea_ffa_acyls,
    "dea_protein":      dea_protein,
    "dea_plasma_metab": dea_plasma_metab,
}

print("── S1 measurements ──────────────────────────────────────────────────────")
for name, df in s1_datasets.items():
    print(f"  {name:20s}  {df.shape[0]:>5} rows × {df.shape[1]:>3} cols")

print("\n── S2 DEA results ───────────────────────────────────────────────────────")
for name, df in s2_dea.items():
    print(f"  {name:20s}  {df.shape[0]:>5} rows × {df.shape[1]:>3} cols")
