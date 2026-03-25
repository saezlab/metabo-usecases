"""
Convert selected sheets from ads2547_data_file_s1.xlsx and ads2547_data_file_s2.xlsx
to separate HDF5 files.

Outputs:
  data/ads2547_data_file_s1.h5  — raw time-course measurements (S1)
  data/ads2547_data_file_s2.h5  — DEA results: AUC, t-half, p/q-values (S2)

Each sheet is stored as a named key: /Sample_information, /metabolite, etc.
"""

import tables  # required by pandas HDFStore
import pandas as pd
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "data"

S1_EXCEL = DATA_DIR / "ads2547_data_file_s1.xlsx"
S1_H5    = DATA_DIR / "ads2547_data_file_s1.h5"
S1_SHEETS = [
    "Sample_information",
    "metabolite",
    "lipid",
    "FFAandAcyls",
    "plasma metabolite",
    "protein",
    "Lipid_all",
    "Acylcarnitine_AcylCoA_all",
]

S2_EXCEL = DATA_DIR / "ads2547_data_file_s2.xlsx"
S2_H5    = DATA_DIR / "ads2547_data_file_s2.h5"
S2_SHEETS = [
    "metabolite",
    "lipid",
    "FFAandAcyls",
    "protein",
    "plasma metabolite",
]


def sheet_to_key(sheet: str) -> str:
    return sheet.replace(" ", "_")


def convert(excel_path: Path, h5_path: Path, sheets: list[str]) -> None:
    print(f"\n{excel_path.name}  →  {h5_path.name}")
    with pd.HDFStore(h5_path, mode="w", complevel=5, complib="blosc") as store:
        for sheet in sheets:
            print(f"  Reading '{sheet}' ...", end=" ", flush=True)
            df = pd.read_excel(excel_path, sheet_name=sheet)
            key = sheet_to_key(sheet)
            store.put(key, df, format="fixed")
            print(f"saved → /{key}  ({len(df):,} rows × {df.shape[1]} cols)")


def main():
    convert(S1_EXCEL, S1_H5, S1_SHEETS)
    convert(S2_EXCEL, S2_H5, S2_SHEETS)


if __name__ == "__main__":
    main()
