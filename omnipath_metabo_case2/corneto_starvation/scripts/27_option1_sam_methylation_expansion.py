"""Option 1, third expansion round: the SAM "methylation economy" module
group, added on top of Script 24's 95-reaction network. Chosen because it
couples directly to the Mat1a/Cbs/Cdo1 genotype differences Scripts 25/26
just surfaced (the methionine cycle's committed step and transsulfuration
arm) -- these modules are the main CONSUMERS of SAM, so if WT and ob/ob
really differ in methionine-cycle capacity, the downstream methylation
sinks are where that difference should propagate. All modules are grounded
in our own tissue metabolomics panel (Guanidinoacetate, Creatine,
Creatinine, Putrescine, Spermidine, Spermine, 5-Methylthioadenosine,
N1-Acetylspermidine, Choline, Betaine, N,N-Dimethylglycine, Sarcosine,
Nicotinamide, 1-Methylnicotinamide are all measured), and 12/21 of the
new GPR genes are in the proteomics panel (Gnmt, Chdh, Dmgdh, Sardh,
Nnmt, Nampt, Nmnat1, Oat, Srm, Gldc, Amt, Parp1 -- checked directly).

What's added:

1. GNMT made real (replaces Script 24's generic METHYLATION_SINK, which
   was `SAM_c -> SAH_c` with GPR "GNMT" but no methyl acceptor). Real
   glycine N-methyltransferase: `GLY_c + SAM_c -> SARC_c + SAH_c`, plus
   SARDH (sarcosine -> glycine, FAD-linked) closing the loop, plus the
   glycine cleavage system (GLDC/AMT, lumped) as a true glycine sink
   releasing NH3 into the (real, Script 24) urea cycle. GNMT is the
   liver's SAM overflow valve -- the single most abundant methyl-
   transferase in hepatocytes.
2. Creatine synthesis: GATM (Arg + Gly -> guanidinoacetate + ornithine)
   + GAMT (GAA + SAM -> creatine + SAH) + creatine export. Quantitatively
   the LARGEST SAM consumer in mammals; both intermediates measured.
3. Polyamine synthesis: AMD1 (SAM -> decarboxylated SAM) + ODC1
   (ornithine -> putrescine) + SRM/SMS (-> spermidine -> spermine, each
   consuming one dcSAM and releasing MTA). Requires net ornithine
   synthesis, so ORN_SYN (glutamate -> ornithine, ALDH18A1+OAT lumped,
   ATP+NADPH cost) is added too -- the urea cycle alone is a closed
   ornithine loop with no net production.
4. Betaine remethylation arm: choline uptake -> CHDH/ALDH7A1 (choline ->
   betaine, lumped) -> BHMT (betaine + homocysteine -> methionine +
   dimethylglycine) -> DMGDH (DMG -> sarcosine), joining the GNMT/SARDH
   glycine loop. Gives homocysteine a second fate besides CBS
   transsulfuration -- the real fasting-liver branch point.
5. NAD salvage + NNMT: NAD consumers (SIRT1/PARP1 lumped, NAD -> NAM) +
   NAMPT/NMNAT1 resynthesis (NAM + 2 ATP -> NAD) + NNMT (NAM + SAM ->
   1-methylnicotinamide + SAH, exported). Links the SAM economy to the
   NAD_c pool; 1-methylnicotinamide is measured.

5 new tasks (20 total): Creatine synthesis, Polyamine synthesis,
Homocysteine remethylation via betaine (forces GNMT flux with NO
methionine input -- proves the remethylation arm can sustain the
methylation cycle on choline alone), Nicotinamide clearance via NNMT,
NAD salvage cycle.

Network: Script 24's 95 reactions - 1 (METHYLATION_SINK removed) + 23 new
= 117 reactions.

Same two-stage MultiSampleFBA(lambda_reg=1.0) mechanism as Scripts 19-24.

Outputs (corneto_starvation/result/option1_sam_methylation_expansion/):
  task_flux.csv, selected_reactions.csv, task_flux_gallery.html,
  build_summary.txt, figures/union_network.{svg,png}
"""

from __future__ import annotations

import sys
from pathlib import Path

PKG_SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(PKG_SRC))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import corneto as cn  # noqa: E402
from corneto.backend import VarType  # noqa: E402
from corneto.methods import MultiSampleFBA  # noqa: E402

from corneto_mcmt.model_example import plot_task_flux_gallery  # noqa: E402

PKG_ROOT = Path(__file__).resolve().parents[1]
OUTDIR = PKG_ROOT / "result" / "option1_sam_methylation_expansion"

INPUTS = (
    "EX_GLC_IN", "EX_FA_IN", "EX_LAC_IN", "EX_ALA_IN", "EX_GLYCOGEN_IN",
    "EX_GLY_IN", "EX_GLN_IN",
    "EX_ASP_IN", "EX_LEU_IN", "EX_ILE_IN", "EX_VAL_IN",
    "EX_MET_IN", "EX_SER_IN", "EX_GLYAA_IN",
    # NEW
    "EX_CHOL_IN", "EX_NAM_IN",
)
OUTPUTS = (
    "EX_GLC_OUT", "EX_LAC_OUT", "EX_LIPID_OUT", "EX_KET_OUT",
    "EX_GLYCOGEN_OUT", "EX_HEXOSAMINE_OUT", "EX_UREA_OUT",
    "EX_TAU_OUT", "EX_GSH_OUT",
    # NEW
    "EX_CREAT_OUT", "EX_SPD_OUT", "EX_SPM_OUT", "EX_MTA_OUT", "EX_MNA_OUT",
)
BLOCKED_CYCLES = ("FUTILE_ATPASE",)


def build_full_metabolic_network() -> cn.Graph:
    """Script 24's network (minus the generic METHYLATION_SINK) plus the
    SAM methylation economy modules described in the module docstring."""
    specs = [
        # --- Exchange & Transport ---
        ("EX_GLC_IN", {}, {"GLC_c": 1.0}, 0.0, 10.0, "Transport", "SLC2A2 or SLC2A1"),
        ("EX_GLC_OUT", {"GLC_c": 1.0}, {}, 0.0, 10.0, "Transport", "SLC2A2"),
        ("EX_LAC_IN", {}, {"LAC_c": 1.0}, 0.0, 10.0, "Transport", "SLC16A1 or SLC16A7"),
        ("EX_LAC_OUT", {"LAC_c": 1.0}, {}, 0.0, 10.0, "Transport", "SLC16A1"),
        ("EX_ALA_IN", {}, {"ALA_c": 1.0}, 0.0, 10.0, "Transport", "SLC38A2 or SLC38A4"),
        ("EX_FA_IN", {}, {"FA_c": 1.0}, 0.0, 10.0, "Transport", "CD36 or SLC27A2"),
        ("EX_LIPID_OUT", {"PALM_c": 1.0}, {}, 0.0, 10.0, "Lipogenesis", "MTTP and APOB"),
        ("EX_KET_OUT", {"BHB_m": 1.0}, {}, 0.0, 10.0, "Ketogenesis", "SLC16A6 or SLC16A1"),
        ("EX_GLYCOGEN_OUT", {"GLYG_c": 1.0}, {}, 0.0, 10.0, "Glycogen", ""),
        ("EX_HEXOSAMINE_OUT", {"UDPGNAc_c": 1.0}, {}, 0.0, 10.0, "Decoy_Hexosamine", ""),
        ("ATP_MAINT", {"ATP_c": 1.0}, {"ADP_c": 1.0}, 1.0, 10.0, "Maintenance", ""),
        ("EX_GLYCOGEN_IN", {}, {"GLYG_c": 1.0}, 0.0, 10.0, "Glycogen", ""),

        # --- Upper Glycolysis & Gluconeogenesis ---
        ("GK", {"GLC_c": 1.0, "ATP_c": 1.0}, {"G6P_c": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Glycolysis", "GCK or HK1"),
        ("G6PC", {"G6P_c": 1.0}, {"GLC_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "G6PC and SLC37A4"),
        ("GPI_F", {"G6P_c": 1.0}, {"F6P_c": 1.0}, 0.0, 10.0, "Glycolysis", "GPI"),
        ("GPI_R", {"F6P_c": 1.0}, {"G6P_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "GPI"),
        ("PFK1", {"F6P_c": 1.0, "ATP_c": 1.0}, {"FBP_c": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Glycolysis", "PFKL or PFKM"),
        ("FBP1", {"FBP_c": 1.0}, {"F6P_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "FBP1"),
        ("ALDO_F", {"FBP_c": 1.0}, {"GAP_c": 2.0}, 0.0, 10.0, "Glycolysis", "ALDOB or ALDOA"),
        ("ALDO_R", {"GAP_c": 2.0}, {"FBP_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "ALDOB"),

        # --- Lower Glycolysis & Gluconeogenesis ---
        ("GAPDH_PGK_F", {"GAP_c": 1.0, "NAD_c": 1.0, "ADP_c": 1.0}, {"PGA3_c": 1.0, "NADH_c": 1.0, "ATP_c": 1.0}, 0.0, 10.0, "Glycolysis", "GAPDH and PGK1"),
        ("GAPDH_PGK_R", {"PGA3_c": 1.0, "NADH_c": 1.0, "ATP_c": 1.0}, {"GAP_c": 1.0, "NAD_c": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "GAPDH and PGK1"),
        ("ENO_F", {"PGA3_c": 1.0}, {"PEP_c": 1.0}, 0.0, 10.0, "Glycolysis", "ENO1 or ENO3"),
        ("ENO_R", {"PEP_c": 1.0}, {"PGA3_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "ENO1"),
        ("PK", {"PEP_c": 1.0, "ADP_c": 1.0}, {"PYR_c": 1.0, "ATP_c": 1.0}, 0.0, 10.0, "Glycolysis", "PKLR or PKM"),

        # --- Lactate / Alanine / Pyruvate ---
        ("LDH_F", {"PYR_c": 1.0, "NADH_c": 1.0}, {"LAC_c": 1.0, "NAD_c": 1.0}, 0.0, 10.0, "Glycolysis", "LDHA"),
        ("LDH_R", {"LAC_c": 1.0, "NAD_c": 1.0}, {"PYR_c": 1.0, "NADH_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "LDHB or LDHA"),
        ("ALT_IN", {"ALA_c": 1.0, "AKG_c": 1.0}, {"PYR_c": 1.0, "GLU_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "GPT or GPT2"),
        ("MPC", {"PYR_c": 1.0}, {"PYR_m": 1.0}, -10.0, 10.0, "Transport", "MPC1 and MPC2"),

        # --- Anaplerosis & PEPCK ---
        ("PC", {"PYR_m": 1.0, "ATP_m": 1.0}, {"OAA_m": 1.0, "ADP_m": 1.0}, 0.0, 10.0, "Gluconeogenesis", "PC"),
        ("MDH_SHUTTLE_OUT", {"OAA_m": 1.0, "NADH_m": 1.0}, {"OAA_c": 1.0, "NADH_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "MDH2 and MDH1 and SLC25A11"),
        ("PEPCK_C", {"OAA_c": 1.0, "ATP_c": 1.0}, {"PEP_c": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "PCK1"),
        ("PEPCK_M", {"OAA_m": 1.0, "ATP_m": 1.0}, {"PEP_c": 1.0, "ADP_m": 1.0}, 0.0, 10.0, "Gluconeogenesis", "PCK2"),

        # --- PPP ---
        ("G6PDH", {"G6P_c": 1.0, "NADP_c": 1.0}, {"P6G_c": 1.0, "NADPH_c": 1.0}, 0.0, 10.0, "PPP", "G6PDX and PGLS"),
        ("PGD", {"P6G_c": 1.0, "NADP_c": 1.0}, {"RU5P_c": 1.0, "NADPH_c": 1.0}, 0.0, 10.0, "PPP", "PGD"),
        ("TKT_TALDO", {"RU5P_c": 1.0}, {"F6P_c": 1.0}, 0.0, 10.0, "PPP", "TKT and TALDO1"),
        ("ME1", {"PYR_c": 1.0, "NADP_c": 1.0}, {"ACCOA_c": 1.0, "NADPH_c": 1.0}, 0.0, 10.0, "Lipogenesis", "ME1"),

        # --- PDH, Citrate Shuttle & Lipogenesis ---
        ("PDH", {"PYR_m": 1.0, "NAD_m": 1.0}, {"ACCOA_m": 1.0, "NADH_m": 1.0}, 0.0, 10.0, "TCA_OxPhos", "PDHA1 and PDHB and DLD"),
        ("CS", {"ACCOA_m": 1.0, "OAA_m": 1.0}, {"CIT_m": 1.0}, 0.0, 10.0, "TCA_OxPhos", "CS"),
        ("CL_SHUTTLE", {"CIT_m": 1.0, "ATP_c": 1.0}, {"ACCOA_c": 1.0, "OAA_m": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Lipogenesis", "SLC25A1 and ACLY"),
        ("ACC", {"ACCOA_c": 1.0, "ATP_c": 1.0}, {"MALCOA_c": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Lipogenesis", "ACACA or ACACB"),
        ("FASN", {"ACCOA_c": 1.0, "MALCOA_c": 1.0, "NADPH_c": 2.0}, {"PALM_c": 1.0, "NADP_c": 2.0}, 0.0, 10.0, "Lipogenesis", "FASN"),

        # --- Beta-Oxidation & Ketogenesis ---
        ("CPT1_CPT2", {"FA_c": 1.0, "ATP_c": 1.0}, {"FACOA_m": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Beta_Oxidation", "CPT1A and CPT2"),
        ("ACADM_STEP", {"FACOA_m": 1.0, "FAD_m": 1.0}, {"ENOYL_m": 1.0, "FADH2_m": 1.0}, 0.0, 10.0, "Beta_Oxidation", "ACADVL or ACADM or ACADS"),
        ("MTP_STEP", {"ENOYL_m": 1.0, "NAD_m": 1.0}, {"ACCOA_m": 2.0, "NADH_m": 1.0}, 0.0, 10.0, "Beta_Oxidation", "HADHA and HADHB"),
        ("HMGCS2", {"ACCOA_m": 2.0}, {"HMGCOA_m": 1.0}, 0.0, 10.0, "Ketogenesis", "HMGCS2"),
        ("HMGCL_BDH1", {"HMGCOA_m": 1.0, "NADH_m": 1.0}, {"BHB_m": 1.0, "NAD_m": 1.0}, 0.0, 10.0, "Ketogenesis", "HMGCL and BDH1"),

        # --- TCA & OxPhos ---
        ("TCA_OX", {"CIT_m": 1.0, "NAD_m": 2.0, "FAD_m": 1.0, "ADP_m": 1.0}, {"OAA_m": 1.0, "NADH_m": 2.0, "FADH2_m": 1.0, "ATP_m": 1.0}, 0.0, 10.0, "TCA_OxPhos", "IDH2 and OGDH and SDHA"),
        ("ETC_NADH", {"NADH_m": 1.0, "ADP_m": 2.0}, {"NAD_m": 1.0, "ATP_m": 2.0}, 0.0, 15.0, "TCA_OxPhos", "NDUFS1 and UQCRC1 and COX4I1"),
        ("ETC_FADH2", {"FADH2_m": 1.0, "ADP_m": 1.0}, {"FAD_m": 1.0, "ATP_m": 1.0}, 0.0, 15.0, "TCA_OxPhos", "SDHA and SDHB and COX4I1"),
        ("ANT_ATP", {"ATP_m": 1.0, "ADP_c": 1.0}, {"ATP_c": 1.0, "ADP_m": 1.0}, -15.0, 15.0, "Transport", "SLC25A4 or SLC25A5"),
        ("MAS_SHUTTLE", {"NADH_c": 1.0, "NAD_m": 1.0}, {"NAD_c": 1.0, "NADH_m": 1.0}, 0.0, 10.0, "Transport", "GOT1 and GOT2 and SLC25A12"),
        ("G3P_SHUTTLE", {"NADH_c": 1.0, "FAD_m": 1.0}, {"NAD_c": 1.0, "FADH2_m": 1.0}, 0.0, 10.0, "Transport", "GPD1 and GPD2"),

        # --- Decoy Branches ---
        ("GLYG_SYN", {"G6P_c": 1.0, "ATP_c": 1.0}, {"GLYG_c": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Glycogen", "GYS2 and UGP2"),
        ("GLYG_DEG", {"GLYG_c": 1.0}, {"G6P_c": 1.0}, 0.0, 10.0, "Glycogen", "PYGL"),
        ("GFAT_HBP", {"F6P_c": 1.0, "ATP_c": 1.0}, {"UDPGNAc_c": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Decoy_Hexosamine", "GFPT1 and UAP1"),
        ("FUTILE_ATPASE", {"ATP_c": 1.0}, {"ADP_c": 1.0}, 0.0, 10.0, "Decoy_Futile", "ENTPD5"),
        ("NNT_FORWARD", {"NADH_m": 1.0, "NADP_c": 1.0}, {"NAD_m": 1.0, "NADPH_c": 1.0}, 0.0, 10.0, "Decoy_Futile", "NNT"),
        ("NNT_REVERSE", {"NADPH_c": 1.0, "NAD_m": 1.0}, {"NADP_c": 1.0, "NADH_m": 1.0}, 0.0, 10.0, "Decoy_Futile", "NNT_REV"),

        # --- Glutamate handling ---
        ("GLUD", {"GLU_c": 1.0, "NAD_c": 1.0}, {"AKG_c": 1.0, "NH3_c": 1.0, "NADH_c": 1.0}, 0.0, 10.0, "Ureagenesis", "GLUD1"),
        ("EX_UREA_OUT", {"UREA_c": 1.0}, {}, 0.0, 10.0, "Ureagenesis", "SLC14A1"),

        # --- Glycerol gluconeogenesis ---
        ("EX_GLY_IN", {}, {"GLYR_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "AQP9"),
        ("GLY_KIN", {"GLYR_c": 1.0, "ATP_c": 1.0}, {"G3P_c": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "GK"),
        ("G3PDH_CYT", {"G3P_c": 1.0, "NAD_c": 1.0}, {"DHAP_c": 1.0, "NADH_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "GPD1"),
        ("TPI", {"DHAP_c": 1.0}, {"GAP_c": 1.0}, -10.0, 10.0, "Gluconeogenesis", "TPI1"),

        # --- Glutamine ureagenesis entry ---
        ("EX_GLN_IN", {}, {"GLN_c": 1.0}, 0.0, 10.0, "Ureagenesis", "SLC38A3 or SLC1A5"),
        ("GLS", {"GLN_c": 1.0}, {"GLU_c": 1.0, "NH3_c": 1.0}, 0.0, 10.0, "Ureagenesis", "GLS2"),

        ("AKG_TO_OAA", {"AKG_c": 1.0, "NAD_m": 1.0, "FAD_m": 1.0, "ADP_m": 1.0}, {"OAA_c": 1.0, "NADH_m": 1.0, "FADH2_m": 1.0, "ATP_m": 1.0}, 0.0, 10.0, "TCA_OxPhos", "OGDH and SUCLA2 and SDHA and FH and MDH2"),

        # --- Full urea cycle (Script 24) ---
        ("EX_ASP_IN", {}, {"ASP_c": 1.0}, 0.0, 10.0, "Ureagenesis", "SLC1A3 or SLC25A13"),
        ("CPS1_RXN", {"NH3_c": 1.0, "ATP_m": 2.0}, {"CARBP_m": 1.0, "ADP_m": 2.0}, 0.0, 10.0, "Ureagenesis", "CPS1"),
        ("OTC_RXN", {"CARBP_m": 1.0, "ORN_m": 1.0}, {"CIT_m_UREA": 1.0}, 0.0, 10.0, "Ureagenesis", "OTC"),
        ("CIT_UREA_EXPORT", {"CIT_m_UREA": 1.0}, {"CIT_c_UREA": 1.0}, 0.0, 10.0, "Ureagenesis", "SLC25A15"),
        ("ORN_IMPORT", {"ORN_c": 1.0}, {"ORN_m": 1.0}, 0.0, 10.0, "Ureagenesis", "SLC25A15"),
        ("ASS1_RXN", {"CIT_c_UREA": 1.0, "ASP_c": 1.0, "ATP_c": 1.0}, {"ARGSUC_c": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Ureagenesis", "ASS1"),
        ("ASL_RXN", {"ARGSUC_c": 1.0}, {"ARG_c": 1.0, "FUM_c": 1.0}, 0.0, 10.0, "Ureagenesis", "ASL"),
        ("ARG1_RXN", {"ARG_c": 1.0}, {"ORN_c": 1.0, "UREA_c": 1.0}, 0.0, 10.0, "Ureagenesis", "ARG1"),
        ("FUM_TO_OAA", {"FUM_c": 1.0, "NAD_c": 1.0}, {"OAA_c": 1.0, "NADH_c": 1.0}, 0.0, 10.0, "Ureagenesis", "FH and MDH1"),

        # --- BCAA catabolism (Script 24) ---
        ("EX_LEU_IN", {}, {"LEU_c": 1.0}, 0.0, 10.0, "AA_Catabolism", "SLC7A5 or SLC3A2"),
        ("LEU_OX", {"LEU_c": 1.0, "NAD_m": 1.0}, {"ACCOA_m": 2.0, "NADH_m": 1.0, "NH3_c": 1.0}, 0.0, 10.0, "AA_Catabolism", "BCAT2 and BCKDHA and BCKDHB and IVD and HMGCL"),
        ("EX_ILE_IN", {}, {"ILE_c": 1.0}, 0.0, 10.0, "AA_Catabolism", "SLC7A5 or SLC3A2"),
        ("ILE_OX", {"ILE_c": 1.0, "NAD_m": 1.0, "FAD_m": 1.0}, {"ACCOA_m": 1.0, "OAA_m": 1.0, "NADH_m": 1.0, "FADH2_m": 1.0, "NH3_c": 1.0}, 0.0, 10.0, "AA_Catabolism", "BCAT2 and BCKDHA and BCKDHB and ACADSB"),
        ("EX_VAL_IN", {}, {"VAL_c": 1.0}, 0.0, 10.0, "AA_Catabolism", "SLC7A5 or SLC3A2"),
        ("VAL_OX", {"VAL_c": 1.0, "NAD_m": 1.0, "FAD_m": 1.0}, {"OAA_m": 1.0, "NADH_m": 1.0, "FADH2_m": 1.0, "NH3_c": 1.0}, 0.0, 10.0, "AA_Catabolism", "BCAT2 and BCKDHA and BCKDHB and MUT"),

        # --- Methionine cycle / transsulfuration (Script 24) ---
        ("EX_MET_IN", {}, {"MET_c": 1.0}, 0.0, 10.0, "Methionine_Cycle", "SLC7A5 or SLC3A2"),
        ("MAT1A_RXN", {"MET_c": 1.0, "ATP_c": 1.0}, {"SAM_c": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Methionine_Cycle", "MAT1A"),
        ("SAHH_RXN", {"SAH_c": 1.0}, {"HCY_c": 1.0}, 0.0, 10.0, "Methionine_Cycle", "AHCY"),
        ("EX_SER_IN", {}, {"SER_c": 1.0}, 0.0, 10.0, "Methionine_Cycle", "SLC1A4 or SLC1A5"),
        ("CBS_RXN", {"HCY_c": 1.0, "SER_c": 1.0}, {"CYSTA_c": 1.0}, 0.0, 10.0, "Methionine_Cycle", "CBS"),
        ("CTH_RXN", {"CYSTA_c": 1.0}, {"CYS_c": 1.0}, 0.0, 10.0, "Methionine_Cycle", "CTH"),
        ("CYS_TO_TAU", {"CYS_c": 1.0, "NAD_c": 1.0}, {"TAU_c": 1.0, "NADH_c": 1.0}, 0.0, 10.0, "Methionine_Cycle", "CDO1 and CSAD"),
        ("EX_TAU_OUT", {"TAU_c": 1.0}, {}, 0.0, 10.0, "Methionine_Cycle", "SLC6A6"),
        ("EX_GLYAA_IN", {}, {"GLY_c": 1.0}, 0.0, 10.0, "Methionine_Cycle", "SLC6A9"),
        ("GSH_SYN", {"CYS_c": 1.0, "GLU_c": 1.0, "GLY_c": 1.0, "ATP_c": 2.0}, {"GSH_c": 1.0, "ADP_c": 2.0}, 0.0, 10.0, "Methionine_Cycle", "GCLC and GCLM and GSS"),
        ("EX_GSH_OUT", {"GSH_c": 1.0}, {}, 0.0, 10.0, "Methionine_Cycle", "ABCC2"),

        # ============ NEW: real GNMT + glycine handling ============
        # Replaces the generic METHYLATION_SINK (SAM_c -> SAH_c): GNMT is
        # the liver's actual SAM overflow valve and needs a real methyl
        # acceptor (glycine). SARDH closes the sarcosine->glycine loop;
        # the glycine cleavage system is the true net glycine sink,
        # releasing NH3 into the urea cycle.
        ("GNMT_RXN", {"GLY_c": 1.0, "SAM_c": 1.0}, {"SARC_c": 1.0, "SAH_c": 1.0}, 0.0, 10.0, "Methylation", "GNMT"),
        ("SARDH_RXN", {"SARC_c": 1.0, "FAD_m": 1.0}, {"GLY_c": 1.0, "FADH2_m": 1.0}, 0.0, 10.0, "Methylation", "SARDH"),
        ("GLY_CLEAV", {"GLY_c": 1.0, "NAD_m": 1.0}, {"NH3_c": 1.0, "NADH_m": 1.0}, 0.0, 10.0, "Methylation", "GLDC and AMT"),

        # ============ NEW: creatine synthesis ============
        ("GATM_RXN", {"ARG_c": 1.0, "GLY_c": 1.0}, {"GAA_c": 1.0, "ORN_c": 1.0}, 0.0, 10.0, "Creatine", "GATM"),
        ("GAMT_RXN", {"GAA_c": 1.0, "SAM_c": 1.0}, {"CREAT_c": 1.0, "SAH_c": 1.0}, 0.0, 10.0, "Creatine", "GAMT"),
        ("EX_CREAT_OUT", {"CREAT_c": 1.0}, {}, 0.0, 10.0, "Creatine", ""),

        # ============ NEW: polyamine synthesis ============
        # ORN_SYN (ALDH18A1+OAT lumped, glutamate -> ornithine) provides
        # the NET ornithine the polyamine drain requires -- the urea cycle
        # alone is a closed ornithine loop with no net production.
        ("ORN_SYN", {"GLU_c": 1.0, "ATP_c": 1.0, "NADPH_c": 1.0}, {"ORN_c": 1.0, "ADP_c": 1.0, "NADP_c": 1.0}, 0.0, 10.0, "Polyamine", "ALDH18A1 and OAT"),
        ("SAMDC_RXN", {"SAM_c": 1.0}, {"DCSAM_c": 1.0}, 0.0, 10.0, "Polyamine", "AMD1"),
        ("ODC_RXN", {"ORN_c": 1.0}, {"PUTR_c": 1.0}, 0.0, 10.0, "Polyamine", "ODC1"),
        ("SRM_RXN", {"PUTR_c": 1.0, "DCSAM_c": 1.0}, {"SPD_c": 1.0, "MTA_c": 1.0}, 0.0, 10.0, "Polyamine", "SRM"),
        ("SMS_RXN", {"SPD_c": 1.0, "DCSAM_c": 1.0}, {"SPM_c": 1.0, "MTA_c": 1.0}, 0.0, 10.0, "Polyamine", "SMS"),
        ("EX_SPD_OUT", {"SPD_c": 1.0}, {}, 0.0, 10.0, "Polyamine", ""),
        ("EX_SPM_OUT", {"SPM_c": 1.0}, {}, 0.0, 10.0, "Polyamine", ""),
        ("EX_MTA_OUT", {"MTA_c": 1.0}, {}, 0.0, 10.0, "Polyamine", ""),

        # ============ NEW: betaine remethylation arm ============
        ("EX_CHOL_IN", {}, {"CHOL_c": 1.0}, 0.0, 10.0, "Remethylation", "SLC44A1"),
        ("CHDH_RXN", {"CHOL_c": 1.0, "NAD_c": 1.0}, {"BET_c": 1.0, "NADH_c": 1.0}, 0.0, 10.0, "Remethylation", "CHDH and ALDH7A1"),
        ("BHMT_RXN", {"BET_c": 1.0, "HCY_c": 1.0}, {"MET_c": 1.0, "DMG_c": 1.0}, 0.0, 10.0, "Remethylation", "BHMT"),
        ("DMGDH_RXN", {"DMG_c": 1.0, "FAD_m": 1.0}, {"SARC_c": 1.0, "FADH2_m": 1.0}, 0.0, 10.0, "Remethylation", "DMGDH"),

        # ============ NEW: NAD salvage + NNMT ============
        ("EX_NAM_IN", {}, {"NAM_c": 1.0}, 0.0, 10.0, "NAD_Salvage", ""),
        ("NNMT_RXN", {"NAM_c": 1.0, "SAM_c": 1.0}, {"MNA_c": 1.0, "SAH_c": 1.0}, 0.0, 10.0, "NAD_Salvage", "NNMT"),
        ("EX_MNA_OUT", {"MNA_c": 1.0}, {}, 0.0, 10.0, "NAD_Salvage", ""),
        ("NAMPT_NMNAT", {"NAM_c": 1.0, "ATP_c": 2.0}, {"NAD_c": 1.0, "ADP_c": 2.0}, 0.0, 10.0, "NAD_Salvage", "NAMPT and NMNAT1"),
        ("NAD_CONSUM", {"NAD_c": 1.0}, {"NAM_c": 1.0}, 0.0, 10.0, "NAD_Salvage", "SIRT1 or PARP1"),
    ]

    graph = cn.Graph()
    for rxn_id, reactants, products, lb, ub, pathway, gpr in specs:
        graph.add_edge(
            reactants, products, id=rxn_id, name=rxn_id,
            default_lb=float(lb), default_ub=float(ub),
            pathway=pathway, GPR=gpr, reversible=bool(lb < 0),
        )
    return graph


def build_task_bounds() -> dict[str, dict[str, tuple[float, float]]]:
    base_bounds = {reaction: (0.0, 0.0) for reaction in (*INPUTS, *OUTPUTS, *BLOCKED_CYCLES)}
    base_bounds["ATP_MAINT"] = (1.0, 1.0)
    return {
        # --- Script 24's 15 tasks, unchanged in intent. The transsulfuration
        # and glutathione tasks still work with GNMT replacing the generic
        # sink: GNMT's glycine is fully recycled through SARDH (GNMT makes
        # 1 SARC per GLY consumed, SARDH converts it back), so no new
        # glycine input is needed there.
        "Gluconeogenesis from lactate": {
            **base_bounds, "EX_LAC_IN": (0.0, 10.0), "EX_GLC_OUT": (1.0, 10.0),
        },
        "Gluconeogenesis from alanine": {
            **base_bounds, "EX_ALA_IN": (0.0, 10.0), "EX_GLC_OUT": (1.0, 10.0),
            "EX_UREA_OUT": (0.0, 10.0), "EX_ASP_IN": (0.0, 10.0),
        },
        "Glycogenolysis": {
            **base_bounds, "EX_GLYCOGEN_IN": (0.0, 10.0), "EX_GLC_OUT": (1.0, 10.0),
        },
        "Ketogenesis from fatty acid": {
            **base_bounds, "EX_FA_IN": (0.0, 10.0), "EX_KET_OUT": (1.0, 10.0),
        },
        "Fatty acid oxidation for energy": {
            **base_bounds, "EX_FA_IN": (0.0, 10.0), "ATP_MAINT": (5.0, 10.0),
        },
        "Ureagenesis from alanine": {
            **base_bounds, "EX_ALA_IN": (0.0, 10.0), "EX_UREA_OUT": (1.0, 10.0),
            "EX_ASP_IN": (0.0, 10.0),
        },
        "Gluconeogenesis from glycerol": {
            **base_bounds, "EX_GLY_IN": (0.0, 10.0), "EX_GLC_OUT": (1.0, 10.0),
        },
        "Ureagenesis from glutamine": {
            **base_bounds, "EX_GLN_IN": (0.0, 10.0), "EX_UREA_OUT": (1.0, 10.0),
            "EX_GLC_OUT": (0.0, 10.0), "EX_ASP_IN": (0.0, 10.0),
        },
        "Combined fasted state": {
            **base_bounds,
            "EX_FA_IN": (0.0, 10.0), "EX_LAC_IN": (0.0, 10.0), "EX_ALA_IN": (0.0, 10.0),
            "EX_GLY_IN": (0.0, 10.0), "EX_GLN_IN": (0.0, 10.0), "EX_ASP_IN": (0.0, 10.0),
            "EX_GLC_OUT": (1.0, 1.0), "EX_KET_OUT": (1.0, 1.0), "EX_UREA_OUT": (1.0, 1.0),
            "ATP_MAINT": (2.0, 2.0),
        },
        "Glucose to glycogen (storage)": {
            **base_bounds, "EX_GLC_IN": (0.0, 10.0), "EX_GLYCOGEN_OUT": (1.0, 10.0),
        },
        "Glucose to lipid (lipogenesis)": {
            **base_bounds, "EX_GLC_IN": (0.0, 10.0), "EX_LIPID_OUT": (1.0, 10.0),
        },
        "Combined fed state": {
            **base_bounds, "EX_GLC_IN": (0.0, 10.0),
            "EX_GLYCOGEN_OUT": (1.0, 1.0), "EX_LIPID_OUT": (1.0, 1.0),
            "ATP_MAINT": (2.0, 2.0),
        },
        "BCAA oxidation for energy": {
            **base_bounds,
            "EX_LEU_IN": (0.0, 10.0), "EX_ILE_IN": (0.0, 10.0), "EX_VAL_IN": (0.0, 10.0),
            "EX_UREA_OUT": (1.0, 10.0), "EX_ASP_IN": (0.0, 10.0),
            "ATP_MAINT": (2.0, 10.0),
        },
        "Transsulfuration: methionine to taurine": {
            **base_bounds,
            "EX_MET_IN": (0.0, 10.0), "EX_SER_IN": (0.0, 10.0), "EX_TAU_OUT": (1.0, 10.0),
        },
        "Glutathione synthesis": {
            **base_bounds,
            "EX_MET_IN": (0.0, 10.0), "EX_SER_IN": (0.0, 10.0), "EX_GLYAA_IN": (0.0, 10.0),
            "EX_GLN_IN": (0.0, 10.0), "EX_GSH_OUT": (1.0, 10.0),
            "EX_UREA_OUT": (0.0, 10.0), "EX_ASP_IN": (0.0, 10.0),
        },

        # --- NEW tasks exercising the SAM methylation economy ---
        # Net reaction: Gly + (SAM methyl) -> guanidinoacetate -> creatine.
        # Arg is regenerated by the urea cycle (ORN -> ... -> ARG), which
        # needs NH3 (from Gln via GLS) and Asp (ASS1); urea is an optional
        # co-product, glucose the carbon side-exit for Gln's skeleton.
        # GAMT's SAH byproduct ends as homocysteine, whose disposal is
        # obligate (same lesson as every nitrogen task): the
        # transsulfuration arm (Ser in, taurine out) is opened for it.
        "Creatine synthesis": {
            **base_bounds,
            "EX_MET_IN": (0.0, 10.0), "EX_GLYAA_IN": (0.0, 10.0),
            "EX_ASP_IN": (0.0, 10.0), "EX_GLN_IN": (0.0, 10.0),
            "EX_SER_IN": (0.0, 10.0), "EX_TAU_OUT": (0.0, 10.0),
            "EX_UREA_OUT": (0.0, 10.0), "EX_GLC_OUT": (0.0, 10.0),
            "EX_CREAT_OUT": (1.0, 10.0),
        },
        # Spermine export forced so BOTH SRM and SMS must run (spermidine
        # is the intermediate). Ornithine comes from ORN_SYN (glutamate),
        # dcSAM from methionine; glucose opens the PPP for ORN_SYN's NADPH.
        "Polyamine synthesis": {
            **base_bounds,
            "EX_MET_IN": (0.0, 10.0), "EX_GLN_IN": (0.0, 10.0),
            "EX_ASP_IN": (0.0, 10.0), "EX_GLC_IN": (0.0, 10.0),
            "EX_UREA_OUT": (0.0, 10.0), "EX_GLC_OUT": (0.0, 10.0),
            "EX_MTA_OUT": (0.0, 10.0), "EX_SPM_OUT": (1.0, 10.0),
        },
        # NO methionine input: forces GNMT to run on a methylation cycle
        # sustained entirely by choline-derived betaine remethylating
        # homocysteine back to methionine. The net glycine surplus (DMG ->
        # sarcosine -> glycine exceeds GNMT's consumption) exits via the
        # glycine cleavage system into the urea cycle.
        "Homocysteine remethylation via betaine": {
            **base_bounds,
            "EX_CHOL_IN": (0.0, 10.0), "GNMT_RXN": (1.0, 10.0),
            "EX_ASP_IN": (0.0, 10.0), "EX_UREA_OUT": (0.0, 10.0),
            "EX_GLC_OUT": (0.0, 10.0),
        },
        # Methyl group arrives on SAM (from Met) and leaves on
        # 1-methylnicotinamide; homocysteine disposal may use either arm
        # (transsulfuration via Ser, or remethylation via choline).
        "Nicotinamide clearance via NNMT": {
            **base_bounds,
            "EX_NAM_IN": (0.0, 10.0), "EX_MET_IN": (0.0, 10.0),
            "EX_SER_IN": (0.0, 10.0), "EX_CHOL_IN": (0.0, 10.0),
            "EX_TAU_OUT": (0.0, 10.0), "EX_ASP_IN": (0.0, 10.0),
            "EX_UREA_OUT": (0.0, 10.0), "EX_GLC_OUT": (0.0, 10.0),
            "EX_MNA_OUT": (1.0, 10.0),
        },
        # NAD consumer flux (SIRT1/PARP1) forced; NAMPT/NMNAT1 must
        # resynthesize NAD from the released nicotinamide (2 ATP/turn),
        # with fatty acid oxidation covering the energy cost.
        "NAD salvage cycle": {
            **base_bounds,
            "NAD_CONSUM": (1.0, 10.0), "EX_FA_IN": (0.0, 10.0),
            "EX_KET_OUT": (0.0, 10.0),
        },
    }


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "figures").mkdir(exist_ok=True)

    graph = build_full_metabolic_network()
    print(f"Full network: {graph.num_edges} reactions")

    task_bounds = build_task_bounds()
    print(f"Tasks ({len(task_bounds)}):", ", ".join(task_bounds))

    task_method = MultiSampleFBA(lambda_reg=1.0, beta_reg=0.0)
    task_problem = task_method.build_many(graph, reaction_bounds=task_bounds)
    task_problem.add_objective(
        task_problem.expr.edge_has_flux_OR.sum(), name="explicit_minimum_network_size"
    )
    first_solve = task_problem.solve(solver="GUROBI", verbosity=1)
    status_1 = str(first_solve.status).lower()
    print(f"Stage 1 (minimize union) solve status: {status_1}")
    if status_1 != "optimal":
        raise RuntimeError(f"The tasks could not be satisfied: {first_solve.status}")

    minimum_size = int(np.count_nonzero(np.asarray(task_problem.expr.edge_has_flux_OR.value) > 0.5))
    flow = task_problem.expr.flow
    absolute_flow = task_method.backend.Variable(
        name="task_absolute_flux", shape=flow.shape, lb=0.0, vartype=VarType.CONTINUOUS
    )
    task_problem.add_constraints([
        task_problem.expr.edge_has_flux_OR.sum() == minimum_size,
        absolute_flow >= flow,
        absolute_flow >= -flow,
    ])
    task_problem.add_objective(absolute_flow.sum(), name="parsimonious_task_flux")
    second_solve = task_problem.solve(solver="GUROBI", verbosity=1)
    status_2 = str(second_solve.status).lower()
    print(f"Stage 2 (minimize flux at fixed size) solve status: {status_2}")
    if status_2 != "optimal":
        raise RuntimeError(f"Flux refinement failed: {second_solve.status}")

    reaction_ids = list(graph.get_attr_from_edges("id"))
    task_flux = pd.DataFrame(
        np.asarray(task_problem.expr.flow.value), index=reaction_ids, columns=task_bounds
    )
    task_flux.to_csv(OUTDIR / "task_flux.csv")

    selected_reactions = np.asarray(task_problem.expr.edge_has_flux_OR.value).ravel() > 0.5
    assert selected_reactions.sum() == minimum_size
    assert np.array_equal(selected_reactions, (np.abs(task_flux.to_numpy()) > 1e-6).any(axis=1))

    print(f"Smallest joint network: {minimum_size} of {graph.num_edges} reactions")

    selected_df = pd.DataFrame({
        "reaction": np.asarray(reaction_ids)[selected_reactions],
    }).set_index("reaction")
    selected_df.to_csv(OUTDIR / "selected_reactions.csv")

    html = plot_task_flux_gallery(
        graph, task_flux.to_numpy(), list(task_bounds), selected_reactions,
    )
    (OUTDIR / "task_flux_gallery.html").write_text(f"<html><body>{html}</body></html>")

    summary = (
        f"Full network: {graph.num_edges} reactions total "
        f"(Script 24's 95 - 1 [METHYLATION_SINK replaced by real GNMT] + 23 new)\n"
        f"Tasks ({len(task_bounds)}): {', '.join(task_bounds)}\n"
        f"Stage 1 (minimize union) solve status: {status_1}\n"
        f"Stage 2 (minimize flux at fixed size) solve status: {status_2}\n"
        f"Smallest joint network: {minimum_size}/{graph.num_edges} reactions\n"
    )
    (OUTDIR / "build_summary.txt").write_text(summary)
    print("\n" + summary)


if __name__ == "__main__":
    main()
