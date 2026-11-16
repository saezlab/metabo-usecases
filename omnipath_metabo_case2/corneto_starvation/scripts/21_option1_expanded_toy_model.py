"""Option 1 (metabolic-task network curation), staying faithful to the
tutorial's approach (small hand-curated network + task feasibility,
NEVER the full Mouse-GEM -- confirmed directly: the notebook's own text
says "Option 2 uses the entire Mouse-GEM instead", meaning Option 1 is
specifically the small-network method by design, not a scale we chose to
avoid) -- but the tutorial's bare 58-reaction network (Script 19/20) is
biochemically too thin to represent several real fasting-liver functions.
This script EXTENDS that same network (same reaction-spec style, same
`cn.Graph()` construction) with real biochemistry to close those gaps,
rather than patching bounds around them.

What's added, and why (not arbitrary patches -- each fixes a specific gap
found in Script 20):

1. ALT_IN corrected to real biochemistry. The tutorial's version
   (`Ala_c -> Pyr_c + NADH_c`) is not how alanine aminotransferase (ALT/
   GPT) actually works -- it's a transaminase (PLP cofactor), not an
   oxidoreductase, so it should never produce NADH. Real reaction:
   `Ala_c + AKG_c -> Pyr_c + Glu_c` (alpha-ketoglutarate / glutamate
   couple, no nicotinamide cofactor).
2. GLUD1 (glutamate dehydrogenase) added: `Glu_c + NAD_c -> AKG_c + NH3_c
   + NADH_c`. This regenerates AKG_c for ALT_IN (closed loop, no external
   AKG pool needed) AND is the REAL source of the NADH gluconeogenesis
   needs from amino acid nitrogen disposal -- replacing the tutorial's
   biochemically-wrong shortcut with the actual pathway. This symmetrically
   fixes Script 20's finding (alanine route had no NAD_c sink to match
   GAPDH_PGK_R's production, the lactate route's LDH_R did): GLUD1 now
   consumes exactly 1 NAD_c per alanine processed, mirroring LDH_R.
3. Urea cycle, lumped into one reaction (same simplification style the
   network already uses for TCA_OX/ETC): `NH3_c + 2 ATP_c -> Urea_c +
   2 ADP_c` (2 ATP/urea is the real stoichiometric cost of carbamoyl
   phosphate synthetase + the full cycle) + `EX_UREA_OUT`. Not
   representable at all in Script 20's network (no urea cycle reactions
   existed) -- amino-acid nitrogen disposal is a core fasting-liver
   function, and now follows directly from fix #2's NH3_c byproduct.
4. `EX_GLYCOGEN_IN` added (external glycogen store -> GLYG_c). The
   tutorial's network only has glucose<->glycogen CYCLING reactions
   (GLYG_SYN/GLYG_DEG), no exchange for a pre-existing store to deplete --
   standard steady-state FBA has no accumulation/depletion dynamic, so an
   external "store" is modeled as a boundary exchange, exactly like every
   other substrate here. Also un-blocks GLYG_DEG (Script 19/20 blocked it
   as a "substrate cycle" -- it is one only when GLYG_SYN runs too; opening
   the external store on one side removes that ambiguity).

5. Glycerol gluconeogenesis added (`EX_GLY_IN -> GLY_KIN -> G3PDH_CYT ->
   TPI -> GAP_c`, merging into the existing ALDO_R/gluconeogenic chain).
   Glycerol (from adipose lipolysis during fasting) is a major
   gluconeogenic substrate alongside lactate/alanine -- previously absent
   entirely, not just infeasible.
6. Glutamine ureagenesis added (`EX_GLN_IN -> GLS -> GLU_c`, reusing the
   existing GLU_c -> GLUD -> NH3_c -> UREA_SYN chain from fix #2/#3).
   Glutamine is the other major amino-acid nitrogen carrier for fasting
   nitrogen disposal, alongside alanine.

Tasks (12: 9 fasting + 3 fed/T0-state, all representable with the above
fixes and -- for the fed-state tasks -- no new reactions at all, just
tasking machinery the original 58 already had but never exercised):
  Gluconeogenesis from lactate     (Cori cycle)
  Gluconeogenesis from alanine     (glucose-alanine cycle)
  Gluconeogenesis from glycerol    (adipose lipolysis byproduct)
  Glycogenolysis                   (early-fasting glycogen store -> glucose)
  Ketogenesis from fatty acid      (liver's main fasting fuel export)
  Fatty acid oxidation for energy  (fat alone meeting elevated ATP demand)
  Ureagenesis from alanine         (amino-acid nitrogen disposal)
  Ureagenesis from glutamine       (amino-acid nitrogen disposal)
  Combined fasted state            (FA + lactate + alanine + glycerol +
                                     glutamine all available simultaneously
                                     -> glucose + ketone + urea all at once)
  Glucose to glycogen (storage)    -- T0/fed state
  Glucose to lipid (lipogenesis)   -- T0/fed state
  Combined fed state               -- T0/fed state: glucose -> glycogen
                                       storage + de novo lipogenesis + ATP
                                       maintenance simultaneously

Urea cycle is still simplified to a single nitrogen donor (NH3_c only);
the real cycle's second nitrogen (from aspartate, releasing fumarate back
to the TCA cycle) is not modeled -- a candidate for a further round.

Still NOT representable (would need the full Mouse-GEM, not attempted
here): any pathway this hand-curated network doesn't name at all (e.g.
branched-chain amino acid catabolism, individual fatty acid species,
cholesterol/bile acid synthesis).

Same two-stage MultiSampleFBA(lambda_reg=1.0) mechanism as Scripts 19/20.

Outputs (corneto_starvation/result/option1_expanded_toy_model/): same set
as Scripts 19/20.
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
OUTDIR = PKG_ROOT / "result" / "option1_expanded_toy_model"

INPUTS = ("EX_GLC_IN", "EX_FA_IN", "EX_LAC_IN", "EX_ALA_IN", "EX_GLYCOGEN_IN", "EX_GLY_IN", "EX_GLN_IN")
OUTPUTS = (
    "EX_GLC_OUT", "EX_LAC_OUT", "EX_LIPID_OUT", "EX_KET_OUT",
    "EX_GLYCOGEN_OUT", "EX_HEXOSAMINE_OUT", "EX_UREA_OUT",
)
# GLYG_DEG no longer blocked (see module docstring, fix #4); FUTILE_ATPASE
# stays blocked (genuinely a pointless substrate cycle, no task needs it).
BLOCKED_CYCLES = ("FUTILE_ATPASE",)


def build_expanded_core_metabolic_network() -> cn.Graph:
    """Script 20's 58-reaction network (same specs, copied verbatim from
    corneto_mcmt.flux_benchmark.build_core_metabolic_network) plus the four
    biochemistry fixes described in the module docstring."""
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
        # --- NEW: external glycogen store, for a real glycogenolysis task ---
        ("EX_GLYCOGEN_IN", {}, {"GLYG_c": 1.0}, 0.0, 10.0, "Glycogen", ""),

        # --- Upper Glycolysis & Gluconeogenesis (with substrate cycle pairs) ---
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
        # FIXED: real ALT biochemistry (transamination, no NAD(H) involved).
        ("ALT_IN", {"ALA_c": 1.0, "AKG_c": 1.0}, {"PYR_c": 1.0, "GLU_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "GPT or GPT2"),
        ("MPC", {"PYR_c": 1.0}, {"PYR_m": 1.0}, -10.0, 10.0, "Transport", "MPC1 and MPC2"),

        # --- Anaplerosis & PEPCK (Gluconeogenesis commit) ---
        ("PC", {"PYR_m": 1.0, "ATP_m": 1.0}, {"OAA_m": 1.0, "ADP_m": 1.0}, 0.0, 10.0, "Gluconeogenesis", "PC"),
        ("MDH_SHUTTLE_OUT", {"OAA_m": 1.0, "NADH_m": 1.0}, {"OAA_c": 1.0, "NADH_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "MDH2 and MDH1 and SLC25A11"),
        ("PEPCK_C", {"OAA_c": 1.0, "ATP_c": 1.0}, {"PEP_c": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "PCK1"),
        ("PEPCK_M", {"OAA_m": 1.0, "ATP_m": 1.0}, {"PEP_c": 1.0, "ADP_m": 1.0}, 0.0, 10.0, "Gluconeogenesis", "PCK2"),

        # --- Pentose Phosphate Pathway (NADPH supply for Lipogenesis) ---
        ("G6PDH", {"G6P_c": 1.0, "NADP_c": 1.0}, {"P6G_c": 1.0, "NADPH_c": 1.0}, 0.0, 10.0, "PPP", "G6PDX and PGLS"),
        ("PGD", {"P6G_c": 1.0, "NADP_c": 1.0}, {"RU5P_c": 1.0, "NADPH_c": 1.0}, 0.0, 10.0, "PPP", "PGD"),
        ("TKT_TALDO", {"RU5P_c": 1.0}, {"F6P_c": 1.0}, 0.0, 10.0, "PPP", "TKT and TALDO1"),
        ("ME1", {"PYR_c": 1.0, "NADP_c": 1.0}, {"ACCOA_c": 1.0, "NADPH_c": 1.0}, 0.0, 10.0, "Lipogenesis", "ME1"),

        # --- PDH, Citrate Shuttle & De Novo Lipogenesis ---
        ("PDH", {"PYR_m": 1.0, "NAD_m": 1.0}, {"ACCOA_m": 1.0, "NADH_m": 1.0}, 0.0, 10.0, "TCA_OxPhos", "PDHA1 and PDHB and DLD"),
        ("CS", {"ACCOA_m": 1.0, "OAA_m": 1.0}, {"CIT_m": 1.0}, 0.0, 10.0, "TCA_OxPhos", "CS"),
        ("CL_SHUTTLE", {"CIT_m": 1.0, "ATP_c": 1.0}, {"ACCOA_c": 1.0, "OAA_m": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Lipogenesis", "SLC25A1 and ACLY"),
        ("ACC", {"ACCOA_c": 1.0, "ATP_c": 1.0}, {"MALCOA_c": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Lipogenesis", "ACACA or ACACB"),
        ("FASN", {"ACCOA_c": 1.0, "MALCOA_c": 1.0, "NADPH_c": 2.0}, {"PALM_c": 1.0, "NADP_c": 2.0}, 0.0, 10.0, "Lipogenesis", "FASN"),

        # --- Mitochondrial Beta-Oxidation & Ketogenesis ---
        ("CPT1_CPT2", {"FA_c": 1.0, "ATP_c": 1.0}, {"FACOA_m": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Beta_Oxidation", "CPT1A and CPT2"),
        ("ACADM_STEP", {"FACOA_m": 1.0, "FAD_m": 1.0}, {"ENOYL_m": 1.0, "FADH2_m": 1.0}, 0.0, 10.0, "Beta_Oxidation", "ACADVL or ACADM or ACADS"),
        ("MTP_STEP", {"ENOYL_m": 1.0, "NAD_m": 1.0}, {"ACCOA_m": 2.0, "NADH_m": 1.0}, 0.0, 10.0, "Beta_Oxidation", "HADHA and HADHB"),
        ("HMGCS2", {"ACCOA_m": 2.0}, {"HMGCOA_m": 1.0}, 0.0, 10.0, "Ketogenesis", "HMGCS2"),
        ("HMGCL_BDH1", {"HMGCOA_m": 1.0, "NADH_m": 1.0}, {"BHB_m": 1.0, "NAD_m": 1.0}, 0.0, 10.0, "Ketogenesis", "HMGCL and BDH1"),

        # --- TCA Cycle & Oxidative Phosphorylation ---
        ("TCA_OX", {"CIT_m": 1.0, "NAD_m": 2.0, "FAD_m": 1.0, "ADP_m": 1.0}, {"OAA_m": 1.0, "NADH_m": 2.0, "FADH2_m": 1.0, "ATP_m": 1.0}, 0.0, 10.0, "TCA_OxPhos", "IDH2 and OGDH and SDHA"),
        ("ETC_NADH", {"NADH_m": 1.0, "ADP_m": 2.0}, {"NAD_m": 1.0, "ATP_m": 2.0}, 0.0, 15.0, "TCA_OxPhos", "NDUFS1 and UQCRC1 and COX4I1"),
        ("ETC_FADH2", {"FADH2_m": 1.0, "ADP_m": 1.0}, {"FAD_m": 1.0, "ATP_m": 1.0}, 0.0, 15.0, "TCA_OxPhos", "SDHA and SDHB and COX4I1"),
        ("ANT_ATP", {"ATP_m": 1.0, "ADP_c": 1.0}, {"ATP_c": 1.0, "ADP_m": 1.0}, -15.0, 15.0, "Transport", "SLC25A4 or SLC25A5"),
        ("MAS_SHUTTLE", {"NADH_c": 1.0, "NAD_m": 1.0}, {"NAD_c": 1.0, "NADH_m": 1.0}, 0.0, 10.0, "Transport", "GOT1 and GOT2 and SLC25A12"),
        ("G3P_SHUTTLE", {"NADH_c": 1.0, "FAD_m": 1.0}, {"NAD_c": 1.0, "FADH2_m": 1.0}, 0.0, 10.0, "Transport", "GPD1 and GPD2"),

        # --- Decoy Branches & Potential Futile Loops ---
        ("GLYG_SYN", {"G6P_c": 1.0, "ATP_c": 1.0}, {"GLYG_c": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Glycogen", "GYS2 and UGP2"),
        ("GLYG_DEG", {"GLYG_c": 1.0}, {"G6P_c": 1.0}, 0.0, 10.0, "Glycogen", "PYGL"),
        ("GFAT_HBP", {"F6P_c": 1.0, "ATP_c": 1.0}, {"UDPGNAc_c": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Decoy_Hexosamine", "GFPT1 and UAP1"),
        ("FUTILE_ATPASE", {"ATP_c": 1.0}, {"ADP_c": 1.0}, 0.0, 10.0, "Decoy_Futile", "ENTPD5"),
        ("NNT_FORWARD", {"NADH_m": 1.0, "NADP_c": 1.0}, {"NAD_m": 1.0, "NADPH_c": 1.0}, 0.0, 10.0, "Decoy_Futile", "NNT"),
        ("NNT_REVERSE", {"NADPH_c": 1.0, "NAD_m": 1.0}, {"NADP_c": 1.0, "NADH_m": 1.0}, 0.0, 10.0, "Decoy_Futile", "NNT_REV"),

        # --- NEW: Glutamate/urea handling (fixes #2 and #3) ---
        ("GLUD", {"GLU_c": 1.0, "NAD_c": 1.0}, {"AKG_c": 1.0, "NH3_c": 1.0, "NADH_c": 1.0}, 0.0, 10.0, "Ureagenesis", "GLUD1"),
        ("UREA_SYN", {"NH3_c": 1.0, "ATP_c": 2.0}, {"UREA_c": 1.0, "ADP_c": 2.0}, 0.0, 10.0, "Ureagenesis", "CPS1 and OTC and ASS1 and ASL and ARG1"),
        ("EX_UREA_OUT", {"UREA_c": 1.0}, {}, 0.0, 10.0, "Ureagenesis", "SLC14A1"),

        # --- NEW: glycerol gluconeogenesis (adipose lipolysis byproduct,
        # a major fasting gluconeogenic substrate alongside lactate/alanine,
        # previously entirely absent) ---
        ("EX_GLY_IN", {}, {"GLY_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "AQP9"),
        ("GLY_KIN", {"GLY_c": 1.0, "ATP_c": 1.0}, {"G3P_c": 1.0, "ADP_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "GK"),
        ("G3PDH_CYT", {"G3P_c": 1.0, "NAD_c": 1.0}, {"DHAP_c": 1.0, "NADH_c": 1.0}, 0.0, 10.0, "Gluconeogenesis", "GPD1"),
        ("TPI", {"DHAP_c": 1.0}, {"GAP_c": 1.0}, -10.0, 10.0, "Gluconeogenesis", "TPI1"),

        # --- NEW: glutamine ureagenesis (second major amino-acid nitrogen
        # carrier for fasting nitrogen disposal, alongside alanine; reuses
        # the existing GLU_c -> GLUD -> NH3_c -> UREA_SYN chain) ---
        ("EX_GLN_IN", {}, {"GLN_c": 1.0}, 0.0, 10.0, "Ureagenesis", "SLC38A3 or SLC1A5"),
        ("GLS", {"GLN_c": 1.0}, {"GLU_c": 1.0, "NH3_c": 1.0}, 0.0, 10.0, "Ureagenesis", "GLS2"),

        # NEW: AKG_c sink (lumped alpha-ketoglutarate dehydrogenase ->
        # succinate thiokinase -> succinate dehydrogenase -> fumarase ->
        # malate dehydrogenase, same lumping style as TCA_OX). Needed
        # because GLUD's AKG_c byproduct was previously only consumable by
        # ALT_IN (requires alanine simultaneously present) -- isolated
        # glutamine processing had no AKG_c sink at all without this.
        # GPR lumps alpha-ketoglutarate dehydrogenase -> succinate thiokinase
        # -> succinate dehydrogenase -> fumarase -> malate dehydrogenase.
        ("AKG_TO_OAA", {"AKG_c": 1.0, "NAD_m": 1.0, "FAD_m": 1.0, "ADP_m": 1.0}, {"OAA_c": 1.0, "NADH_m": 1.0, "FADH2_m": 1.0, "ATP_m": 1.0}, 0.0, 10.0, "TCA_OxPhos", "OGDH and SUCLA2 and SDHA and FH and MDH2"),
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
        "Gluconeogenesis from lactate": {
            **base_bounds, "EX_LAC_IN": (0.0, 10.0), "EX_GLC_OUT": (1.0, 10.0),
        },
        # EX_UREA_OUT must be left OPEN (not just un-forced) here: GLUD's NH3_c
        # byproduct can only leave via UREA_SYN -> EX_UREA_OUT (no other NH3_c
        # consumer exists). Closing it (the first attempt, inherited from
        # base_bounds) forces UREA_SYN to 0, which forces GLUD to 0 (NH3_c has
        # nowhere to go), which forces ALT_IN to 0 too (GLU_c has no consumer
        # without GLUD) -- alanine couldn't enter the network at all. Real
        # biochemistry: you cannot do gluconeogenesis from an amino acid
        # without disposing of its amino group somewhere; urea is an obligate
        # co-product here, not an independent switch.
        "Gluconeogenesis from alanine": {
            **base_bounds, "EX_ALA_IN": (0.0, 10.0), "EX_GLC_OUT": (1.0, 10.0),
            "EX_UREA_OUT": (0.0, 10.0),
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
        },
        "Gluconeogenesis from glycerol": {
            **base_bounds, "EX_GLY_IN": (0.0, 10.0), "EX_GLC_OUT": (1.0, 10.0),
        },
        # Same obligate-byproduct reasoning as alanine ureagenesis (NH3_c
        # must reach EX_UREA_OUT) PLUS a second one: glutamine's carbon
        # skeleton (via GLUD's AKG_c byproduct, through the new AKG_TO_OAA
        # sink) has nowhere to go except toward gluconeogenesis -- real
        # biochemistry: disposing of glutamine's nitrogen also obligately
        # processes its carbon skeleton, typically anaplerotically/
        # gluconeogenically. EX_GLC_OUT is opened (allowed, not required).
        "Ureagenesis from glutamine": {
            **base_bounds, "EX_GLN_IN": (0.0, 10.0), "EX_UREA_OUT": (1.0, 10.0),
            "EX_GLC_OUT": (0.0, 10.0),
        },
        "Combined fasted state": {
            **base_bounds,
            "EX_FA_IN": (0.0, 10.0), "EX_LAC_IN": (0.0, 10.0), "EX_ALA_IN": (0.0, 10.0),
            "EX_GLY_IN": (0.0, 10.0), "EX_GLN_IN": (0.0, 10.0),
            "EX_GLC_OUT": (1.0, 1.0), "EX_KET_OUT": (1.0, 1.0), "EX_UREA_OUT": (1.0, 1.0),
            "ATP_MAINT": (2.0, 2.0),
        },
        # --- T0 / fed-state tasks: no new reactions needed, the original 58
        # already included forward glycolysis + glycogen synthesis + de
        # novo lipogenesis machinery (Script 19/20 just never tasked it).
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
    }


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    graph = build_expanded_core_metabolic_network()
    print(f"Expanded network: {graph.num_edges} reactions (tutorial's 58 + 11 new)")

    task_bounds = build_task_bounds()
    print("Tasks:", ", ".join(task_bounds))

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
    exchange_summary = task_flux.loc[
        ["EX_GLC_IN", "EX_FA_IN", "EX_LAC_IN", "EX_ALA_IN", "EX_GLYCOGEN_IN",
         "EX_GLY_IN", "EX_GLN_IN", "EX_GLC_OUT", "EX_KET_OUT", "EX_UREA_OUT", "ATP_MAINT"]
    ].T.round(3)
    print(exchange_summary)

    selected_df = pd.DataFrame({
        "reaction": np.asarray(reaction_ids)[selected_reactions],
    }).set_index("reaction")
    selected_df.to_csv(OUTDIR / "selected_reactions.csv")

    html = plot_task_flux_gallery(
        graph, task_flux.to_numpy(), list(task_bounds), selected_reactions,
    )
    (OUTDIR / "task_flux_gallery.html").write_text(f"<html><body>{html}</body></html>")

    summary = (
        f"Expanded network: {graph.num_edges} reactions total (tutorial's 58 + 11 new)\n"
        f"Tasks: {', '.join(task_bounds)}\n"
        f"Stage 1 (minimize union) solve status: {status_1}\n"
        f"Stage 2 (minimize flux at fixed size) solve status: {status_2}\n"
        f"Smallest joint network: {minimum_size}/{graph.num_edges} reactions\n"
    )
    (OUTDIR / "build_summary.txt").write_text(summary)
    print("\n" + summary)


if __name__ == "__main__":
    main()
