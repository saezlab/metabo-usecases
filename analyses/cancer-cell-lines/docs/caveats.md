# Caveats and open points (data inventory, Phase 1b)

State of 2026-10-28. Evidence for each point is in the audit notebooks
(`audit/01`–`06`, reports in `output/reports/`). Points marked **(you)** need
input, data or a decision from you.

## Data and provenance

1. **Batch, plate and the authors' tissue labels come from an unlicensed
   source** (GitHub `shorthouse-mrc/CellLine_Metabolomics` @ cce3b15, kept
   locally). Culture batch is now a covariate in every model (audit/02), so the
   analysis depends on this file. It is not redistributed; it would be cited by
   repository and commit. **Decided (Checkpoint 1b):** use locally and cite;
   the authors may be asked later.
2. **Only normalised data are available for Shorthouse.** The published
   intensities are Lowess-corrected for injection order and plate confluence;
   ions associated with batch or plate were removed by the authors. The
   normalisation cannot be redone from these files. Less-processed data of the
   same measurements exist in the thesis folder
   (`Data/Cherkaoui_2022/CCLE2_metabolomicsData_notNormalized.xlsx`) but have
   not been inspected. **(you)** Only relevant if the normalisation should be
   questioned.
3. **Culture conditions:** all lines were grown in the same medium (RPMI,
   10 % FCS, glutamine; Cherkaoui et al. 2022) in seven batches. Medium
   differences between sources only matter if Ortmayr or other data sets are
   combined later.
4. **Local-only inputs** (not in git): `cellosaurus.txt` (release 53.0),
   DepMap 24Q4 (about 730 MB), the authors' GitHub files, Ortmayr 2019 and the
   literature full texts. The READMEs describe how to obtain them; a full run
   needs all of them.

## Samples

5. **Lung dominates:** 88 of 175 lines are lung; prostate (2) and brain (4)
   are small, which limits tissue-blocked contrasts in those tissues.
6. **Three lines without tissue** (184A1, NCIH460DNP53, SR8714) are excluded
   from tissue-adjusted models. 184A1 is a non-malignant breast epithelial
   line. **Decided:** excluded from all analyses (`analysis.exclude_lines` in
   `config.yml`).
7. **Seventeen lines were measured in several batches** (mostly repeated and
   QC lines); they are assigned the batch with most of their injections.
8. **Replicate structure:** 3 biological × 2 technical replicates for most
   lines (6 injections); the QC line MDA-MB-231 has 78. Line values are means
   over all injections, as in the original study.

## Metabolites

9. **Annotation is ambiguous:** fewer than half of the ions (519 of 1096) have
   a single candidate metabolite; flow-injection mass spectrometry cannot
   separate isomers. Results are per ion.
10. **Readout metabolites for positive controls:** glutathione-related ions are
    measured unambiguously (KEAP1/NFE2L2 usable); 2-hydroxyglutarate shares its
    ion with five candidates (IDH not specific); 5'-methylthioadenosine (MTAP)
    is not measured.
11. **Some ions carry little line-specific signal:** for 151 of 1096 ions
    (14 %), the variance between replicates of a line is larger than the
    variance between lines; the median ratio is 2.3 (audit/01, section 5).
    **Decided:** no filter; limma's variance moderation handles noisy ions, and
    the ratio is reported alongside results.
12. **Three ions with duplicate annotation are dropped** (as in the previous
    pipeline); the authors kept all 1099 ions. **(you)** Default: keep
    dropping.
13. **Values <= 0** (0.15 %) are imputed with half the metabolite minimum;
    results flag affected metabolites (audit/01).
14. **Identifier expansions with 50 or more candidates are discarded** in the
    identifier processing (step 03, inherited); affects only identifier lists.

## Prior knowledge

15. **Coverage:** 55 % of the ions have an identifier in COSMOS+, 38 % in the
    previous COSMOS network, 27 % in MetalinksDB (audit/05). **OCEAN is not
    installed**; its coverage is open. **(you)** Install OCEAN when the
    integration phase starts.

## Mutations

16. **Source decided:** DepMap 24Q4 for mutation and copy-number status,
    Cellosaurus for identity (audit/04). **Decided (D19):** oncogenes count
    with a likely-driver variant, tumour suppressors with any protein-changing
    variant or a homozygous deletion (see 21).
17. **Done:** DepMap status per line is written by `preprocessing/04_depmap_status.Rmd`.

## Other modalities

18. **No proteomics, phosphoproteomics or transcriptomics are present
    locally.** DepMap 24Q4 also contains RNA expression (e.g.
    `OmicsExpressionAllGenesTPMLogp1Profile.csv`, about 1 GB); CCLE proteomics
    (Nusinow et al. 2020) would be a separate download; no
    phosphoproteomics source has been identified. **Decided:** sources are
    chosen before the integration phase; metabolomics first.

## Other data sets

19. **Ortmayr et al. 2019** (54 NCI-60 lines, 43 shared with Shorthouse) uses a
    different value scale (relative abundances) and only deprotonated ions; it
    could serve as a replication set but would need its own processing
    (audit/06).

## Environment

20. Analyses run on beauty with MetaProViz 4.99.0 (`f76eacf`) and OmnipathR
    4.1.2 (`43056a2`) installed from GitHub into a private library; versions are
    recorded in the manifests. The Docker image planned for the repository will
    pin them.

## Contrasts

21. **Passenger variants in hypermutated lines (audit/07).** The definition of
    "altered" for tumour suppressors (any protein-changing variant, D19) also
    counts passenger variants in hypermutated lines. In DepMap, four of the
    five models with both KEAP1 and NFE2L2 altered are uterine models in the
    top 2 % of mutation burden, and in four of the five the KEAP1 variant is
    not flagged as a likely driver. This matters little in lung but can mislabel lines
    in tissue-blocked contrasts where hypermutated lines (e.g. uterus, colon)
    contribute. **Decided (2026-10-29):** keep D19; check once whether any
    final contrast is affected, and decide then.
22. **Oncogene status (audit/11).** The likely-driver rule (D19) misses
    activating variants without a DepMap hotspot flag (EGFR exon 19
    deletions in PC9 and NCIH1650; BRAF G464V) and counts truncating variants
    in oncogenes as drivers (KM12, MOLT4, NCIH1048, CORL23). **Decided
    (D29):** kept for now; revisit if a final contrast involves EGFR, BRAF,
    NRAS or PIK3CA.
23. **Access to OmniPath Metabo is interim.** The notebooks query the dev5
    Postgres tables directly with SQL (`metabo.figures::pg_connect()`;
    audit/07, audit/12, `R/metsigdb.R`, step 06), because dev5 has no API
    container. OmniPath Metabo, its API and its R and Python clients are
    under active development; once they are final, these queries are to be
    refactored to the client calls, and results re-checked against the
    released build.
