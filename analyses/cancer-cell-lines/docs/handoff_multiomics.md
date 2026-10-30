# Handoff: multi-omics expansion of the cancer cell lines use case

For the next, exploratory branch (to be created from `cancer-cell-lines-v2`
at or after `ef92bad`). Read `use_case_overview.md` first; this document says
where things are, what to carry over, what is missing, and how to start.

## Goal of the new branch

Use more of the data and more of the prior knowledge around the contrasts
that the metabolomics phase established:

- **NRF2 activation in lung lines** (KEAP1/NFE2L2; positive control)
- **Neuroendocrine vs non-small cell lung lines**
- **Colon, kidney and skin vs the other lines**

The branch is exploratory: findings feed back into the choice of final
contrasts and the manuscript, but nothing on it is confirmatory unless
pre-specified as before (contrast chosen before results, positive control
fixed before the fit).

## Starting point

- **Code:** branch `cancer-cell-lines-v2`.
  - Entry point `analyses/cancer-cell-lines/run_analysis.R` (`--list` shows the steps).
  - Helpers in `analyses/cancer-cell-lines/R/`:
    - `mutations.R`: alteration status;
    - `differential.R`: contrast lines, design, limma fit, exploratory and tissue-vs-rest lines;
    - `metsigdb.R`: ion–set links;
    - `plotting.R`, `transform.R`, `cellosaurus.R`, `use_case.R`.
  - Contrast definitions in `config/config.yml` (`differential`, `exploratory`, `enrichment`).
- **Analysis data set:** `data/derived/cancer-cell-lines/shorthouse_v1/shorthouse_analysis_v1.rds`, written by `preprocessing/05`. It contains:
  - log2 values and raw values;
  - the line table, with tissue, batch, histology, DepMap model and modelling flags;
  - the ion table;
  - mutation status, deletions and copy number.

  Ion–set links come from `06_ion_sets.tsv` and `06_sets.tsv`.
- **Full result tables** (all ions, signed statistics) for every contrast are
  in `output/05_*`, `output/06_*`, `output/07_*/exploratory_results.tsv` and
  `output/08_enrichment/`. These are the inputs for footprint methods.

## Environment (beauty)

- **Clone:** `~/dev/metabo-usecases-restructure`. Edit locally, `rsync -a
  --chmod=D755,F644` to beauty, render there. After each push, run `git fetch &&
  git reset --mixed origin/<branch>` on beauty so the manifests record the commit.
- **Render a step:**

  ```bash
  R_LIBS=$PWD/.Rlib nix-shell -p pandoc --run \
    "~/.local/bin/with-metabo-figures-env Rscript analyses/cancer-cell-lines/run_analysis.R <steps>"
  ```

- **Private R library** `.Rlib`:
  - MetaProViz 4.99.0, OmnipathR 4.1.2, `omnipath.client` 0.1.0;
  - decoupleR 2.16.0 and limma come from the Nix environment;
  - duckdb 1.5.0 is in the Nix store: add `/nix/store/j2m4xz2fxihi8vkyzmq0y8hspn2w62kj-r-duckdb-1.5.0/library` to `R_LIBS` when using the client.
  - Never install into the shared library `~/R/x86_64-pc-linux-gnu-library/4.5`.
- **Side effect to know:** loading the fresh OmnipathR clears the OmnipathR
  cache of the `omnipath` user.
- **OmniPath Metabo:**
  - dev5 Postgres via `metabo.figures::pg_connect()`, read only; there is no dev5 API.
  - The R/Python clients are on the parquet-migration branch of omnipath-build (`~/instances/metsig1/omnipath-build/packages/query`).
  - MetSigDB is not yet published as Parquet; until then use SQL (caveat 23).
  - Instances on beauty use `~/instances/lockctl`: claim before touching an instance.
- **Do not touch** the other beauty checkouts (`~/dev/metabo-usecases`,
  `-tables`, `-fig3`).

## Data

Present (local only where marked, never committed):

- **Shorthouse** EV1–EV9 (`data/raw/cancer-cell-lines/`).
- The **authors' GitHub files** (batch mapping, CCLE calls; local only, no licence).
- **DepMap 24Q4** `Model.csv`, somatic mutations, absolute copy number (`depmap_24Q4/`, local only, CC BY 4.0).
- **Cellosaurus 53.0.**
- **Ortmayr 2019** (local only; descriptive audit only).

Needed for multi-omics (each a download, to be approved):

- **Transcriptomics:** DepMap 24Q4 RNA expression, e.g.
  `OmicsExpressionProteinCodingGenesTPMLogp1.csv` or the all-genes profile file
  (about 1 GB, Figshare, same DOI as the other DepMap files, CC BY 4.0).
  Matching via `ModelID` is already in the line table.
- **Proteomics:** CCLE proteomics (Nusinow et al. 2020); source and licence to
  be checked.
- **Phosphoproteomics:** no source identified yet.
- First check for each: **overlap with the 167 lines** and with the lines of
  each contrast. Contrast groups shrink in every modality.

## Ideas and how to approach them

Suggested order.

### 1. Signed enzyme footprints from reactions

- **What:** OmniPath Metabo holds reactions with substrate (`has_input`),
  product (`has_output`) and enzyme (`enabled_by`, `catalyzes`) relations from
  Rhea, KEGG, Reactome, Recon3D and Metatlas/Human-GEM (`relation_evidence`,
  predicates and roles in `vocab_relation_predicate` / `vocab_relation_role`).
- **Footprint per enzyme:** products +1, substrates −1. Score enzymes on the
  ions' moderated t-statistics, e.g. decoupleR ULM with signed `mor`; ions
  linked as in step 06.
- **Why:** substrates and products move in opposite directions, so unsigned
  pathway sets cancel coherent changes (audit/13, analysis/08). This is the
  logic of OCEAN.
- **Decide first:**
  - which reaction resource(s);
  - cofactor exclusion (MetSigDB cofactor lists from the metsig1 instance: "cofactors" and "cofactors_extended");
  - how ion ambiguity enters;
  - which significance statement is possible, given the correlation between ions.
- **Pre-specified positive control:** glutathione synthesis enzymes (GCLC,
  GCLM, GSS) active in NRF2-altered lung lines.

### 2. Product/substrate ratios

- Log ratios of selected product–substrate pairs (e.g. glutathione/glutamate,
  γ-glutamylcysteine/glutamate), tested with the same limma models.
- Ratios cancel line-wide effects.
- Pre-specified reactions only, no genome-wide screen.

### 3. The strong biological axis

- After batch and tissue, one component explains 20 % of the residual
  variance (phosphatidylglycerols, citrate, nucleotides; audit/13).
- Relate it to doubling time (NCI-60; used by Shorthouse et al.), DepMap
  expression signatures, and mutation burden.
- If it reflects proliferation or mitochondrial content, that is a result in
  itself and matters for every contrast.

### 4. Transcriptomics and proteomics for the same contrasts

- Same lines (restricted to the overlap), same designs, limma or voom as fits
  the data.
- **NRF2 contrast:** do NRF2 target genes and glutathione enzymes rise in the
  lines with high glutathione?
- **Enzyme level:** compare enzyme activities from idea 1 with enzyme
  expression and protein.
- **TF activities:** the plan names decoupleR with CollecTRI. In audit/07 you
  preferred not to add CollecTRI as a separate tool, so ask before using it.
  An OmniPath Metabo source may be possible.

### 5. Later

COSMOS-type integration with signed upstream activities and downstream
metabolites (optionally MetalinksDB); metabolite–dependency associations
(DepMap CRISPR).

## What to carry over

- **Contrast definitions and design rules:**
  - culture batch in every metabolomics model;
  - tissue blocking in pan-tissue contrasts;
  - histology blocking in lung;
  - co-mutations as covariates, not exclusion;
  - copy number required for MTAP/CDKN2A.

  Batch and plate structure in other modalities must be checked anew.
- **Status definitions** (D19, D24, D25, D27; known limitations in caveats 21–22).
- **Ion–set linking** (D30–D35): feature side only, prior knowledge unchanged.
- **Imputation guard** for small groups (count of imputed values per ion and group).
- **Working rules** (planning prompt and your feedback):
  - ask before downloads, installs and commits;
  - every decision-relevant analysis gets a notebook with plots and a written interpretation;
  - simplest standard methods;
  - comments above code;
  - scrollable previews of data;
  - groupings of genes need a verified reference and prior-knowledge or data support;
  - prefer the OmniPath Metabo APIs and clients where they exist;
  - no unlicensed data or full texts in commits;
  - no personal absolute paths;
  - stage files explicitly.

## First steps on the new branch

1. Create the branch from `cancer-cell-lines-v2`. Decide whether the
   multi-omics work lives in `analyses/cancer-cell-lines/` (new numbered steps)
   or in a new `analyses/cancer-cell-lines-multiomics/` that reads the
   metabolomics outputs.
2. Approve and download DepMap expression; build the line overlap table per
   contrast (notebook).
3. Write the analysis plan for the branch, with decision points and
   pre-specified positive controls, before fitting anything.
4. Start with idea 1 (footprints) on the NRF2 contrast, since it has a clear
   positive control.

## Open questions for you

- Which modalities first: transcriptomics only, or also proteomics?
- Folder layout: same use-case folder or a new one?
- Reaction resource for footprints: Rhea, Recon3D/Metatlas (genome-scale
  models), or several?
- Is a TF-activity analysis wanted, and with which prior knowledge?
- How far can the branch go before results are shown to the team?
