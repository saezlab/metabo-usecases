# Figure 1 — OmniPath Metabo architecture and database content

**Mode**: mixed (pipeline panels + manual-asset logos vendored in Panel A)

**Caption draft**

> Architecture and content of the OmniPath Metabo database. **(A)**
> Workflow from upstream resources through download / processing
> (`omnipath-resources`), parallel database builds (`omnipath-utils`,
> `omnipath-build`, `omnipath-metabo`), the three web services
> (`utils.omnipathdb.org`, `omnipathdb.org`, `metabo.omnipathdb.org`)
> with the OmniPath Metabo web app, and the clients (OmnipathR for R,
> omnipath-client for Python, `annnet` and Cytoscape for network
> exploration). Inside the metabo column, MetalinksDB, the COSMOS PKN,
> and the server-side cheminformatics layer (RDKit) appear as
> first-class datasets. MetalinksDB is built as a materialized view
> inside the main OmniPath Postgres database during the main build and
> exposed by the Metabo web service. **(B)** Number of molecular
> entities by upstream resource. **(C)** Number of interactions by
> upstream resource. **(D)** Number of interactions by interaction
> type. **(E)** Number of annotation classes by upstream resource.
> **(F)** Number of ontology terms by ontology. Counts come from a
> single pinned OmniPath build, recorded in the figure's provenance
> sidecar.

**Files**

- `build.R` — orchestration. Renders Panels B–F via the package's
  `db_content` plot module, compiles Panel A from TikZ via xelatex,
  composes via the `tex/compose_fig01.tex` template, and emits the
  provenance sidecar.
- `out/` — emitted artifacts: `fig01-overview.{pdf,svg}` (composite),
  `panel{A,B,C,D,E,F}.pdf` (constituents), `fig01-overview.pdf.provenance.json`.
- `manual/` — placeholder for any manual sub-assets specific to this
  figure (currently empty — all logos for Panel A live in the shared
  `inst/extdata/assets/` per FR-005b).
