# Metabolite ID handling with OmniPath Metabo

How this use case links measured ions to MetSigDB prior knowledge through
OmniPath Metabo, written as a reference for rewiring the MetaProViz ID backend
(compatibility check, `equivalent_id()`, ID translation, MetSigDB retrieval,
ambiguity check) onto OmniPath Metabo later.

Status: interim implementation on the dev5 Postgres (build `dac1e0f66ee8`) by
direct SQL (`metabo.figures::pg_connect()`), because dev5 has no API container
and the OmniPath Metabo API and R/Python clients are still under development
(caveat 23). Code: `R/metsigdb.R` (`ion_identifiers()`,
`link_ions_to_sets()`), `preprocessing/06_metsigdb_sets.Rmd`; evidence and
alternatives: `audit/12_metabolite_identifiers.Rmd`; decisions D30–D34.

## Principle

The prior knowledge stays as it is. Only the feature side is extended: each
ion receives the database entities that are equivalent to its identifiers, and
these are matched to the unchanged MetSigDB membership table by exact entity.
This is the role `equivalent_id()` plays in MetaProViz, but the equivalence
comes from OmniPath Metabo's own entity resolution (InChIKey groups) instead
of RaMP or ID-translation tables.

## Steps

1. **Original identifiers.** Per ion, the published HMDB IDs (one per candidate
   metabolite, parallel to the names) and KEGG compound IDs (`C\d{5}`). HMDB is
   normalised to seven digits (`HMDB01659` → `HMDB0001659`). No compatibility
   cleaning: the IDs of an ion are alternative candidates (isomers, isobars),
   not synonyms, so a pairwise compatibility check would remove valid
   candidates (audit/12: 1780 candidates removed from 369 ions by the previous
   pipeline).
2. **Identifier → entities.** `identifier_evidence` (type + value) joined to
   `entity_identifier_lookup`. One identifier often maps to several entities:
   duplicate records, charge variants, name-only records (HMDB0000125
   glutathione → 3 entities; their MetSigDB memberships are spread across
   them: 196, 70 and 0 sets).
3. **Equivalent entities ("equivalent ID").** Entities are grouped with
   `chemical_resolution_group_member` at a chosen level, and the ion receives
   all entities of its groups:
   - level 1: InChIKey block 1 (connectivity; stereoisomers equivalent)
   - **level 2 (used): blocks 1+2 (connectivity + stereochemistry; charge
     variants equivalent, stereoisomers distinct)**
   - level 3: full InChIKey (identical structures only)

   Entities without an InChIKey (name-only records) have no group; they keep
   their own memberships but add no equivalents.
4. **Generic IDs are not expanded.** A KEGG ID that resolves to several
   compounds (counted by skeleton, level 1) is a class ID (e.g. C00157
   phosphatidylcholine → 790 structures). It gets no equivalents and matches
   only sets that list that KEGG compound in their `kegg` column ("rule C").
5. **MetSigDB retrieval.** `metsigdb_membership` (resources KEGG, Reactome,
   WikiPathways, MACdb, ClassyFire; one row per set and member entity), by
   exact `metabolite_entity_id`. Set selection: human pathway sets
   (`organism = 9606`); MACdb and ClassyFire have no organism; ClassyFire sets
   only with a ChemOnt ID.
6. **Ambiguity.** Counted by InChIKey skeleton (level 1), independent of the
   linking level, so that stereo or charge variants of one candidate do not
   make an ion ambiguous. An ion belongs to a set if any candidate does; each
   link records the ion's number of compounds and how many of them are in the
   set; each set reports its share of unambiguous ions (D31, D33).

## SQL (as in `R/metsigdb.R`; `%s` quoted ID lists, `%d` resolution level)

```sql
-- identifiers -> entities
SELECT CASE WHEN t.name LIKE 'Hmdb%' THEN 'HMDB' ELSE 'KEGG' END AS id_type, ie.value,
       l.entity_id::text AS entity_id
FROM omnipath.identifier_evidence ie
JOIN omnipath.vocab_identifier_type t USING (identifier_type_id)
JOIN omnipath.entity_identifier_lookup l ON l.identifier_id = ie.identifier_id
WHERE (t.name = 'Hmdb:OM:0004' AND ie.value IN (%s))
   OR (t.name = 'Kegg Compound:MI:2012' AND ie.value IN (%s));

-- entities -> InChIKey groups, and all entities of those groups (equivalents)
SELECT entity_id::text AS entity_id, group_key FROM omnipath.chemical_resolution_group_member
WHERE level_id = %d AND entity_id::text IN (%s);
SELECT group_key, entity_id::text AS entity_id FROM omnipath.chemical_resolution_group_member
WHERE level_id = %d AND group_key IN (%s);

-- MetSigDB memberships of the ion's entities
SELECT resource, set_type, set_source_id, set_label, set_size, metabolite_entity_id::text AS entity_id
FROM omnipath.metsigdb_membership
WHERE metabolite_entity_id::text IN (%s)
  AND (organism = 9606 OR (organism IS NULL AND resource IN ('ClassyFire', 'MACdb')))
  AND (resource <> 'ClassyFire' OR set_source_id LIKE 'CHEMONTID:%');

-- generic KEGG IDs: only sets that list the class compound itself
SELECT DISTINCT kegg AS value, resource, set_type, set_source_id, set_label, set_size
FROM omnipath.metsigdb_membership
WHERE <same set selection> AND kegg IN (%s);
```

## Correspondence to MetaProViz functions

| MetaProViz | Here, via OmniPath Metabo |
|---|---|
| `seed_id_compatibility_check()` | not used for linking; resolving each ID to entities shows whether IDs of one candidate agree (same InChIKey group) without RaMP |
| `equivalent_id()` | equivalent entities from `chemical_resolution_group_member` at a chosen InChIKey level |
| `translate_id()` / `traverse_ids()` | not needed: entities carry all their identifiers (`entity_identifier_lookup`) |
| `metsigdb_kegg()`, `metsigdb_reactome()`, … | `metsigdb_membership`, matched by entity instead of by one ID type |
| `count_id()` / mapping ambiguity | compounds per ion counted by InChIKey skeleton; per-set share of unambiguous ions |

## Results in this data set (audit/12, step 06)

- 4504 of 4740 original identifiers resolve to entities.
- 698 of 1096 ions are in a KEGG, Reactome or WikiPathways set (529 when
  matching KEGG IDs to KEGG pathways only).
- 519 ions have one candidate compound; 577 have several.

## Database issues found (caveat 24)

KEGG class IDs resolve to hundreds of specific lipids; C00350 also resolves to
hydroxyclomipramine; 82 WikiPathways sets without organism and label; 657 HMDB
disease annotations filed as ClassyFire chemical classes; the per-member `kegg`
column can hold KEGG DRUG IDs (ATP: D08646); the resolver lookup tables
(`resolver_chemical_identifier_lookup*`) are empty in the database, probably
because the Parquet load does not include them (the lookup Parquet files exist
in the build directory).

## Clients (state 2026-10-30)

The parquet-migration branch of omnipath-build (checkout in
`~/instances/metsig1/omnipath-build`, commit 51e7e09, 2026-10-26) has a query
layer with a Python client (`omnipath_client`) and an R client
(`omnipath.client`: `op_init`, `op_connect`, `op_sql`, `op_search`,
`op_relations`, `op_entity`). They download the published Parquet tables via
the hosted development API into a local DuckDB. The published tables are
`entity`, `entity_identifier`, `entity_annotation`, `entity_evidence`,
`relation`, `relation_annotation`, `relation_evidence` and `evidence_payloads`.
Step 2 above (identifier → entity) maps to `entity_identifier`. The InChIKey
groups (`chemical_resolution_group_member`) and the MetSigDB membership are not
among them; MetSigDB is a separate subset product
(`packages/subsets/omnipath_subsets/metsigdb`).

### Test of the R client (2026-10-30)

`omnipath.client` 0.1.0 (parquet-migration, 51e7e09) installed into the
project's private R library on beauty, with duckdb 1.5.0 from nixpkgs;
`op_init(c("hmdb", "kegg"))` from the hosted dev API, release with hmdb
2026.10.19 and kegg 2026.10.16; 0.7 min, 1.7 GB on disk (local, gitignored
`output/omnipath_client_data`).

- **Identifier → entity works.** Of the 4509 distinct original identifiers,
  4491 resolve with only these two resources installed (SQL on dev5: 4504 with
  all resources).
- **The published layer already merges entities by structure.** Chemical
  entities are keyed by full InChIKey (`entity.namespace = 'inchikey'`):
  HMDB0000125 and C00051 give one glutathione entity in both resources, where
  the dev5 database has three. 3256 of 3456 HMDB IDs resolve to one entity;
  KEGG class IDs still resolve to many (up to 411).
- **Groups:** `entity.group_connectivity` and `entity_group` (kinds
  `connectivity`, `reference`, `gene`) give InChIKey level 1. There is no
  level 2 (stereochemistry) group, so the level-2 equivalents used in step 06
  are not available from the client; level 3 (full InChIKey) is the entity
  itself.
- **MetSigDB is not published** (planned). `metsigdb_membership` in the
  database carries `inchikey` and `metabolite_structure_key` (connectivity
  block). Once MetSigDB is published, ions can be linked by InChIKey from the
  client alone.

Issues for the client developers:
1. `op_init()` fails while writing the database if DuckDB cannot autoload the
   `json` extension (duckdb from nixpkgs ships without it; "Extension
   Autoloading Error … json"). `INSTALL json` once fixes it; the client could
   install it explicitly or document it.
2. Calling `DBI::dbGetQuery()` inline through `Rscript -e` failed on a quoted
   query, while the same queries from a script file worked (quoting in the
   shell; probably not a client issue).
3. No stereochemistry-level groups in the published tables.
