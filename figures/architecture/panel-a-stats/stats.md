# Figure 1 Panel A — basic database statistics

Snapshot id: `fca549e16f3c`. Deployment: `dev5`. Generated: 2026-06-14T04:48:56+0200.

## Section 1 — Entities

Contributing resources: **44** (chembl, foodb, intact, stitch, bindingdb, pfocr, kegg, rhea, swisslipids, ptfi, hmdb, refmet, chebi, metatlas, hpo, recon3d, brenda, connectomedb, mirbase, reactome, signor, wikipathways, macdb, cellchat, mondo, tcdb, lipidmaps, guidetopharma, go, uniprot, drugcentral, cellinker, nichenet, corum, phenol_explorer, cellphonedb, icellnet, chemont, mrclinksdb, mebocost, psi_mi, neuronchat, omnipath_ontology, slctables)

| Metric | Value | Definition | State |
|--------|------:|------------|-------|
| `n_metabolites` | 945238 | `chemical_class:metabolite` | populated |
| `n_foods` | 152612 | `chemical_class:food` | populated |
| `n_drugs` | 1422967 | `chemical_class:drug` | populated |
| `n_lipids` | 616937 | `chemical_class:lipid` | populated |
| `n_xenobiotics` | 0 | `chemical_class:xenobiotic` | empty |
| `n_proteins_or_genes` | 460647 | `entity_type:Gene:MI:0250+Protein:MI:0326` | populated |
| `n_genes` | 178768 | `entity_type:Gene:MI:0250` | populated |
| `n_unresolved_proteins` | 281879 | `entity_type:Protein:MI:0326` | populated |

## Section 2 — Metabolite-Protein Interactions

Contributing resources: **23** (chembl, intact, stitch, bindingdb, rhea, kegg, metatlas, connectomedb, wikipathways, signor, cellchat, recon3d, guidetopharma, tcdb, drugcentral, cellinker, nichenet, cellphonedb, mrclinksdb, icellnet, mebocost, neuronchat, reactome)

| Metric | Value | Definition | State |
|--------|------:|------------|-------|
| `n_mpi_relations` | 2858507 | `mpi_base_cte:chemical+protein_or_gene` | populated |
| `n_metabolites` | 1378770 | `mpi_base_cte:distinct_chem_id` | populated |
| `n_proteins_or_genes` | 49356 | `mpi_base_cte:distinct_prot_id` | populated |
| `n_transporters` | 13560 | `transporter:tcdb+slctables` | populated |
| `n_receptors` | 3157 | `receptor:guidetopharma+connectomedb+cellphonedb` | populated |

## Section 3 — Interactions

Contributing resources: **23** (chembl, intact, stitch, bindingdb, rhea, kegg, metatlas, connectomedb, wikipathways, signor, cellchat, recon3d, guidetopharma, tcdb, drugcentral, cellinker, nichenet, cellphonedb, mrclinksdb, icellnet, mebocost, neuronchat, reactome)

| Metric | Value | Definition | State |
|--------|------:|------------|-------|
| `n_interactions` | 5060648 | `relation.category:interaction` | populated |
| `n_proteins_or_genes` | 356753 | `interactions:distinct_protein_or_gene` | populated |
| `n_pathways` | 2870 | `pathway:annotation` | populated |
| `n_reactions` | 39778 | `reaction:Reaction:OM:0015` | populated |

## Section 4 — Structures

Contributing resources: **29** (chembl, swisslipids, bindingdb, hmdb, refmet, chebi, foodb, stitch, lipidmaps, ptfi, rhea, kegg, guidetopharma, pfocr, macdb, wikipathways, metatlas, recon3d, drugcentral, reactome, signor, tcdb, intact, phenol_explorer, cellinker, mrclinksdb, mebocost, cellphonedb, neuronchat)

| Metric | Value | Definition | State |
|--------|------:|------------|-------|
| `n_structures_total` | 3137754 | `metabo_entity_structural_specificity:total` | populated |
| `n_structures_stereospecific` | 1307679 | `metabo_entity_structural_specificity.specificity_level:stereospecific` | populated |
| `n_structures_cis_trans_only` | 111549 | `metabo_entity_structural_specificity.specificity_level:cis_trans_only` | populated |
| `n_structures_constitution_only` | 874378 | `metabo_entity_structural_specificity.specificity_level:constitution_only` | populated |
| `n_structures_variable_constitution` | 11154 | `metabo_entity_structural_specificity.specificity_level:variable_constitution` | populated |
| `n_structures_unknown_constitution` | 282056 | `metabo_entity_structural_specificity.specificity_level:unknown_constitution` | populated |
| `n_structures_no_structure` | 550938 | `metabo_entity_structural_specificity.specificity_level:no_structure` | populated |
| `n_constitutional_skeletons` | 2222638 | `inchikey:first_block(standard_inchikey)` | populated |
| `n_full_inchikey` | 2573321 | `inchikey:distinct(standard_inchikey)` | populated |
| `n_ramp_conflicts` | 22014 | `metabo_ramp_inchikey_conflict:total` | populated |

## Section 5 — Annotation

Contributing resources: **44** (chembl, swisslipids, intact, hmdb, refmet, stitch, kegg, uniprot, chebi, pfocr, connectomedb, ptfi, lipidmaps, hpo, recon3d, metatlas, go, cellchat, mondo, rhea, mirbase, reactome, guidetopharma, foodb, signor, drugcentral, wikipathways, nichenet, cellinker, macdb, phenol_explorer, brenda, tcdb, corum, bindingdb, chemont, cellphonedb, icellnet, mrclinksdb, psi_mi, slctables, neuronchat, mebocost, omnipath_ontology)

| Metric | Value | Definition | State |
|--------|------:|------------|-------|
| `n_annotation_records` | 22196787 | `entity_evidence_annotation:total` | populated |
| `n_organisms` | 484 | `entity_type:Organism:OM:0032+annotated` | populated |
| `n_diseases` | 58902 | `annotation_slot3:diseases` | populated |
| `n_localizations` | 5358 | `entity_ontology_term:uberon_anatomy` | populated |
| `n_pathways_annotation` | 2870 | `entity_ontology_term:pathway_annotations` | populated |

