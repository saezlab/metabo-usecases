# MPI resource baseline snapshots

Vendored baseline snapshots of external metabolite-protein interaction (MPI)
resources for the MetaLinksDB comparison figure
(`figures/fig04-metalinks-versions/`) belong here. Loaded by
`load_vendored_mpi_snapshot()`.

Expected shape:
- `resource`
- `hmdb_id`
- `uniprot_id`
- `source`
- `relation_type`
- `source_count`
- `citation_count`
- `affinity_value`
- `curation_mode`
- `metabolite_class_label`
- `protein_class_label`
- `interaction_definition`
- `interaction_id`

Current planned resources:
- `cellphonedb.csv`
- `scconnect.csv`
- `stitch.csv`
