# database-content

What the OmniPath Metabo database contains: entities, interactions,
identifiers and structures per resource, overlap between resources, and
structural specificity by chemical category.

**Status (2026-11-24): fails** on dev5: `metabo_ramp_inchikey_conflict`
does not exist (nor does `metabo_entity_structural_specificity`); chemical
classes and facet bitmaps exist only in Postgres. To be moved to the Parquet release when updated
(`docs/DB_ACCESS.md`).

## Panels

The composite takes the panels listed in `composition.yaml`; the other
outputs are supplementary.

| Output (`out/`) | Shows | Data (`R/`) |
|---|---|---|
| `resource-overview-total` | Records per resource in six facets (entities, associations, interactions, identifiers, structures, literature) | `data-resource_overview.R` |
| `resource-overview-supplementary` | The same per magnitude band | `data-resource_overview.R` |
| `specificity-by-category` | Structural specificity × chemical category | `data-specificity.R` |
| `resource-overlap-networks` | Shared entities and interactions between resources | `data-resource_overlap.R` |
| `entity-interaction-matrix` | Participant type × interaction class | `data-entity_interaction.R` |
| `coverage-profile` | Items present in at least N resources | `data-coverage_profile.R` |
| `database-content.{pdf,svg}` | The composite | `composition.yaml` |

Each output comes as PDF and SVG, plus `database-content-with-caption.pdf`,
`caption.txt` and a provenance sidecar.

## Run

```sh
./rebuild.sh database-content
```
