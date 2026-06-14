\pset pager off
\pset format unaligned
\pset tuples_only on

WITH
table_inventory AS (
  SELECT jsonb_agg(
    jsonb_build_object(
      'schema', schema_name,
      'relation', relation_name,
      'kind', relation_kind,
      'estimated_rows', estimated_rows,
      'total_size_bytes', total_size_bytes,
      'total_size_pretty', pg_size_pretty(total_size_bytes)
    )
    ORDER BY total_size_bytes DESC, relation_name
  ) AS data
  FROM (
    SELECT
      n.nspname AS schema_name,
      c.relname AS relation_name,
      CASE c.relkind
        WHEN 'r' THEN 'table'
        WHEN 'p' THEN 'partitioned table'
        WHEN 'm' THEN 'materialized view'
        WHEN 'v' THEN 'view'
        ELSE c.relkind::text
      END AS relation_kind,
      COALESCE(s.n_live_tup, c.reltuples)::bigint AS estimated_rows,
      pg_total_relation_size(c.oid)::bigint AS total_size_bytes
    FROM pg_class c
    JOIN pg_namespace n ON n.oid = c.relnamespace
    LEFT JOIN pg_stat_user_tables s ON s.relid = c.oid
    WHERE c.relkind IN ('r', 'p', 'm', 'v')
      AND n.nspname NOT IN ('pg_catalog', 'information_schema')
      AND n.nspname !~ '^pg_toast'
  ) rows
),
resource_contributions AS (
  SELECT jsonb_agg(
    jsonb_build_object(
      'resource_id', resource_id,
      'resource_name', COALESCE(resource_name, resource_id),
      'resource_kind', resource_kind,
      'build_status', build_status,
      'categories', COALESCE(categories, '[]'::jsonb),
      'annotation_ontologies', COALESCE(annotation_ontologies, '[]'::jsonb),
      'entity_count', entity_count,
      'interaction_count', interaction_count,
      'association_count', association_count,
      'identifier_count', identifier_count,
      'ontology_term_count', ontology_term_count,
      'total_size_bytes', total_size_bytes,
      'last_downloaded_at', last_downloaded_at,
      'last_built_at', last_built_at
    )
    ORDER BY total_size_bytes DESC, resource_id
  ) AS data
  FROM resources
),
resource_category_counts AS (
  SELECT COALESCE(
    jsonb_agg(
      jsonb_build_object(
        'category', category,
        'resource_count', resource_count,
        'entity_count', entity_count,
        'interaction_count', interaction_count,
        'association_count', association_count
      )
      ORDER BY resource_count DESC, category
    ),
    '[]'::jsonb
  ) AS data
  FROM (
    SELECT
      category,
      count(*)::bigint AS resource_count,
      sum(entity_count)::bigint AS entity_count,
      sum(interaction_count)::bigint AS interaction_count,
      sum(association_count)::bigint AS association_count
    FROM resources
    CROSS JOIN LATERAL jsonb_array_elements_text(COALESCE(categories, '[]'::jsonb)) AS category
    GROUP BY category
  ) rows
),
resource_annotation_ontologies AS (
  SELECT COALESCE(
    jsonb_agg(
      jsonb_build_object(
        'ontology', ontology,
        'resource_count', resource_count,
        'ontology_term_count', ontology_term_count
      )
      ORDER BY resource_count DESC, ontology
    ),
    '[]'::jsonb
  ) AS data
  FROM (
    SELECT
      ontology,
      count(*)::bigint AS resource_count,
      sum(ontology_term_count)::bigint AS ontology_term_count
    FROM resources
    CROSS JOIN LATERAL jsonb_array_elements_text(COALESCE(annotation_ontologies, '[]'::jsonb)) AS ontology
    GROUP BY ontology
  ) rows
),
entity_facets AS (
  SELECT COALESCE(jsonb_object_agg(facet_name, facet_values), '{}'::jsonb) AS data
  FROM (
    SELECT
      facet_name,
      jsonb_agg(
        jsonb_build_object(
          'label', facet_value,
          'count', entity_count
        )
        ORDER BY entity_count DESC, facet_value
      ) AS facet_values
    FROM facet_entity_bitmap
    WHERE facet_name IN ('entity_type', 'source', 'taxonomy_id', 'ontology_id')
    GROUP BY facet_name
  ) rows
),
relation_facets AS (
  SELECT COALESCE(jsonb_object_agg(facet_name, facet_values), '{}'::jsonb) AS data
  FROM (
    SELECT
      facet_name,
      jsonb_agg(
        jsonb_build_object(
          'label', facet_value,
          'category', facet_category,
          'count', relation_count
        )
        ORDER BY relation_count DESC, facet_value
      ) AS facet_values
    FROM facet_relation_bitmap
    WHERE facet_name IN ('predicate', 'source', 'participant_type')
    GROUP BY facet_name
  ) rows
),
subject_predicate_object_types AS (
  SELECT COALESCE(
    jsonb_agg(
      jsonb_build_object(
        'subject_type', subject_type,
        'predicate', predicate,
        'object_type', object_type,
        'category', relation_category,
        'relation_count', relation_count,
        'evidence_count', evidence_count
      )
      ORDER BY relation_count DESC, subject_type, predicate, object_type
    ),
    '[]'::jsonb
  ) AS data
  FROM (
    SELECT
      COALESCE(subject_entity.entity_type, 'Unspecified') AS subject_type,
      er.predicate,
      COALESCE(object_entity.entity_type, 'Unspecified') AS object_type,
      er.relation_category,
      count(*)::bigint AS relation_count,
      sum(er.evidence_count)::bigint AS evidence_count
    FROM entity_relation er
    JOIN entity subject_entity ON subject_entity.entity_id = er.subject_entity_id
    JOIN entity object_entity ON object_entity.entity_id = er.object_entity_id
    GROUP BY
      COALESCE(subject_entity.entity_type, 'Unspecified'),
      er.predicate,
      COALESCE(object_entity.entity_type, 'Unspecified'),
      er.relation_category
    ORDER BY relation_count DESC, subject_type, er.predicate, object_type
    LIMIT 90
  ) rows
),
ontology_prefix_counts AS (
  SELECT COALESCE(
    jsonb_agg(
      jsonb_build_object(
        'ontology_prefix', ontology_prefix,
        'term_count', term_count
      )
      ORDER BY term_count DESC, ontology_prefix
    ),
    '[]'::jsonb
  ) AS data
  FROM (
    SELECT
      COALESCE(ontology_prefix, 'unknown') AS ontology_prefix,
      count(*)::bigint AS term_count
    FROM ontology_terms
    GROUP BY COALESCE(ontology_prefix, 'unknown')
  ) rows
),
annotation_terms AS (
  SELECT jsonb_build_object(
    'entity',
    (
      SELECT COALESCE(
        jsonb_agg(
          jsonb_build_object(
            'term_entity_id', term_entity_id,
            'term_id', term_id,
            'label', label,
            'ontology_prefix', ontology_prefix,
            'global_count', global_count
          )
          ORDER BY global_count DESC, term_id
        ),
        '[]'::jsonb
      )
      FROM (
        SELECT
          b.term_entity_id,
          COALESCE(ot.term_id, b.term_entity_id::text) AS term_id,
          COALESCE(ot.label, ot.term_id, b.term_entity_id::text) AS label,
          COALESCE(ot.ontology_prefix, 'unknown') AS ontology_prefix,
          b.global_count
        FROM annotation_term_entity_bitmap b
        LEFT JOIN ontology_terms ot USING (term_entity_id)
        ORDER BY b.global_count DESC, COALESCE(ot.term_id, b.term_entity_id::text)
        LIMIT 50
      ) rows
    ),
    'relation',
    (
      SELECT COALESCE(
        jsonb_agg(
          jsonb_build_object(
            'term_entity_id', term_entity_id,
            'term_id', term_id,
            'label', label,
            'ontology_prefix', ontology_prefix,
            'global_count', global_count
          )
          ORDER BY global_count DESC, term_id
        ),
        '[]'::jsonb
      )
      FROM (
        SELECT
          b.term_entity_id,
          COALESCE(ot.term_id, b.term_entity_id::text) AS term_id,
          COALESCE(ot.label, ot.term_id, b.term_entity_id::text) AS label,
          COALESCE(ot.ontology_prefix, 'unknown') AS ontology_prefix,
          b.global_count
        FROM annotation_term_relation_bitmap b
        LEFT JOIN ontology_terms ot USING (term_entity_id)
        ORDER BY b.global_count DESC, COALESCE(ot.term_id, b.term_entity_id::text)
        LIMIT 50
      ) rows
    )
  ) AS data
),
top_connected_entities AS (
  SELECT COALESCE(
    jsonb_agg(
      jsonb_build_object(
        'entity_id', entity_id,
        'entity_key', entity_key,
        'canonical_identifier', canonical_identifier,
        'canonical_identifier_type', canonical_identifier_type,
        'label', display_label,
        'entity_type', entity_type,
        'taxonomy_id', taxonomy_id,
        'relation_count', relation_count
      )
      ORDER BY relation_count DESC, entity_id
    ),
    '[]'::jsonb
  ) AS data
  FROM (
    SELECT
      h.entity_id,
      h.entity_key,
      h.canonical_identifier,
      h.canonical_identifier_type,
      COALESCE(NULLIF(left(label_attr.value, 90), ''), h.canonical_identifier) AS display_label,
      h.entity_type,
      h.taxonomy_id,
      h.relation_count
    FROM (
      SELECT
        e.entity_id,
        e.entity_key,
        e.canonical_identifier,
        e.canonical_identifier_type,
        e.entity_attributes,
        e.entity_type,
        e.taxonomy_id,
        c.relation_count
      FROM entity_relation_counts c
      JOIN entity e USING (entity_id)
      WHERE e.entity_type IS DISTINCT FROM 'OM:0012:Cv Term'
      ORDER BY c.relation_count DESC, e.entity_id
      LIMIT 50
    ) h
    LEFT JOIN LATERAL (
      SELECT attr->>'value' AS value
      FROM jsonb_array_elements(COALESCE(h.entity_attributes, '[]'::jsonb)) attr
      WHERE attr->>'term' IN ('OM:0202:Name', 'OM:0613:Description', 'OM:0203:Synonym')
      ORDER BY CASE attr->>'term'
        WHEN 'OM:0202:Name' THEN 1
        WHEN 'OM:0613:Description' THEN 2
        WHEN 'OM:0203:Synonym' THEN 3
        ELSE 4
      END
      LIMIT 1
    ) label_attr ON true
  ) rows
),
totals AS (
  SELECT jsonb_build_object(
    'entity_count', COALESCE((SELECT sum(entity_count)::bigint FROM facet_entity_bitmap WHERE facet_name = 'entity_type'), 0),
    'relation_count', COALESCE((SELECT sum(relation_count)::bigint FROM facet_relation_bitmap WHERE facet_name = 'predicate'), 0),
    'interaction_count', COALESCE((SELECT sum(relation_count)::bigint FROM facet_relation_bitmap WHERE facet_name = 'predicate' AND facet_category = 'interaction'), 0),
    'association_count', COALESCE((SELECT sum(relation_count)::bigint FROM facet_relation_bitmap WHERE facet_name = 'predicate' AND facet_category = 'association'), 0),
    'resource_count', (SELECT count(*)::bigint FROM resources),
    'built_resource_count', (SELECT count(*)::bigint FROM resources WHERE build_status = 'success'),
    'ontology_term_count', (SELECT count(*)::bigint FROM ontology_terms),
    'database_size_bytes', (
      SELECT COALESCE(sum((item->>'total_size_bytes')::bigint), 0)
      FROM table_inventory, jsonb_array_elements(table_inventory.data) item
    )
  ) AS data
)
SELECT jsonb_pretty(
  jsonb_build_object(
    'generated_at', to_char(now() AT TIME ZONE 'UTC', 'YYYY-MM-DD"T"HH24:MI:SS"Z"'),
    'totals', totals.data,
    'table_inventory', table_inventory.data,
    'resources', resource_contributions.data,
    'resource_categories', resource_category_counts.data,
    'resource_annotation_ontologies', resource_annotation_ontologies.data,
    'entity_facets', entity_facets.data,
    'relation_facets', relation_facets.data,
    'subject_predicate_object_types', subject_predicate_object_types.data,
    'ontology_prefix_counts', ontology_prefix_counts.data,
    'annotation_terms', annotation_terms.data,
    'top_connected_entities', top_connected_entities.data
  )
)
FROM totals
CROSS JOIN table_inventory
CROSS JOIN resource_contributions
CROSS JOIN resource_category_counts
CROSS JOIN resource_annotation_ontologies
CROSS JOIN entity_facets
CROSS JOIN relation_facets
CROSS JOIN subject_predicate_object_types
CROSS JOIN ontology_prefix_counts
CROSS JOIN annotation_terms
CROSS JOIN top_connected_entities;
