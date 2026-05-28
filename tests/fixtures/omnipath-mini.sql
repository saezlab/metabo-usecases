-- Minimal OmniPath schema + synthetic data for CI / R CMD check.
-- Mirrors a thin slice of the live schema on dev3 (see tests/fixtures/README.md).

DROP TABLE IF EXISTS entity_evidence CASCADE;
DROP TABLE IF EXISTS entity_annotation_relation CASCADE;
DROP TABLE IF EXISTS relation CASCADE;
DROP TABLE IF EXISTS annotation CASCADE;
DROP TABLE IF EXISTS entity CASCADE;
DROP TABLE IF EXISTS ontology_terms CASCADE;
DROP TABLE IF EXISTS dataset CASCADE;
DROP TABLE IF EXISTS data_source CASCADE;

CREATE TABLE data_source (
    source_id BIGSERIAL PRIMARY KEY,
    name      TEXT NOT NULL UNIQUE
);

CREATE TABLE dataset (
    dataset_id BIGSERIAL PRIMARY KEY,
    source_id  BIGINT NOT NULL REFERENCES data_source(source_id),
    name       TEXT NOT NULL,
    UNIQUE (source_id, name)
);

CREATE TABLE entity (
    entity_id BIGSERIAL PRIMARY KEY,
    label     TEXT NOT NULL
);

CREATE TABLE entity_evidence (
    entity_evidence_id BIGSERIAL PRIMARY KEY,
    entity_id BIGINT NOT NULL REFERENCES entity(entity_id),
    source_id BIGINT NOT NULL REFERENCES data_source(source_id)
);

CREATE TABLE relation (
    relation_id BIGSERIAL PRIMARY KEY,
    source_id   BIGINT NOT NULL REFERENCES data_source(source_id),
    predicate   TEXT NOT NULL,
    subject_id  BIGINT NOT NULL REFERENCES entity(entity_id),
    object_id   BIGINT NOT NULL REFERENCES entity(entity_id)
);

CREATE TABLE annotation (
    annotation_id BIGSERIAL PRIMARY KEY,
    label         TEXT NOT NULL,
    ontology      TEXT
);

CREATE TABLE entity_annotation_relation (
    entity_annotation_id BIGSERIAL PRIMARY KEY,
    entity_id     BIGINT NOT NULL REFERENCES entity(entity_id),
    annotation_id BIGINT NOT NULL REFERENCES annotation(annotation_id),
    source_id     BIGINT NOT NULL REFERENCES data_source(source_id)
);

CREATE TABLE ontology_terms (
    term_id BIGSERIAL PRIMARY KEY,
    ontology TEXT NOT NULL,
    label    TEXT NOT NULL
);

INSERT INTO data_source(name) VALUES
    ('signor'), ('chebi'), ('omnipath_metabo'), ('wikipathways');

INSERT INTO entity(label) VALUES
    ('TP53'), ('KRAS'), ('EGFR'), ('azelaic_acid'), ('citrate'), ('MAPK1');

INSERT INTO entity_evidence(entity_id, source_id) VALUES
    (1, 1), (2, 1), (3, 1),
    (4, 2), (5, 2),
    (6, 3),
    (1, 4), (3, 4);

INSERT INTO relation(source_id, predicate, subject_id, object_id) VALUES
    (1, 'phosphorylates', 1, 2),
    (1, 'binds',          2, 3),
    (1, 'inhibits',       3, 6),
    (3, 'regulates',      4, 5),
    (3, 'allosteric',     5, 6),
    (4, 'participates_in', 1, 6);

INSERT INTO annotation(label, ontology) VALUES
    ('membrane',    'GO_CC'),
    ('cytoplasm',   'GO_CC'),
    ('apoptosis',   'GO_BP'),
    ('metabolism',  'GO_BP');

INSERT INTO entity_annotation_relation(entity_id, annotation_id, source_id) VALUES
    (1, 3, 1), (2, 1, 1), (3, 1, 1),
    (4, 4, 2), (5, 4, 2),
    (6, 1, 3);

INSERT INTO ontology_terms(ontology, label) VALUES
    ('GO_CC', 'membrane'), ('GO_CC', 'cytoplasm'),
    ('GO_BP', 'apoptosis'), ('GO_BP', 'metabolism'),
    ('CHEBI', 'azelaic acid');
