# Link measured ions to MetSigDB sets through OmniPath Metabo (decisions D30, D31; audit/12).
#
# Interim access: direct SQL on the dev5 Postgres tables via metabo.figures::pg_connect().
# OmniPath Metabo, its API and its R/Python clients are under active development; once they
# are final, these queries should be replaced by the client calls (see docs/caveats.md, 23).

#' Original identifiers of the ions
#'
#' HMDB IDs (one per candidate) in the seven-digit form and KEGG compound IDs,
#' as published (Dataset EV2), without the compatibility cleaning of step 03.
#'
#' @param ions Ion table with `ionIdx`, `hmdb` and `kegg` (";"-separated).
#' @return One row per ion, identifier type and value.
ion_identifiers <- function(ions) {
    split_ids <- function(x) lapply(strsplit(dplyr::coalesce(x, ""), ";"), function(v) {
        v <- trimws(v)
        v[v != ""]
    })
    digits <- function(x) sub("^HMDB", "", toupper(x))
    hmdb <- tibble::tibble(ionIdx = ions$ionIdx, value = split_ids(ions$hmdb)) |>
        tidyr::unnest(value) |>
        dplyr::filter(grepl("^[0-9]+$", digits(value))) |>
        dplyr::transmute(ionIdx, id_type = "HMDB",
                         value = paste0("HMDB", formatC(as.integer(digits(value)), width = 7, flag = "0")))
    kegg <- tibble::tibble(ionIdx = ions$ionIdx, value = split_ids(ions$kegg)) |>
        tidyr::unnest(value) |>
        dplyr::filter(grepl("^C[0-9]{5}$", value)) |>
        dplyr::transmute(ionIdx, id_type = "KEGG", value)
    dplyr::distinct(dplyr::bind_rows(hmdb, kegg))
}

#' Ion-set links from OmniPath Metabo
#'
#' The MetSigDB sets are used as they are; only the feature side is extended.
#' Each identifier of an ion is resolved to OmniPath Metabo entities, and the
#' ion receives all entities equivalent to these at chemical resolution
#' `level` (InChIKey groups; level 2 = connectivity and stereochemistry, i.e.
#' charge variants are equivalent, D32). The ion's entities are then looked up
#' in the membership table by exact entity match. KEGG IDs that resolve to
#' several structures (classes such as phosphatidylcholine) get no
#' equivalents; they match only sets that list that KEGG compound (rule C).
#' Pathway sets are restricted to human; ClassyFire and MACdb sets have no
#' organism.
#'
#' With `cross_references`, the ion also reaches records without a structure
#' (no InChIKey) that carry a specific identifier (HMDB, KEGG, ChEBI, PubChem)
#' of one of its structure entities, e.g. WikiPathways' "GSH", which carries
#' only KEGG C00051 (D35). An identifier is specific if all entities with a
#' structure that carry it share one skeleton. Such records are not merged
#' into the structure entity in this build (caveat 24).
#'
#' Different compounds are counted by InChIKey skeleton (level 1), whatever
#' `level` is: the number of skeletons per ion measures its ambiguity, and a
#' KEGG ID is generic if it resolves to several skeletons. Stereo or charge
#' variants of one compound thus do not count as separate candidates (D33).
#'
#' @param con DBI connection from `pg_connect()`.
#' @param ids Output of `ion_identifiers()`.
#' @param level Chemical resolution level (1 connectivity, 2 + stereochemistry, 3 full InChIKey).
#' @param cross_references Also link records without a structure through specific identifiers.
#' @return List with `ion_set` (one row per ion and set, with the number of
#'   the ion's skeletons in the set and the linking routes: identifier,
#'   equivalent, cross-reference, KEGG class), `ion_structures`
#'   (skeletons per ion) and `sets` (set labels and sizes).
link_ions_to_sets <- function(con, ids, level = 2, cross_references = TRUE) {
    # quoted SQL list; NULL for an empty vector, so that "IN (NULL)" matches nothing
    quote_list <- function(x) {
        x <- unique(stats::na.omit(x))
        if (!length(x)) return("NULL")
        paste(DBI::dbQuoteString(con, x), collapse = ", ")
    }

    # identifiers -> entities
    entities <- pg_query(con, sprintf("
        SELECT CASE WHEN t.name LIKE 'Hmdb%%' THEN 'HMDB' ELSE 'KEGG' END AS id_type, ie.value,
               l.entity_id::text AS entity_id
        FROM omnipath.identifier_evidence ie
        JOIN omnipath.vocab_identifier_type t USING (identifier_type_id)
        JOIN omnipath.entity_identifier_lookup l ON l.identifier_id = ie.identifier_id
        WHERE (t.name = 'Hmdb:OM:0004' AND ie.value IN (%s))
           OR (t.name = 'Kegg Compound:MI:2012' AND ie.value IN (%s))",
        quote_list(ids$value[ids$id_type == "HMDB"]), quote_list(ids$value[ids$id_type == "KEGG"])))

    # entities -> structures; entities without an InChIKey keep their own memberships but are no structure
    groups <- pg_query(con, sprintf("
        SELECT entity_id::text AS entity_id, group_key AS structure FROM omnipath.chemical_resolution_group_member
        WHERE level_id = %d AND entity_id::text IN (%s)", level, quote_list(entities$entity_id)))
    entities <- dplyr::left_join(entities, groups, by = "entity_id")
    members <- pg_query(con, sprintf("
        SELECT group_key AS structure, entity_id::text AS entity_id FROM omnipath.chemical_resolution_group_member
        WHERE level_id = %d AND group_key IN (%s)", level, quote_list(stats::na.omit(groups$structure))))

    # skeleton (level 1) of an entity, used to count different compounds (D33)
    skeleton_of <- function(entity_ids) {
        pg_query(con, sprintf("
            SELECT entity_id::text AS entity_id, group_key AS skeleton FROM omnipath.chemical_resolution_group_member
            WHERE level_id = 1 AND entity_id::text IN (%s)", quote_list(entity_ids)))
    }

    # KEGG IDs resolving to several skeletons are generic (rule C)
    generic_kegg <- entities |>
        dplyr::filter(id_type == "KEGG") |>
        dplyr::left_join(skeleton_of(entities$entity_id), by = "entity_id") |>
        dplyr::group_by(value) |>
        dplyr::summarise(structures = dplyr::n_distinct(stats::na.omit(skeleton)), .groups = "drop") |>
        dplyr::filter(structures > 1) |>
        dplyr::pull(value)

    # ion -> entities: direct, plus the equivalent entities of the ion's structures
    direct <- ids |>
        dplyr::filter(!(id_type == "KEGG" & value %in% generic_kegg)) |>
        dplyr::inner_join(entities, by = c("id_type", "value"), relationship = "many-to-many")
    ion_entities <- dplyr::bind_rows(
        dplyr::distinct(direct, ionIdx, entity_id, structure) |> dplyr::mutate(via = "identifier"),
        direct |>
            dplyr::filter(!is.na(structure)) |>
            dplyr::distinct(ionIdx, structure) |>
            dplyr::inner_join(members, by = "structure", relationship = "many-to-many") |>
            dplyr::mutate(via = "equivalent")
    ) |>
        dplyr::distinct(ionIdx, entity_id, .keep_all = TRUE)
    ion_entities <- dplyr::left_join(ion_entities, skeleton_of(ion_entities$entity_id), by = "entity_id")

    # cross-references (D35): specific identifiers of the ion's structure entities lead to
    # records without a structure that carry them
    if (cross_references) {
        source_skeleton <- ion_entities |>
            dplyr::filter(!is.na(skeleton)) |>
            dplyr::distinct(source_entity = entity_id, source_skeleton = skeleton)
        carriers <- pg_query(con, sprintf("
            WITH src AS (
                SELECT l.entity_id, l.identifier_id
                FROM omnipath.entity_identifier_lookup l
                JOIN omnipath.identifier_evidence ie ON ie.identifier_id = l.identifier_id
                JOIN omnipath.vocab_identifier_type t ON t.identifier_type_id = ie.identifier_type_id
                WHERE l.entity_id::text IN (%s)
                  AND t.name IN ('Hmdb:OM:0004', 'Kegg Compound:MI:2012', 'Chebi:MI:0474', 'Pubchem Compound:OM:0002'))
            SELECT src.entity_id::text AS source_entity, src.identifier_id::text AS identifier_id,
                   l.entity_id::text AS entity_id
            FROM src JOIN omnipath.entity_identifier_lookup l ON l.identifier_id = src.identifier_id
            WHERE l.entity_id <> src.entity_id", quote_list(source_skeleton$source_entity)))
        carriers <- carriers |>
            dplyr::select(source_entity, identifier_id, entity_id) |>
            dplyr::left_join(skeleton_of(carriers$entity_id), by = "entity_id") |>
            dplyr::left_join(source_skeleton, by = "source_entity", relationship = "many-to-many")
        # an identifier is specific if all structures carrying it share one skeleton
        specific <- carriers |>
            dplyr::group_by(identifier_id) |>
            dplyr::summarise(skeletons = dplyr::n_distinct(c(source_skeleton, stats::na.omit(skeleton))),
                             .groups = "drop") |>
            dplyr::filter(skeletons == 1) |>
            dplyr::pull(identifier_id)
        cross <- carriers |>
            dplyr::filter(identifier_id %in% specific, is.na(skeleton)) |>
            dplyr::distinct(source_entity, entity_id, source_skeleton) |>
            dplyr::inner_join(dplyr::distinct(ion_entities, ionIdx, source_entity = entity_id),
                              by = "source_entity", relationship = "many-to-many") |>
            dplyr::transmute(ionIdx, entity_id, structure = NA_character_, via = "cross-reference",
                             skeleton = source_skeleton) |>
            dplyr::anti_join(ion_entities, by = c("ionIdx", "entity_id"))
        ion_entities <- dplyr::bind_rows(ion_entities, dplyr::distinct(cross, ionIdx, entity_id, .keep_all = TRUE))
    }

    # human pathway sets; ClassyFire and MACdb have no organism (82 WikiPathways sets without
    # organism and label are left out); ClassyFire sets need a ChemOnt ID (657 HMDB disease
    # annotations are filed under ClassyFire in this build, caveat 24)
    organism_filter <- paste(
        "(organism = 9606 OR (organism IS NULL AND resource IN ('ClassyFire', 'MACdb')))",
        "AND (resource <> 'ClassyFire' OR set_source_id LIKE 'CHEMONTID:%')"
    )
    membership <- pg_query(con, sprintf("
        SELECT resource, set_type, set_source_id, set_label, set_size, metabolite_entity_id::text AS entity_id
        FROM omnipath.metsigdb_membership
        WHERE metabolite_entity_id::text IN (%s) AND %s",
        quote_list(ion_entities$entity_id), organism_filter))
    kegg_listed <- pg_query(con, sprintf("
        SELECT DISTINCT kegg AS value, resource, set_type, set_source_id, set_label, set_size
        FROM omnipath.metsigdb_membership
        WHERE %s AND kegg IN (%s)", organism_filter, quote_list(generic_kegg)))

    ion_structures <- ion_entities |>
        dplyr::filter(!is.na(skeleton)) |>
        dplyr::distinct(ionIdx, skeleton) |>
        dplyr::count(ionIdx, name = "ion_structures")

    by_structure <- ion_entities |>
        dplyr::inner_join(membership, by = "entity_id", relationship = "many-to-many") |>
        dplyr::group_by(ionIdx, resource, set_source_id) |>
        dplyr::summarise(structures_in_set = dplyr::n_distinct(stats::na.omit(skeleton)),
                         route = paste(sort(unique(via)), collapse = "+"), .groups = "drop")
    by_kegg_class <- ids |>
        dplyr::filter(id_type == "KEGG", value %in% generic_kegg) |>
        dplyr::inner_join(kegg_listed, by = "value", relationship = "many-to-many") |>
        dplyr::distinct(ionIdx, resource, set_source_id) |>
        dplyr::mutate(structures_in_set = NA_integer_, route = "KEGG class")

    ion_set <- dplyr::bind_rows(by_structure, by_kegg_class) |>
        dplyr::group_by(ionIdx, resource, set_source_id) |>
        dplyr::summarise(structures_in_set = max(structures_in_set, -1L, na.rm = TRUE),
                         route = paste(sort(unique(unlist(strsplit(route, "+", fixed = TRUE)))), collapse = "+"),
                         .groups = "drop") |>
        dplyr::mutate(structures_in_set = dplyr::na_if(structures_in_set, -1L)) |>
        dplyr::left_join(ion_structures, by = "ionIdx")

    sets <- dplyr::bind_rows(
        dplyr::distinct(membership, resource, set_type, set_source_id, set_label, set_size),
        dplyr::distinct(kegg_listed, resource, set_type, set_source_id, set_label, set_size)
    ) |>
        dplyr::distinct(resource, set_source_id, .keep_all = TRUE)

    list(ion_set = ion_set, ion_structures = ion_structures, sets = sets, generic_kegg = generic_kegg)
}
