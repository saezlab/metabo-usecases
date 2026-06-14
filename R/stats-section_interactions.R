#' Section 3 — Interactions (FR-043c)
#'
#' Reports the four headline interaction counts. The section
#' EXPLICITLY includes the metabolite-protein interactions reported
#' in Section 2 (per FR-043c paragraph 1 — keeps the relative scale
#' visible on the Panel A diagram). Pathway and reaction metrics
#' carry a preferred-vs-fallback definition switch; the chosen
#' branch is recorded in each metric's \code{definition_label} and
#' the orchestrator mirrors the choice into
#' \code{metadata.definitions.pathway} / \code{.reaction}.
#'
#' Routes to \code{dev5} via \code{\link{pg_query_panel}}
#' (post-2026-06-14 dev5 integrated-build promotion).
#'
#' @param panel_id Character: panel id for deployment routing.
#'     Defaults to \code{"fig01-architecture"}.
#' @param definitions Named list: from
#'     \code{\link{digest_definitions}}; \code{NULL} → load fresh.
#' @param runtime Named list: from \code{\link{digest_runtime}};
#'     \code{NULL} → load fresh.
#'
#' @return Tibble of DigestMetric rows. Carries
#'     \code{attr(., "queries")} and \code{attr(., "branches")} (a
#'     named list with \code{pathway} and \code{reaction} chosen
#'     strings the orchestrator uses to build
#'     \code{metadata.definitions}).
#'
#' @importFrom logger log_info
#' @importFrom tibble tibble
#' @export
section_interactions <- function(
    panel_id = "fig01-architecture",
    definitions = NULL,
    runtime = NULL
) {

    defs <- definitions %||% digest_definitions()
    rt   <- runtime     %||% digest_runtime()
    facet <- "panel_a_stats_interactions"

    totals_sql <- interactions_totals_sql()
    totals_rows <- pg_query_panel(panel_id, totals_sql, facet = facet)

    n_int <- as.integer(totals_rows$n_interactions %||% 0L)
    n_pog <- as.integer(totals_rows$n_proteins_or_genes %||% 0L)

    pathway_result <- resolve_pathway_count(
        panel_id,
        defs,
        rt,
        facet
    )
    reaction_result <- resolve_reaction_count(
        panel_id,
        defs,
        rt,
        facet
    )

    metrics <- dplyr::bind_rows(
        tibble::tibble(
            section_id       = 3L,
            metric_name      = "n_interactions",
            value            = n_int,
            definition_label = "relation.category:interaction",
            state            = if (n_int > 0L) "populated" else "empty",
            deployment       = attr(totals_rows, "deployment"),
            sql_hash         = sql_sha256(totals_sql)
        ),
        tibble::tibble(
            section_id       = 3L,
            metric_name      = "n_proteins_or_genes",
            value            = n_pog,
            definition_label = "interactions:distinct_protein_or_gene",
            state            = if (n_pog > 0L) "populated" else "empty",
            deployment       = attr(totals_rows, "deployment"),
            sql_hash         = sql_sha256(totals_sql)
        ),
        tibble::tibble(
            section_id       = 3L,
            metric_name      = "n_pathways",
            value            = pathway_result$value,
            definition_label = pathway_result$definition_label,
            state            = if (pathway_result$value > 0L) {
                "populated"
            } else {
                "empty"
            },
            deployment       = pathway_result$deployment,
            sql_hash         = pathway_result$sql_hash
        ),
        tibble::tibble(
            section_id       = 3L,
            metric_name      = "n_reactions",
            value            = reaction_result$value,
            definition_label = reaction_result$definition_label,
            state            = if (reaction_result$value > 0L) {
                "populated"
            } else {
                "empty"
            },
            deployment       = reaction_result$deployment,
            sql_hash         = reaction_result$sql_hash
        )
    )

    queries <- list(
        totals = list(
            metric_names = c("n_interactions", "n_proteins_or_genes"),
            sql          = totals_sql,
            sql_hash     = sql_sha256(totals_sql),
            deployment   = attr(totals_rows, "deployment"),
            result_hash  = attr(totals_rows, "result_hash"),
            row_count    = nrow(totals_rows)
        ),
        pathway = pathway_result$query,
        reaction = reaction_result$query
    )

    branches <- list(
        pathway  = pathway_result$branch,
        reaction = reaction_result$branch
    )

    attr(metrics, "queries")  <- queries
    attr(metrics, "branches") <- branches

    logger::log_info(
        "section_interactions: n_int={n_int}, n_pog={n_pog}, ",
        "pathway[{pathway_result$branch}]={pathway_result$value}, ",
        "reaction[{reaction_result$branch}]={reaction_result$value}"
    )

    metrics
}


#' Totals SQL — distinct interaction relations + distinct protein/gene
#' partners
#'
#' @return Character.
#'
#' @keywords internal
#' @noRd
interactions_totals_sql <- function() {
    # Materialize the interaction_cat lookup + pog set so Postgres
    # can use the relation_category_subject_idx /
    # relation_category_object_idx indexes for the count + the
    # protein-participation probe. The protein-participation count
    # unions subject and object roles then dedups, mirroring the
    # MPI base-set pattern in stats-section_mpi.R.
    paste(
        "WITH interaction_cat AS (",
        "    SELECT relation_category_id FROM vocab_relation_category",
        "     WHERE name = 'interaction'",
        "), pog AS MATERIALIZED (",
        "    SELECT entity_id FROM entity e",
        "      JOIN vocab_entity_type vet USING (entity_type_id)",
        "     WHERE vet.name IN ('Gene:MI:0250','Protein:MI:0326')",
        "), participating_pog AS MATERIALIZED (",
        "    SELECT DISTINCT r.subject_entity_id AS entity_id",
        "      FROM relation r",
        "      JOIN pog ON pog.entity_id = r.subject_entity_id",
        "     WHERE r.relation_category_id IN",
        "             (SELECT relation_category_id FROM interaction_cat)",
        "    UNION",
        "    SELECT DISTINCT r.object_entity_id AS entity_id",
        "      FROM relation r",
        "      JOIN pog ON pog.entity_id = r.object_entity_id",
        "     WHERE r.relation_category_id IN",
        "             (SELECT relation_category_id FROM interaction_cat)",
        ")",
        "SELECT",
        "    (SELECT COUNT(*)::bigint FROM relation r",
        "      WHERE r.relation_category_id IN",
        "              (SELECT relation_category_id FROM interaction_cat))",
        "        AS n_interactions,",
        "    (SELECT COUNT(*)::bigint FROM participating_pog)",
        "        AS n_proteins_or_genes",
        sep = "\n"
    )
}


#' Resolve the pathway-count metric (preferred branch with fallback)
#'
#' Issues the preferred query first; if the returned count is below
#' \code{runtime.pathway.min_rows_for_preferred} the fallback query
#' fires. The returned record carries the chosen value, the chosen
#' definition label, the sql hash of the chosen query, the chosen
#' branch string (the value mirrored into
#' \code{metadata.definitions.pathway}), the deployment, and the
#' query record for sidecar transcription.
#'
#' @param panel_id Character.
#' @param definitions Named list.
#' @param runtime Named list.
#' @param facet Character: facet id for \code{pg_query_panel}.
#'
#' @return Named list:
#'     \code{value, definition_label, sql_hash, deployment,
#'     branch, query}.
#'
#' @keywords internal
#' @noRd
resolve_pathway_count <- function(panel_id, definitions, runtime, facet) {

    pref_sql <- pathway_preferred_sql()
    pref_rows <- pg_query_panel(panel_id, pref_sql, facet = facet)
    pref_value <- as.integer(pref_rows$n_pathways %||% 0L)

    min_rows <- as.integer(runtime$pathway$min_rows_for_preferred %||% 10L)

    if (pref_value >= min_rows) {
        return(list(
            value            = pref_value,
            definition_label = sprintf(
                "pathway:%s",
                definitions$pathway$preferred
            ),
            sql_hash         = sql_sha256(pref_sql),
            deployment       = attr(pref_rows, "deployment"),
            branch           = definitions$pathway$preferred,
            query = list(
                metric_names = "n_pathways",
                sql          = pref_sql,
                sql_hash     = sql_sha256(pref_sql),
                deployment   = attr(pref_rows, "deployment"),
                result_hash  = attr(pref_rows, "result_hash"),
                row_count    = nrow(pref_rows),
                branch       = definitions$pathway$preferred
            )
        ))
    }

    fb_sql <- pathway_annotation_sql()
    fb_rows <- pg_query_panel(panel_id, fb_sql, facet = facet)
    fb_value <- as.integer(fb_rows$n_pathways %||% 0L)

    list(
        value            = fb_value,
        definition_label = sprintf(
            "pathway:%s",
            definitions$pathway$fallback
        ),
        sql_hash         = sql_sha256(fb_sql),
        deployment       = attr(fb_rows, "deployment"),
        branch           = definitions$pathway$fallback,
        query = list(
            metric_names = "n_pathways",
            sql          = fb_sql,
            sql_hash     = sql_sha256(fb_sql),
            deployment   = attr(fb_rows, "deployment"),
            result_hash  = attr(fb_rows, "result_hash"),
            row_count    = nrow(fb_rows),
            branch       = definitions$pathway$fallback
        )
    )
}


#' Pathway preferred SQL — distinct Pathway:OM:0014 entities that
#' appear as either subject or object in any interaction relation
#'
#' Simpler shape: count distinct \code{Pathway:OM:0014} entities
#' that participate (as subject or object) in any relation of
#' category \code{interaction}. The original "with >= 1 member
#' protein/gene" criterion is dropped because the live schema does
#' not maintain a separate pathway-membership relation type — every
#' interaction-class relation involving a pathway entity already
#' implies a protein/gene-anchored connection.
#'
#' @return Character.
#'
#' @keywords internal
#' @noRd
pathway_preferred_sql <- function() {
    # The OR in EXISTS forces a sequential scan because Postgres
    # can't pick between the (cat_id, subject) and (cat_id, object)
    # indexes for one predicate. Materialize "entities that
    # participate in any interaction" once (two index scans + a
    # UNION) and then intersect with Pathway:OM:0014.
    paste(
        "WITH ic AS (",
        "    SELECT relation_category_id FROM vocab_relation_category",
        "     WHERE name = 'interaction'",
        "), interaction_participants AS MATERIALIZED (",
        "    SELECT DISTINCT subject_entity_id AS entity_id FROM relation",
        "     WHERE relation_category_id IN (SELECT relation_category_id FROM ic)",
        "    UNION",
        "    SELECT DISTINCT object_entity_id AS entity_id FROM relation",
        "     WHERE relation_category_id IN (SELECT relation_category_id FROM ic)",
        ")",
        "SELECT COUNT(*)::bigint AS n_pathways",
        "  FROM entity e",
        "  JOIN vocab_entity_type vet USING (entity_type_id)",
        " WHERE vet.name = 'Pathway:OM:0014'",
        "   AND e.entity_id IN (SELECT entity_id FROM interaction_participants)",
        sep = "\n"
    )
}


#' Pathway fallback SQL — distinct pathway-annotation terms
#'
#' Counts distinct \code{entity_ontology_term.term_id} values whose
#' \code{ontology_prefix} matches a pathway-bearing ontology
#' (Reactome IDs all start with \code{r-hsa} / \code{r-mmu} / etc.
#' on dev5; KEGG appears as a data_source not an ontology prefix
#' so it is not counted here).
#'
#' @return Character.
#'
#' @keywords internal
#' @noRd
pathway_annotation_sql <- function() {
    paste(
        "SELECT COUNT(DISTINCT term_id)::bigint AS n_pathways",
        "  FROM entity_ontology_term",
        " WHERE ontology_prefix LIKE 'r-%'",
        sep = "\n"
    )
}


#' Resolve the reaction-count metric (preferred branch with fallback)
#'
#' Mirrors \code{\link{resolve_pathway_count}} for the reaction
#' definition.
#'
#' @param panel_id Character.
#' @param definitions Named list.
#' @param runtime Named list.
#' @param facet Character.
#'
#' @return Named list.
#'
#' @keywords internal
#' @noRd
resolve_reaction_count <- function(panel_id, definitions, runtime, facet) {

    pref_sql <- reaction_preferred_sql()
    pref_rows <- pg_query_panel(panel_id, pref_sql, facet = facet)
    pref_value <- as.integer(pref_rows$n_reactions %||% 0L)

    min_rows <- as.integer(
        runtime$reaction$min_rows_for_preferred %||% 10L
    )

    if (pref_value >= min_rows) {
        return(list(
            value            = pref_value,
            definition_label = sprintf(
                "reaction:%s",
                definitions$reaction$preferred
            ),
            sql_hash         = sql_sha256(pref_sql),
            deployment       = attr(pref_rows, "deployment"),
            branch           = definitions$reaction$preferred,
            query = list(
                metric_names = "n_reactions",
                sql          = pref_sql,
                sql_hash     = sql_sha256(pref_sql),
                deployment   = attr(pref_rows, "deployment"),
                result_hash  = attr(pref_rows, "result_hash"),
                row_count    = nrow(pref_rows),
                branch       = definitions$reaction$preferred
            )
        ))
    }

    fb_sql <- reaction_predicate_sql()
    fb_rows <- pg_query_panel(panel_id, fb_sql, facet = facet)
    fb_value <- as.integer(fb_rows$n_reactions %||% 0L)

    list(
        value            = fb_value,
        definition_label = sprintf(
            "reaction:%s",
            definitions$reaction$fallback
        ),
        sql_hash         = sql_sha256(fb_sql),
        deployment       = attr(fb_rows, "deployment"),
        branch           = definitions$reaction$fallback,
        query = list(
            metric_names = "n_reactions",
            sql          = fb_sql,
            sql_hash     = sql_sha256(fb_sql),
            deployment   = attr(fb_rows, "deployment"),
            result_hash  = attr(fb_rows, "result_hash"),
            row_count    = nrow(fb_rows),
            branch       = definitions$reaction$fallback
        )
    )
}


#' Reaction preferred SQL — distinct Reaction:OM:0015 entities
#' participating in any interaction relation
#'
#' Same simplification as
#' \code{\link{pathway_preferred_sql}}.
#'
#' @return Character.
#'
#' @keywords internal
#' @noRd
reaction_preferred_sql <- function() {
    # Same UNION strategy as pathway_preferred_sql.
    paste(
        "WITH ic AS (",
        "    SELECT relation_category_id FROM vocab_relation_category",
        "     WHERE name = 'interaction'",
        "), interaction_participants AS MATERIALIZED (",
        "    SELECT DISTINCT subject_entity_id AS entity_id FROM relation",
        "     WHERE relation_category_id IN (SELECT relation_category_id FROM ic)",
        "    UNION",
        "    SELECT DISTINCT object_entity_id AS entity_id FROM relation",
        "     WHERE relation_category_id IN (SELECT relation_category_id FROM ic)",
        ")",
        "SELECT COUNT(*)::bigint AS n_reactions",
        "  FROM entity e",
        "  JOIN vocab_entity_type vet USING (entity_type_id)",
        " WHERE vet.name = 'Reaction:OM:0015'",
        "   AND e.entity_id IN (SELECT entity_id FROM interaction_participants)",
        sep = "\n"
    )
}


#' Reaction fallback SQL — distinct relations whose predicate
#' classifies as a metabolic reaction (orthosteric / enzyme-substrate
#' edges from BRENDA / KEGG / Rhea)
#'
#' @return Character.
#'
#' @keywords internal
#' @noRd
reaction_predicate_sql <- function() {
    paste(
        "SELECT COUNT(DISTINCT r.relation_id)::bigint AS n_reactions",
        "  FROM relation r",
        "  JOIN vocab_relation_category vrc",
        "    ON vrc.relation_category_id = r.relation_category_id",
        "  JOIN relation_evidence_relation rer",
        "    ON rer.relation_id = r.relation_id",
        "  JOIN data_source ds",
        "    ON ds.source_id = rer.source_id",
        " WHERE vrc.name = 'interaction'",
        "   AND ds.name IN ('brenda','kegg','rhea','recon3d')",
        sep = "\n"
    )
}
