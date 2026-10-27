#' SQL CTE for the metabolite-protein interaction (MPI) base set
#'
#' Returns the SQL fragment that defines the \code{chem}, \code{pog},
#' and \code{mpi} CTEs used by every Section-2 metric. The MPI base
#' set is the subset of \code{relation} rows whose
#' \code{vocab_relation_category = 'interaction'} AND one participant
#' is a \code{Chemical:OM:0037} entity AND the other is a protein /
#' gene entity (\code{Gene:MI:0250} ∪ \code{Protein:MI:0326}).
#'
#' The two orientations (subject=chem ∧ object=pog and
#' subject=pog ∧ object=chem) are captured as a UNION rather than a
#' WHERE-OR so Postgres can use the
#' \code{relation_category_subject_idx} / \code{...object_idx}
#' indexes as two index scans. The \code{mpi} CTE is \code{MATERIALIZED}
#' so the totals + transporter + receptor metrics share one scan
#' rather than re-deriving the same row set three times.
#'
#' @return Character: SQL fragment ending after the CTE definitions
#'     (no trailing comma; the caller appends additional CTEs or the
#'     terminal SELECT).
#'
#' @keywords internal
#' @noRd
mpi_base_cte_sql <- function() {
    paste(
        "WITH chem AS MATERIALIZED (",
        "    SELECT entity_id FROM entity e",
        "      JOIN vocab_entity_type vet USING (entity_type_id)",
        "     WHERE vet.name = 'Chemical:OM:0037'",
        "), pog AS MATERIALIZED (",
        "    SELECT entity_id FROM entity e",
        "      JOIN vocab_entity_type vet USING (entity_type_id)",
        "     WHERE vet.name IN ('Gene:MI:0250','Protein:MI:0326')",
        "), interaction_cat AS (",
        "    SELECT relation_category_id FROM vocab_relation_category",
        "     WHERE name = 'interaction'",
        "), mpi AS MATERIALIZED (",
        "    SELECT r.relation_id,",
        "           r.subject_entity_id AS chem_id,",
        "           r.object_entity_id  AS prot_id,",
        "           r.predicate_id",
        "      FROM relation r",
        "      JOIN chem ON chem.entity_id = r.subject_entity_id",
        "      JOIN pog  ON pog.entity_id  = r.object_entity_id",
        "     WHERE r.relation_category_id IN (SELECT relation_category_id",
        "                                        FROM interaction_cat)",
        "    UNION ALL",
        "    SELECT r.relation_id,",
        "           r.object_entity_id  AS chem_id,",
        "           r.subject_entity_id AS prot_id,",
        "           r.predicate_id",
        "      FROM relation r",
        "      JOIN chem ON chem.entity_id = r.object_entity_id",
        "      JOIN pog  ON pog.entity_id  = r.subject_entity_id",
        "     WHERE r.relation_category_id IN (SELECT relation_category_id",
        "                                        FROM interaction_cat)",
        ")",
        sep = "\n"
    )
}


#' Build the combined Section-2 SQL — all five metrics in one query
#'
#' Bundles the totals (n_mpi_relations / n_metabolites /
#' n_proteins_or_genes) and the transporter / receptor sub-counts
#' into a single materialized-CTE query so the digest builds the
#' MPI base set ONCE per rebuild rather than three times. The
#' transporter / receptor classifications use the curated
#' \code{data_source.name} lists (\code{definitions.transporter.resources}
#' and \code{definitions.receptor.resources}) plus optional UniProt
#' keywords; empty UniProt-keyword lists disable that branch
#' cleanly.
#'
#' @param tra_resources Character vector.
#' @param tra_keywords Character vector (possibly empty).
#' @param rec_resources Character vector.
#' @param rec_keywords Character vector (possibly empty).
#'
#' @return Character: complete SQL string.
#'
#' @keywords internal
#' @noRd
mpi_combined_sql <- function(
    tra_resources,
    tra_keywords,
    rec_resources,
    rec_keywords
) {

    # Pre-compute the small (~33K rows) transporter / receptor
    # candidate-entity sets ONCE so the 2.8M-row MPI scan can
    # intersect against them via in-set lookup rather than a
    # per-row EXISTS scan into entity_evidence_resolution.
    paste(
        mpi_base_cte_sql(),
        ", transporter_resource_set AS MATERIALIZED (",
        "    SELECT DISTINCT eer.entity_id",
        "      FROM entity_evidence_resolution eer",
        "      JOIN entity_evidence ee USING (entity_evidence_id)",
        "      JOIN data_source ds ON ds.source_id = ee.source_id",
        sprintf("     WHERE ds.name IN (%s)",
                sql_in_list(tra_resources)),
        "), receptor_resource_set AS MATERIALIZED (",
        "    SELECT DISTINCT eer.entity_id",
        "      FROM entity_evidence_resolution eer",
        "      JOIN entity_evidence ee USING (entity_evidence_id)",
        "      JOIN data_source ds ON ds.source_id = ee.source_id",
        sprintf("     WHERE ds.name IN (%s)",
                sql_in_list(rec_resources)),
        "), transport_predicates AS MATERIALIZED (",
        "    SELECT vrp.relation_predicate_id",
        "      FROM vocab_relation_predicate vrp",
        "      JOIN vocab_interaction_class vic",
        "        ON vic.interaction_class_id = vrp.interaction_class_id",
        "     WHERE vic.name = 'Transport'",
        ")",
        "SELECT",
        "    COUNT(DISTINCT mpi.relation_id)::bigint AS n_mpi_relations,",
        "    COUNT(DISTINCT mpi.chem_id)::bigint     AS n_metabolites,",
        "    COUNT(DISTINCT mpi.prot_id)::bigint     AS n_proteins_or_genes,",
        "    COUNT(DISTINCT mpi.prot_id) FILTER (WHERE",
        "        mpi.predicate_id IN (SELECT relation_predicate_id",
        "                              FROM transport_predicates)",
        "     OR mpi.prot_id IN (SELECT entity_id FROM transporter_resource_set)",
        "    )::bigint AS n_transporters,",
        "    COUNT(DISTINCT mpi.prot_id) FILTER (WHERE",
        "        mpi.prot_id IN (SELECT entity_id FROM receptor_resource_set)",
        "    )::bigint AS n_receptors",
        "  FROM mpi",
        sep = "\n"
    )
}


#' Build the totals SQL for the Section-2 MPI counts
#'
#' Reuses \code{\link{mpi_base_cte_sql}} for the CTE; produces three
#' headline counts: total MPI relations, distinct metabolites, and
#' distinct protein/gene partners. The transporter / receptor
#' sub-counts are computed by separate queries
#' (\code{\link{mpi_transporter_sql}}, \code{\link{mpi_receptor_sql}})
#' so the curated-resource list can be parameterised cleanly.
#'
#' @return Character: complete SQL string.
#'
#' @keywords internal
#' @noRd
mpi_totals_sql <- function() {
    paste(
        mpi_base_cte_sql(),
        "SELECT",
        "    COUNT(DISTINCT relation_id)::bigint AS n_mpi_relations,",
        "    COUNT(DISTINCT chem_id)::bigint     AS n_metabolites,",
        "    COUNT(DISTINCT prot_id)::bigint     AS n_proteins_or_genes",
        "  FROM mpi",
        sep = "\n"
    )
}


#' Build the transporter sub-count SQL
#'
#' Counts distinct protein/gene partners in the MPI base set that
#' satisfy at least one of:
#' \enumerate{
#'   \item the partnering predicate's
#'         \code{vocab_relation_predicate.interaction_class_id}
#'         resolves to \code{Transport} (via
#'         \code{vocab_interaction_class.name}).
#'   \item the protein partner is sourced by a
#'         \code{data_source.name} in
#'         \code{definitions.transporter.resources}.
#'   \item the protein partner carries a UniProt keyword in
#'         \code{definitions.transporter.uniprot_keywords}
#'         (via the entity-ontology join — empty list disables this
#'         branch).
#' }
#'
#' @param resources Character vector: curated transporter
#'     \code{data_source.name} list.
#' @param uniprot_keywords Character vector: UniProt keyword codes
#'     (e.g. \code{c("KW-0813", "KW-0769", "KW-0917")}). Empty
#'     vector disables this branch.
#'
#' @return Character: complete SQL string.
#'
#' @importFrom DBI dbQuoteString ANSI
#' @keywords internal
#' @noRd
mpi_transporter_sql <- function(resources, uniprot_keywords) {
    paste(
        mpi_base_cte_sql(),
        "SELECT COUNT(DISTINCT mpi.prot_id)::bigint AS n_transporters",
        "  FROM mpi",
        " WHERE",
        "      EXISTS (",
        "          SELECT 1 FROM vocab_relation_predicate vrp",
        "            JOIN vocab_interaction_class vic",
        "              ON vic.interaction_class_id = vrp.interaction_class_id",
        "           WHERE vrp.relation_predicate_id = mpi.predicate_id",
        "             AND vic.name = 'Transport'",
        "      )",
        "   OR EXISTS (",
        "          SELECT 1 FROM entity_evidence_resolution eer",
        "            JOIN entity_evidence ee USING (entity_evidence_id)",
        "            JOIN data_source ds ON ds.source_id = ee.source_id",
        "           WHERE eer.entity_id = mpi.prot_id",
        sprintf(
            "             AND ds.name IN (%s)",
            sql_in_list(resources)
        ),
        "      )",
        uniprot_keyword_clause(uniprot_keywords),
        sep = "\n"
    )
}


#' Build the receptor sub-count SQL
#'
#' Counts distinct protein/gene partners in the MPI base set that
#' satisfy at least one of:
#' \enumerate{
#'   \item the protein partner is sourced by a
#'         \code{data_source.name} in
#'         \code{definitions.receptor.resources}.
#'   \item the protein partner carries a UniProt keyword in
#'         \code{definitions.receptor.uniprot_keywords}.
#' }
#'
#' Per FR-043b, \code{interaction_class_id} alone is NOT sufficient
#' for receptors — the cycle-001 predicate vocabulary collapsed
#' Ligand-receptor into \code{Other}; the dev5 integrated build may
#' expand this but the curated-resource path is the primary
#' classification route for digest stability.
#'
#' @param resources Character vector: curated receptor
#'     \code{data_source.name} list.
#' @param uniprot_keywords Character vector: UniProt keyword codes.
#'
#' @return Character: complete SQL string.
#'
#' @keywords internal
#' @noRd
mpi_receptor_sql <- function(resources, uniprot_keywords) {
    paste(
        mpi_base_cte_sql(),
        "SELECT COUNT(DISTINCT mpi.prot_id)::bigint AS n_receptors",
        "  FROM mpi",
        " WHERE",
        "      EXISTS (",
        "          SELECT 1 FROM entity_evidence_resolution eer",
        "            JOIN entity_evidence ee USING (entity_evidence_id)",
        "            JOIN data_source ds ON ds.source_id = ee.source_id",
        "           WHERE eer.entity_id = mpi.prot_id",
        sprintf(
            "             AND ds.name IN (%s)",
            sql_in_list(resources)
        ),
        "      )",
        uniprot_keyword_clause(uniprot_keywords),
        sep = "\n"
    )
}


#' Append the UniProt-keyword EXISTS branch to a transporter/receptor SQL
#'
#' Returns an empty string when no keywords are configured (the
#' caller still includes the preceding \code{OR} via the resource
#' branch, so the result is a syntactically valid WHERE clause).
#' The keyword join targets the entity-ontology shape used by the
#' integrated build; if the build labels UniProt keywords under a
#' different ontology source the orchestrator's first run on beauty
#' surfaces it as a clear missing-table error rather than a silent
#' undercount.
#'
#' @param uniprot_keywords Character vector.
#'
#' @return Character: SQL fragment (possibly empty).
#'
#' @keywords internal
#' @noRd
uniprot_keyword_clause <- function(uniprot_keywords) {
    if (length(uniprot_keywords) == 0L) {
        return("")
    }
    # dev5 stores UniProt keywords under ontology_prefix = 'uniprot'
    # (no separate keyword ontology); the term_id carries the
    # KW-NNNN code.
    paste(
        "   OR EXISTS (",
        "          SELECT 1 FROM entity_ontology_term eot",
        "           WHERE eot.term_entity_id = mpi.prot_id",
        "             AND eot.ontology_prefix = 'uniprot'",
        sprintf(
            "             AND eot.term_id IN (%s)",
            sql_in_list(uniprot_keywords)
        ),
        "      )",
        sep = "\n"
    )
}


#' Render a character vector as a parenthesised SQL IN-list literal
#'
#' Empty vector → \code{NULL} sentinel so the surrounding
#' \code{<col> IN (NULL)} clause returns false rather than failing
#' parse.
#'
#' @param values Character vector.
#'
#' @return Character: comma-separated quoted-string list.
#'
#' @keywords internal
#' @noRd
sql_in_list <- function(values) {
    if (length(values) == 0L) {
        return("NULL")
    }
    paste(
        sprintf("'%s'", gsub("'", "''", values, fixed = TRUE)),
        collapse = ", "
    )
}


#' Build a human-readable definition label for a curated classifier
#'
#' Encodes the resource list + keyword presence flag into a stable
#' string the orchestrator stamps onto each metric's
#' \code{definition_label} so a reader of the digest can recover the
#' exact classification recipe.
#'
#' @param prefix \code{"transporter"} or \code{"receptor"}.
#' @param resources Character vector.
#' @param uniprot_keywords Character vector.
#'
#' @return Character.
#'
#' @keywords internal
#' @noRd
classifier_definition_label <- function(
    prefix,
    resources,
    uniprot_keywords
) {
    parts <- c(
        paste(resources, collapse = "+"),
        if (length(uniprot_keywords) > 0L) "uniprot_kw" else NULL
    )
    sprintf("%s:%s", prefix, paste(parts, collapse = "+"))
}


#' Section 2 — Metabolite-Protein Interactions (FR-043b)
#'
#' Builds the MPI base set, derives the five Section-2 metrics, and
#' asserts the FR-043b subset constraint
#' (\code{n_transporters ≤ n_proteins_or_genes} AND
#' \code{n_receptors ≤ n_proteins_or_genes}). The chosen
#' transporter / receptor definition labels are mirrored into each
#' metric's \code{definition_label} field.
#'
#' Routes to \code{dev5} via \code{\link{pg_query_panel}}
#' (post-2026-06-14 dev5 integrated-build promotion).
#'
#' @param panel_id Character: panel id for deployment routing.
#'     Defaults to \code{"architecture"}.
#' @param definitions Named list: the \code{definitions:} namespace
#'     from \code{\link{digest_definitions}}. \code{NULL} (default) →
#'     loaded via \code{\link{digest_definitions}()}.
#'
#' @return Tibble of DigestMetric rows (cols
#'     \code{section_id, metric_name, value, definition_label,
#'     state, deployment, sql_hash}). Carries
#'     \code{attr(., "queries")} for sidecar transcription.
#'
#' @importFrom logger log_info
#' @importFrom rlang abort
#' @importFrom tibble tibble
#' @export
section_mpi <- function(
    panel_id = "architecture",
    definitions = NULL
) {

    defs <- definitions %||% digest_definitions()
    facet <- "panel_a_stats_mpi"

    combined_sql <- mpi_combined_sql(
        defs$transporter$resources,
        defs$transporter$uniprot_keywords,
        defs$receptor$resources,
        defs$receptor$uniprot_keywords
    )

    rows <- pg_query_panel(panel_id, combined_sql, facet = facet)

    n_mpi <- as.integer(rows$n_mpi_relations %||% 0L)
    n_met <- as.integer(rows$n_metabolites %||% 0L)
    n_pog <- as.integer(rows$n_proteins_or_genes %||% 0L)
    n_tra <- as.integer(rows$n_transporters %||% 0L)
    n_rec <- as.integer(rows$n_receptors %||% 0L)

    # FR-043b subset assertion.
    if (n_tra > n_pog || n_rec > n_pog) {
        rlang::abort(sprintf(
            paste0(
                "FR-043b subset constraint violated: ",
                "n_transporters=%d, n_receptors=%d, ",
                "n_proteins_or_genes=%d. Both subset counts MUST be ",
                "<= the protein/gene count."
            ),
            n_tra,
            n_rec,
            n_pog
        ))
    }

    tra_label <- classifier_definition_label(
        "transporter",
        defs$transporter$resources,
        defs$transporter$uniprot_keywords
    )
    rec_label <- classifier_definition_label(
        "receptor",
        defs$receptor$resources,
        defs$receptor$uniprot_keywords
    )

    deployment <- attr(rows, "deployment")
    sql_h <- sql_sha256(combined_sql)

    metrics <- dplyr::bind_rows(
        tibble::tibble(
            section_id       = 2L,
            metric_name      = "n_mpi_relations",
            value            = n_mpi,
            definition_label = "mpi_base_cte:chemical+protein_or_gene",
            state            = if (n_mpi > 0L) "populated" else "empty",
            deployment       = deployment,
            sql_hash         = sql_h
        ),
        tibble::tibble(
            section_id       = 2L,
            metric_name      = "n_metabolites",
            value            = n_met,
            definition_label = "mpi_base_cte:distinct_chem_id",
            state            = if (n_met > 0L) "populated" else "empty",
            deployment       = deployment,
            sql_hash         = sql_h
        ),
        tibble::tibble(
            section_id       = 2L,
            metric_name      = "n_proteins_or_genes",
            value            = n_pog,
            definition_label = "mpi_base_cte:distinct_prot_id",
            state            = if (n_pog > 0L) "populated" else "empty",
            deployment       = deployment,
            sql_hash         = sql_h
        ),
        tibble::tibble(
            section_id       = 2L,
            metric_name      = "n_transporters",
            value            = n_tra,
            definition_label = tra_label,
            state            = if (n_tra > 0L) "populated" else "empty",
            deployment       = deployment,
            sql_hash         = sql_h
        ),
        tibble::tibble(
            section_id       = 2L,
            metric_name      = "n_receptors",
            value            = n_rec,
            definition_label = rec_label,
            state            = if (n_rec > 0L) "populated" else "empty",
            deployment       = deployment,
            sql_hash         = sql_h
        )
    )

    queries <- list(
        combined = list(
            metric_names = c(
                "n_mpi_relations",
                "n_metabolites",
                "n_proteins_or_genes",
                "n_transporters",
                "n_receptors"
            ),
            sql          = combined_sql,
            sql_hash     = sql_h,
            deployment   = deployment,
            result_hash  = attr(rows, "result_hash"),
            row_count    = nrow(rows)
        )
    )

    attr(metrics, "queries") <- queries

    logger::log_info(
        "section_mpi: n_mpi={n_mpi}, n_met={n_met}, ",
        "n_pog={n_pog}, n_tra={n_tra}, n_rec={n_rec}"
    )

    metrics
}
