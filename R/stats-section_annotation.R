#' Section 5 — Annotation (FR-043e)
#'
#' Five-slot list of annotation statistics. Slot 3 carries a
#' fallback chain — default \code{diseases} (MONDO + DOID terms);
#' if the disease count falls below
#' \code{runtime.annotation.disease_min_rows} the orchestrator walks
#' \code{runtime.annotation.slot3_fallback_order} until one
#' category meets the threshold. The chosen slot is mirrored into
#' \code{metadata.definitions.annotation_slot3} and the slot 3 metric
#' carries \code{state = "fallback"} when it is not the default.
#'
#' Routes to \code{dev5} via \code{\link{pg_query_panel}}
#' (post-2026-06-14 dev5 integrated-build promotion).
#'
#' @param panel_id Character.
#' @param runtime Named list: from \code{\link{digest_runtime}};
#'     \code{NULL} → load fresh.
#'
#' @return Tibble of DigestMetric rows. Carries
#'     \code{attr(., "queries")} and \code{attr(., "branches")}
#'     (named list with \code{annotation_slot3} string for
#'     \code{metadata.definitions} mirror).
#'
#' @importFrom logger log_info
#' @importFrom tibble tibble
#' @export
section_annotation <- function(
    panel_id = "architecture",
    runtime = NULL
) {

    rt <- runtime %||% digest_runtime()
    facet <- "panel_a_stats_annotation"

    records_sql <- annotation_records_sql()
    records_rows <- pg_query_panel(panel_id, records_sql, facet = facet)
    n_records <- as.integer(records_rows$n_annotation_records %||% 0L)

    organisms_sql <- annotation_organisms_sql()
    organisms_rows <- pg_query_panel(panel_id, organisms_sql, facet = facet)
    n_organisms <- as.integer(organisms_rows$n_organisms %||% 0L)

    slot3 <- resolve_annotation_slot3(panel_id, rt, facet)

    loc_sql <- localization_sql()
    loc_rows <- pg_query_panel(panel_id, loc_sql, facet = facet)
    n_loc <- as.integer(loc_rows$n_localizations %||% 0L)

    pw_sql <- pathway_annotation_sql()
    pw_rows <- pg_query_panel(panel_id, pw_sql, facet = facet)
    n_pw <- as.integer(pw_rows$n_pathways %||% 0L)

    metrics <- dplyr::bind_rows(
        tibble::tibble(
            section_id       = 5L,
            metric_name      = "n_annotation_records",
            value            = n_records,
            definition_label = "entity_evidence_annotation:total",
            state            = if (n_records > 0L) "populated" else "empty",
            deployment       = attr(records_rows, "deployment"),
            sql_hash         = sql_sha256(records_sql)
        ),
        tibble::tibble(
            section_id       = 5L,
            metric_name      = "n_organisms",
            value            = n_organisms,
            definition_label = "entity_type:Organism:OM:0032+annotated",
            state            = if (n_organisms > 0L) "populated" else "empty",
            deployment       = attr(organisms_rows, "deployment"),
            sql_hash         = sql_sha256(organisms_sql)
        ),
        tibble::tibble(
            section_id       = 5L,
            metric_name      = slot3$metric_name,
            value            = slot3$value,
            definition_label = slot3$definition_label,
            state            = slot3$state,
            deployment       = slot3$deployment,
            sql_hash         = slot3$sql_hash
        ),
        tibble::tibble(
            section_id       = 5L,
            metric_name      = "n_localizations",
            value            = n_loc,
            definition_label = "entity_ontology_term:uberon_anatomy",
            state            = if (n_loc > 0L) "populated" else "empty",
            deployment       = attr(loc_rows, "deployment"),
            sql_hash         = sql_sha256(loc_sql)
        ),
        tibble::tibble(
            section_id       = 5L,
            metric_name      = "n_pathways_annotation",
            value            = n_pw,
            definition_label = "entity_ontology_term:pathway_annotations",
            state            = if (n_pw > 0L) "populated" else "empty",
            deployment       = attr(pw_rows, "deployment"),
            sql_hash         = sql_sha256(pw_sql)
        )
    )

    queries <- list(
        records = list(
            metric_names = "n_annotation_records",
            sql          = records_sql,
            sql_hash     = sql_sha256(records_sql),
            deployment   = attr(records_rows, "deployment"),
            result_hash  = attr(records_rows, "result_hash"),
            row_count    = nrow(records_rows)
        ),
        organisms = list(
            metric_names = "n_organisms",
            sql          = organisms_sql,
            sql_hash     = sql_sha256(organisms_sql),
            deployment   = attr(organisms_rows, "deployment"),
            result_hash  = attr(organisms_rows, "result_hash"),
            row_count    = nrow(organisms_rows)
        ),
        slot3 = slot3$query,
        localization = list(
            metric_names = "n_localizations",
            sql          = loc_sql,
            sql_hash     = sql_sha256(loc_sql),
            deployment   = attr(loc_rows, "deployment"),
            result_hash  = attr(loc_rows, "result_hash"),
            row_count    = nrow(loc_rows)
        ),
        pathway_annotation = list(
            metric_names = "n_pathways_annotation",
            sql          = pw_sql,
            sql_hash     = sql_sha256(pw_sql),
            deployment   = attr(pw_rows, "deployment"),
            result_hash  = attr(pw_rows, "result_hash"),
            row_count    = nrow(pw_rows)
        )
    )

    branches <- list(annotation_slot3 = slot3$branch)

    attr(metrics, "queries")  <- queries
    attr(metrics, "branches") <- branches

    logger::log_info(
        "section_annotation: records={n_records}, ",
        "organisms={n_organisms}, slot3={slot3$branch}={slot3$value}, ",
        "localizations={n_loc}, pathways={n_pw}"
    )

    metrics
}


#' Totals SQL — annotation records + distinct organism entities
#'
#' \code{entity_evidence_annotation} has \code{entity_evidence_id}
#' (not \code{entity_id}) so the organism count joins through
#' \code{entity_evidence_resolution} (which carries
#' \code{entity_evidence_id} + \code{entity_id}) to filter for
#' entities of type \code{Organism:OM:0032}.
#'
#' @return Character.
#'
#' @keywords internal
#' @noRd
annotation_records_sql <- function() {
    paste(
        "SELECT COUNT(*)::bigint AS n_annotation_records",
        "  FROM entity_evidence_annotation",
        sep = "\n"
    )
}


#' Distinct organism entities with at least one annotation
#'
#' Starts from \code{entity_evidence_annotation} (the side with the
#' filter) and joins outward through
#' \code{entity_evidence_resolution} → \code{entity} → vocab so the
#' query plan rides the
#' \code{entity_evidence_resolution_entity_idx} +
#' \code{entity_type_taxonomy_idx} indexes instead of scanning all
#' 460K entities.
#'
#' @return Character.
#'
#' @keywords internal
#' @noRd
annotation_organisms_sql <- function() {
    paste(
        "SELECT COUNT(DISTINCT eer.entity_id)::bigint AS n_organisms",
        "  FROM entity_evidence_annotation eea",
        "  JOIN entity_evidence_resolution eer USING (entity_evidence_id)",
        "  JOIN entity e ON e.entity_id = eer.entity_id",
        "  JOIN vocab_entity_type vet USING (entity_type_id)",
        " WHERE vet.name = 'Organism:OM:0032'",
        sep = "\n"
    )
}


#' Resolve the slot-3 annotation metric per FR-043e fallback chain
#'
#' Walks \code{runtime.annotation.slot3_fallback_order} — typically
#' \code{c("diseases", "phenotypes", "tissues", "go_biological_process")} —
#' issuing the probe SQL for each candidate; the first slot whose
#' count meets \code{runtime.annotation.disease_min_rows} is chosen.
#' If no slot meets the threshold, the first slot in the list is
#' kept with whatever count it has and \code{state = "fallback"}.
#'
#' @param panel_id Character.
#' @param runtime Named list.
#' @param facet Character.
#'
#' @return Named list.
#'
#' @keywords internal
#' @noRd
resolve_annotation_slot3 <- function(panel_id, runtime, facet) {

    min_rows <- as.integer(runtime$annotation$disease_min_rows %||% 100L)
    order <- runtime$annotation$slot3_fallback_order
    if (length(order) == 0L) {
        order <- c(
            "diseases",
            "phenotypes",
            "tissues",
            "go_biological_process"
        )
    }
    default_branch <- order[[1L]]

    probes <- purrr::map(order, function(slot) {
        sql <- slot3_sql(slot)
        rows <- pg_query_panel(panel_id, sql, facet = facet)
        value <- as.integer(rows[[1L]] %||% 0L)
        list(
            slot      = slot,
            value     = value,
            sql       = sql,
            sql_hash  = sql_sha256(sql),
            deployment = attr(rows, "deployment"),
            result_hash = attr(rows, "result_hash"),
            row_count = nrow(rows)
        )
    })

    chosen <- NULL
    for (probe in probes) {
        if (probe$value >= min_rows) {
            chosen <- probe
            break
        }
    }
    fellback <- is.null(chosen)
    if (fellback) {
        chosen <- probes[[1L]]
    }

    metric_name <- slot3_metric_name(chosen$slot)

    list(
        metric_name      = metric_name,
        value            = chosen$value,
        definition_label = sprintf(
            "annotation_slot3:%s%s",
            chosen$slot,
            if (fellback) "+fallback_chain_exhausted" else ""
        ),
        state = if (
            chosen$slot == default_branch && !fellback
        ) {
            "populated"
        } else {
            "fallback"
        },
        deployment = chosen$deployment,
        sql_hash   = chosen$sql_hash,
        branch     = chosen$slot,
        query = list(
            metric_names = metric_name,
            sql          = chosen$sql,
            sql_hash     = chosen$sql_hash,
            deployment   = chosen$deployment,
            result_hash  = chosen$result_hash,
            row_count    = chosen$row_count,
            branch       = chosen$slot
        )
    )
}


#' Slot-3 metric name for a given slot category
#'
#' @param slot Character.
#'
#' @return Character.
#'
#' @keywords internal
#' @noRd
slot3_metric_name <- function(slot) {
    map <- c(
        diseases              = "n_diseases",
        phenotypes            = "n_phenotypes",
        tissues               = "n_tissues",
        go_biological_process = "n_go_biological_process"
    )
    map[[slot]] %||% sprintf("n_%s", slot)
}


#' SQL — probe count for a slot-3 candidate category
#'
#' @param slot Character.
#'
#' @return Character.
#'
#' @keywords internal
#' @noRd
slot3_sql <- function(slot) {

    if (slot == "diseases") {
        # MONDO is the unified disease ontology on dev5;
        # DOID is not loaded separately (most DOID terms have
        # MONDO equivalents). HP also carries disease-adjacent
        # phenotype terms but is reserved for the phenotypes slot.
        return(paste(
            "SELECT COUNT(DISTINCT term_id)::bigint AS n_diseases",
            "  FROM entity_ontology_term",
            " WHERE ontology_prefix = 'mondo'",
            sep = "\n"
        ))
    }
    if (slot == "phenotypes") {
        # The Human Phenotype Ontology uses prefix 'hp' on dev5
        # (not 'hpo').
        return(paste(
            "SELECT COUNT(DISTINCT term_id)::bigint AS n_phenotypes",
            "  FROM entity_ontology_term",
            " WHERE ontology_prefix = 'hp'",
            sep = "\n"
        ))
    }
    if (slot == "tissues") {
        # Uberon carries anatomical / tissue terms; the
        # Tissue:OM:0034 entity_type also exists but Uberon
        # ontology terms are the higher-coverage source.
        return(paste(
            "SELECT COUNT(DISTINCT term_id)::bigint AS n_tissues",
            "  FROM entity_ontology_term",
            " WHERE ontology_prefix = 'uberon'",
            sep = "\n"
        ))
    }
    if (slot == "go_biological_process") {
        # The schema does not separate GO sub-ontologies on dev5;
        # this counts all GO terms attached to entities and the
        # caption notes the approximation.
        return(paste(
            "SELECT COUNT(DISTINCT term_id)::bigint",
            "    AS n_go_biological_process",
            "  FROM entity_ontology_term",
            " WHERE ontology_prefix = 'go'",
            sep = "\n"
        ))
    }

    rlang::abort(sprintf(
        "FR-043e slot-3 fallback list contains unknown slot '%s'",
        slot
    ))
}


#' SQL — distinct subcellular localizations
#'
#' Uberon is the anatomy/tissue/compartment ontology on dev5; the
#' GO sub-ontology split (CC / BP / MF) is not present in the schema
#' so the digest uses Uberon as the localization proxy. The figure
#' caption MUST note this when the digest is the source of the
#' localization headline.
#'
#' @return Character.
#'
#' @keywords internal
#' @noRd
localization_sql <- function() {
    paste(
        "SELECT COUNT(DISTINCT term_id)::bigint AS n_localizations",
        "  FROM entity_ontology_term",
        " WHERE ontology_prefix = 'uberon'",
        sep = "\n"
    )
}
