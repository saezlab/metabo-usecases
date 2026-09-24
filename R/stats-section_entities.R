#' Section 1 — Entities (FR-043a)
#'
#' Reports unique entity counts for metabolites, foods, drugs, lipids,
#' xenobiotics, and proteins/genes (the union plus the gene-anchor
#' and unresolved-protein subtotals). Chemical-class membership reads
#' from \code{entity.chemical_class_id} joined to
#' \code{vocab_chemical_class}. The proteins/genes split reads from
#' \code{entity.entity_type_id} joined to \code{vocab_entity_type}.
#'
#' Routes to \code{dev5} via
#' \code{\link{pg_query_panel}} (post-2026-06-14 dev5 integrated-build
#' promotion).
#'
#' @param panel_id Character: panel id for deployment routing.
#'     Defaults to \code{"architecture"}.
#'
#' @return A tibble of DigestMetric rows with columns
#'     \code{section_id}, \code{metric_name}, \code{value},
#'     \code{definition_label}, \code{state}, \code{deployment},
#'     \code{sql_hash}. Also carries \code{attr(., "queries")} — a
#'     list of \code{{metric_name, sql, sql_hash, deployment,
#'     result_hash, row_count}} records the orchestrator transcribes
#'     into the digest sidecar.
#'
#' @importFrom logger log_info
#' @importFrom tibble tibble
#' @export
section_entities <- function(panel_id = "architecture") {

    facet <- "panel_a_stats_entities"

    chemical_sql <- paste(
        "SELECT vcc.name AS chemical_class,",
        "       COUNT(DISTINCT e.entity_id) AS n_entities",
        "  FROM entity e",
        "  JOIN vocab_chemical_class vcc USING (chemical_class_id)",
        " WHERE vcc.name IN ",
        "       ('metabolite','lipid','drug','food','xenobiotic')",
        " GROUP BY vcc.name",
        sep = "\n"
    )

    pog_sql <- paste(
        "SELECT",
        "    SUM(CASE WHEN vet.name = 'Gene:MI:0250'",
        "             THEN 1 ELSE 0 END)::bigint AS n_genes,",
        "    SUM(CASE WHEN vet.name = 'Protein:MI:0326'",
        "             THEN 1 ELSE 0 END)::bigint AS n_unresolved_proteins,",
        "    COUNT(*)::bigint AS n_protein_or_gene_union",
        "  FROM entity e",
        "  JOIN vocab_entity_type vet USING (entity_type_id)",
        " WHERE vet.name IN ('Gene:MI:0250','Protein:MI:0326')",
        sep = "\n"
    )

    chem_rows <- pg_query_panel(panel_id, chemical_sql, facet = facet)
    pog_rows  <- pg_query_panel(panel_id, pog_sql,      facet = facet)

    chem_lookup <- as.list(
        setNames(
            as.integer(chem_rows$n_entities),
            chem_rows$chemical_class
        )
    )

    chemical_classes <- c(
        "metabolite", "food", "drug", "lipid", "xenobiotic"
    )
    metric_names <- c(
        metabolite = "n_metabolites",
        food       = "n_foods",
        drug       = "n_drugs",
        lipid      = "n_lipids",
        xenobiotic = "n_xenobiotics"
    )

    chem_metrics <- purrr::map(chemical_classes, function(cls) {
        value <- chem_lookup[[cls]] %||% 0L
        tibble::tibble(
            section_id       = 1L,
            metric_name      = metric_names[[cls]],
            value            = as.integer(value),
            definition_label = sprintf("chemical_class:%s", cls),
            state            = if (value > 0L) "populated" else "empty",
            deployment       = attr(chem_rows, "deployment"),
            sql_hash         = sql_sha256(chemical_sql)
        )
    })

    n_genes <- as.integer(pog_rows$n_genes %||% 0L)
    n_unresolved <- as.integer(pog_rows$n_unresolved_proteins %||% 0L)
    n_union <- as.integer(pog_rows$n_protein_or_gene_union %||% 0L)

    pog_metrics <- list(
        tibble::tibble(
            section_id       = 1L,
            metric_name      = "n_proteins_or_genes",
            value            = n_union,
            definition_label = "entity_type:Gene:MI:0250+Protein:MI:0326",
            state            = "populated",
            deployment       = attr(pog_rows, "deployment"),
            sql_hash         = sql_sha256(pog_sql)
        ),
        tibble::tibble(
            section_id       = 1L,
            metric_name      = "n_genes",
            value            = n_genes,
            definition_label = "entity_type:Gene:MI:0250",
            state            = "populated",
            deployment       = attr(pog_rows, "deployment"),
            sql_hash         = sql_sha256(pog_sql)
        ),
        tibble::tibble(
            section_id       = 1L,
            metric_name      = "n_unresolved_proteins",
            value            = n_unresolved,
            definition_label = "entity_type:Protein:MI:0326",
            state            = "populated",
            deployment       = attr(pog_rows, "deployment"),
            sql_hash         = sql_sha256(pog_sql)
        )
    )

    metrics <- dplyr::bind_rows(c(chem_metrics, pog_metrics))

    queries <- list(
        chemical = list(
            metric_names = unname(metric_names),
            sql          = chemical_sql,
            sql_hash     = sql_sha256(chemical_sql),
            deployment   = attr(chem_rows, "deployment"),
            result_hash  = attr(chem_rows, "result_hash"),
            row_count    = nrow(chem_rows)
        ),
        protein_or_gene = list(
            metric_names = c(
                "n_proteins_or_genes",
                "n_genes",
                "n_unresolved_proteins"
            ),
            sql          = pog_sql,
            sql_hash     = sql_sha256(pog_sql),
            deployment   = attr(pog_rows, "deployment"),
            result_hash  = attr(pog_rows, "result_hash"),
            row_count    = nrow(pog_rows)
        )
    )

    attr(metrics, "queries") <- queries

    logger::log_info(
        "section_entities: ",
        "{nrow(metrics)} metrics emitted"
    )

    metrics
}
