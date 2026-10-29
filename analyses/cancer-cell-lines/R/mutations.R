# Mutation and deletion status per line and gene (DepMap 24Q4, decision D10).

#' Alteration status for a set of lines and genes
#'
#' @param data The analysis data set (`shorthouse_analysis_v1.rds`).
#' @param lines Line names to include.
#' @param genes Gene symbols to include.
#' @return One row per line and gene with
#'   `any_mutation` (any protein-changing variant),
#'   `likely_driver` (variant flagged by DepMap as hotspot, likely loss of
#'   function, Hess driver or high-impact oncogene/tumour-suppressor variant),
#'   `deleted` (homozygous deletion; FALSE for lines without copy number),
#'   `protein_changes`, and the line's `tissue` and `batch`.
alteration_status <- function(data, lines, genes) {
    mutations <- data$mutation_status |>
        dplyr::filter(line %in% lines, gene %in% genes) |>
        dplyr::transmute(line, gene, any_mutation = n_variants > 0, likely_driver, protein_changes)
    deletions <- data$homozygous_deletions |>
        dplyr::filter(line %in% lines, gene %in% genes) |>
        dplyr::mutate(deleted = TRUE)

    tidyr::expand_grid(line = lines, gene = genes) |>
        dplyr::left_join(mutations, by = c("line", "gene")) |>
        dplyr::left_join(deletions, by = c("line", "gene")) |>
        dplyr::mutate(
            any_mutation = dplyr::coalesce(any_mutation, FALSE),
            likely_driver = dplyr::coalesce(likely_driver, FALSE),
            deleted = dplyr::coalesce(deleted, FALSE)
        ) |>
        dplyr::left_join(dplyr::select(data$lines, line, tissue, batch), by = "line")
}
