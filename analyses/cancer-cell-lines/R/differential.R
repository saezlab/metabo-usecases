# Lines, design and limma fit of one mutation contrast (decisions D18, D19, D22, D25).

#' Altered status of one gene or pathway group per line
#'
#' Oncogenes count with a likely-driver variant, tumour suppressors with any
#' protein-changing variant or a homozygous deletion (D19), genes defined by
#' deletion alone with a homozygous deletion (D25), other genes with a
#' likely-driver variant or a homozygous deletion (definition A of
#' analysis/03). A pathway group is altered if any member is.
#'
#' @param data The analysis data set (`shorthouse_analysis_v1.rds`).
#' @param lines Line names.
#' @param name Gene symbol or name of a group in `config$contrast$pathway_groups`.
#' @param config Configuration from `load_config()`.
#' @return Logical vector named by line; NA where the status is unknown
#'   (genes defined by deletion, lines without copy number).
altered_status <- function(data, lines, name, config) {
    params <- config$contrast
    genes <- params$pathway_groups[[name]] %||% name
    status <- alteration_status(data, lines, genes) |>
        dplyr::mutate(altered = dplyr::case_when(
            gene %in% params$deletion_genes & !has_copy_number ~ NA,
            gene %in% params$deletion_genes ~ deleted,
            gene %in% params$oncogenes ~ likely_driver,
            gene %in% params$tumour_suppressors ~ any_mutation | deleted,
            TRUE ~ likely_driver | deleted
        )) |>
        dplyr::group_by(line) |>
        dplyr::summarise(altered = if (any(altered %in% TRUE)) TRUE else if (anyNA(altered)) NA else FALSE,
                         .groups = "drop")
    stats::setNames(status$altered, status$line)[lines]
}

#' Lines and covariates of one contrast
#'
#' Starts from the lines used in models, restricts them to the tissues of the
#' contrast and drops lines with unknown group status. With
#' `tissues: contributing`, only tissues with at least
#' `min_per_group_per_tissue` altered and unaltered lines are kept (D20).
#'
#' @param data The analysis data set.
#' @param contrast One entry of `config$differential$contrasts`.
#' @param config Configuration from `load_config()`.
#' @return Data frame with line, tissue, batch, `group` and one logical column
#'   per covariate.
contrast_lines <- function(data, contrast, config) {
    lines <- dplyr::filter(data$lines, in_models)
    if (!identical(contrast$tissues, "contributing")) lines <- dplyr::filter(lines, tissue %in% contrast$tissues)
    d <- dplyr::transmute(lines, line, tissue, batch)
    d$group <- altered_status(data, d$line, contrast$group, config)
    for (cov in contrast$covariates) d[[cov]] <- altered_status(data, d$line, cov, config)
    d <- dplyr::filter(d, !is.na(group))

    if (identical(contrast$tissues, "contributing")) {
        n_min <- config$contrast$min_per_group_per_tissue
        keep <- d |>
            dplyr::group_by(tissue) |>
            dplyr::summarise(ok = sum(group) >= n_min & sum(!group) >= n_min, .groups = "drop") |>
            dplyr::filter(ok) |>
            dplyr::pull(tissue)
        d <- dplyr::filter(d, tissue %in% keep)
    }
    d
}

#' Design matrix of one contrast
#'
#' `~ batch + tissue + covariates + group`; batch, tissue and covariates enter
#' only when they vary among the lines, so the design stays estimable.
#'
#' @param d Output of `contrast_lines()`.
#' @param covariates Covariate names.
#' @return Design matrix; the last column is the group effect (`groupTRUE`).
contrast_design <- function(d, covariates) {
    varies <- function(x) dplyr::n_distinct(x) > 1
    terms <- c(
        if (varies(d$batch)) "batch",
        if (varies(d$tissue)) "tissue",
        Filter(function(cov) varies(d[[cov]]), covariates),
        "group"
    )
    X <- stats::model.matrix(stats::reformulate(terms), data = d)
    if (qr(X)$rank < ncol(X)) stop("Design is not of full rank: ", paste(terms, collapse = " + "))
    X
}

#' limma fit of one contrast on all ions
#'
#' @param data The analysis data set (log2 values, lines x ions).
#' @param d Output of `contrast_lines()`.
#' @param X Output of `contrast_design()`.
#' @return List with the `fit` (eBayes) and the full result `table` for the
#'   group coefficient, joined to the ion annotation.
fit_contrast <- function(data, d, X) {
    # limma expects features in rows
    y <- t(data$log2[d$line, , drop = FALSE])
    fit <- limma::eBayes(limma::lmFit(y, X))
    table <- limma::topTable(fit, coef = ncol(X), number = Inf, sort.by = "none") |>
        tibble::rownames_to_column("feature_id") |>
        dplyr::left_join(
            dplyr::select(data$ions, feature_id, name, hmdb, kegg, has_imputed_values, replicate_signal_ratio),
            by = "feature_id"
        )
    list(fit = fit, table = table, y = y)
}
