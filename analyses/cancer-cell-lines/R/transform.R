# Transformation of the Shorthouse intensities (decision: audit/01_scale_transform.Rmd).

#' Half-minimum imputation of non-positive values, then log2
#'
#' Values <= 0 cannot be log-transformed. They are rare and occur in
#' low-abundance metabolites, so they are treated as below the detection limit:
#' set to NA and replaced by half of the metabolite's minimum positive value
#' (Wei et al. 2018) with MetaProViz::mvi_imputation(), then log2-transformed.
#'
#' @param x Matrix of intensities, samples (rows) × metabolites (columns).
#' @return Matrix of the same shape on log2 scale.
log2_half_minimum <- function(x) {
    d <- as.data.frame(x, check.names = FALSE)
    d[d <= 0] <- NA

    # mvi_imputation() imputes within each condition; one condition for all
    # samples gives half of the minimum over all samples of a metabolite
    meta <- data.frame(row.names = rownames(d), condition = rep("all", nrow(d)))
    imputed <- suppressMessages(MetaProViz::mvi_imputation(
        data = d, metadata_sample = meta, metadata_info = c(Conditions = "condition")
    ))

    out <- as.matrix(imputed[rownames(x), colnames(x)])
    stopifnot(!anyNA(out), all(out > 0))
    log2(out)
}
