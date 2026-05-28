#' The project-wide random seed
#'
#' A single fixed integer used by every pipeline script that invokes
#' stochastic operations (constitution Principle I and research.md
#' R-10). Using one constant keeps reproducibility auditable: any
#' deterministic difference between rebuilds against the same snapshot
#' is a bug.
#'
#' @return Integer scalar, \code{2026L}.
#'
#' @examples
#' set.seed(pipeline_seed())
#'
#' @export
pipeline_seed <- function() {
    2026L
}
