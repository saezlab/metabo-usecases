#' Load the old COSMOS prior-knowledge network from cosmosR
#'
#' Installs \pkg{cosmosR} from \code{saezlab/cosmosR} on GitHub if not
#' already present, then loads the \code{meta_network} dataset and
#' categorises edges by node-prefix pattern (FR-011, FR-011a):
#'
#' \describe{
#'   \item{PPI}{Both source and target are plain gene symbols — mixture
#'     of OmniPath-LR and Collectri; indistinguishable without metadata.}
#'   \item{metabolite_protein}{Source starts with \code{Metab__}, target
#'     is a plain gene symbol (STITCH and Recon3D bidirectional).}
#'   \item{metabolic_reaction}{Plain gene symbol source, target starts
#'     with \code{Metab__} (Recon3D GEM reactions).}
#' }
#'
#' Edges where either node carries a \code{Gene} prefix are bridge nodes
#' (Recon3D internal connectors) and are removed before counting.
#'
#' @return A tibble with columns \code{interaction_type} and
#'     \code{n_edges}, ordered by \code{n_edges} descending.
#'     Carries a \code{"source_pkg"} attribute (package name + version)
#'     and a \code{"fingerprint"} attribute (package version string) for
#'     sidecar recording.
#'
#' @importFrom rlang abort
#' @importFrom dplyr mutate filter count rename arrange desc case_when
#' @importFrom tibble as_tibble
#' @export
cosmos_old_pkn <- function() {

    if (!requireNamespace("devtools", quietly = TRUE)) {
        install.packages("devtools")
    }
    if (!requireNamespace("cosmosR", quietly = TRUE)) {
        devtools::install_github("saezlab/cosmosR")
    }
    library(cosmosR)
    data("meta_network")

    if (!exists("meta_network")) {
        rlang::abort(
            "cosmosR::meta_network dataset not found after loading the package."
        )
    }

    pkn <- tibble::as_tibble(meta_network)

    # Categorise by node-prefix pattern.
    # Gene__ prefix = Recon3D gene node (not a bridge per se).
    # Bridges are specifically plain-gene → Gene__ rows (a gene symbol
    # connected to its own Recon3D duplicate); those are removed.
    pkn <- dplyr::mutate(
        pkn,
        is_metab_src = grepl("^Metab__", source),
        is_metab_tgt = grepl("^Metab__", target),
        is_gene_src  = grepl("^Gene",    source),
        is_gene_tgt  = grepl("^Gene",    target),
        is_plain_src = !is_metab_src & !is_gene_src,
        is_plain_tgt = !is_metab_tgt & !is_gene_tgt,
        interaction_type = dplyr::case_when(
            is_plain_src & is_gene_tgt    ~ NA_character_,  # bridge: remove
            is_plain_src & is_plain_tgt   ~ "PPI",
            is_metab_src & is_plain_tgt   ~ "metabolite_protein",
            is_metab_src & is_gene_tgt    ~ "metabolic_reaction",
            is_gene_src  & is_metab_tgt   ~ "metabolic_reaction",
            .default = NA_character_
        )
    )

    pkn <- dplyr::filter(pkn, !is.na(interaction_type))

    counts <- pkn |>
        dplyr::count(interaction_type) |>
        dplyr::rename(n_edges = n) |>
        dplyr::arrange(dplyr::desc(n_edges))

    pkg_version <- as.character(utils::packageVersion("cosmosR"))
    attr(counts, "source_pkg")  <- paste0("cosmosR@", pkg_version)
    attr(counts, "fingerprint") <- pkg_version

    counts
}
