#' Canonicalised SHA-256 of a SQL string
#'
#' The orchestrator hashes each metric's verbatim SQL into a 64-hex
#' digest that becomes the metric's \code{sql_hash} field; the
#' \code{stats.provenance.json} sidecar carries the verbatim SQL
#' itself. The canonicalisation is just \code{trimws()}; minor
#' whitespace differences within a query are intentional and stay
#' in the hash (so a reformat is a deliberate change a reader can
#' detect).
#'
#' @param sql Character.
#'
#' @return 64-character hex string.
#'
#' @importFrom digest digest
#' @keywords internal
#' @noRd
sql_sha256 <- function(sql) {
    digest::digest(trimws(sql), algo = "sha256", serialize = FALSE)
}


#' Default digest-config path under the Figure 1 folder
#'
#' Returns the standard location of \code{digest-config.yaml} relative
#' to the repository root. Used by \code{digest_config()} when the
#' caller does not pass an explicit path.
#'
#' @return Character.
#'
#' @keywords internal
#' @noRd
default_digest_config_path <- function() {
    "figures/fig01-overview/panel-a-stats/digest-config.yaml"
}


#' Load the Panel A digest configuration
#'
#' Reads \code{figures/fig01-overview/panel-a-stats/digest-config.yaml}
#' (FR-043; see \code{contracts/panel-a-stats.schema.json}). The
#' returned list carries two top-level namespaces:
#'
#' - \code{definitions}: values that ARE mirrored into
#'   \code{stats.json -> metadata.definitions} (per the schema's
#'   \code{additionalProperties: false} contract over a fixed key
#'   set).
#' - \code{runtime}: implementation tunables (column names, fallback
#'   thresholds, …) that are NOT mirrored into the metadata block.
#'
#' The SHA-256 fingerprint of the on-disk YAML is also returned so the
#' digest sidecar can pin the config revision actually used.
#'
#' @param path Optional path to the YAML file. \code{NULL} (default)
#'     → \code{\link{default_digest_config_path}}.
#'
#' @return Named list with elements \code{definitions}, \code{runtime},
#'     \code{path} (absolute), \code{sha256}.
#'
#' @examples
#' \dontrun{
#' cfg <- digest_config()
#' cfg$definitions$transporter$resources    # [tcdb, slctables]
#' cfg$runtime$structures$inchikey_column   # "inchikey"
#' cfg$sha256                               # 64-hex sha-256 of the YAML
#' }
#'
#' @importFrom digest digest
#' @importFrom logger log_trace
#' @importFrom rlang abort
#' @importFrom yaml read_yaml
#' @export
digest_config <- function(path = NULL) {

    cfg_path <- path %||% default_digest_config_path()
    cfg_path <- normalizePath(cfg_path, mustWork = TRUE)

    cfg <- yaml::read_yaml(cfg_path)
    sha256 <- digest::digest(file = cfg_path, algo = "sha256")

    if (!is.list(cfg$definitions)) {
        rlang::abort(sprintf(
            "digest-config.yaml at %s is missing the 'definitions:' block",
            cfg_path
        ))
    }
    if (!is.list(cfg$runtime)) {
        rlang::abort(sprintf(
            "digest-config.yaml at %s is missing the 'runtime:' block",
            cfg_path
        ))
    }

    validate_definitions_block(cfg$definitions, cfg_path)
    validate_runtime_block(cfg$runtime, cfg_path)

    logger::log_trace(
        "digest_config loaded from {cfg_path} (sha256={substr(sha256, 1L, 12L)}…)"
    )

    list(
        definitions = cfg$definitions,
        runtime     = cfg$runtime,
        path        = cfg_path,
        sha256      = sha256
    )
}


#' Convenience accessor for the \code{definitions:} namespace
#'
#' @param config Optional pre-loaded config (from
#'     \code{\link{digest_config}}); \code{NULL} → load fresh.
#'
#' @return Named list (the \code{definitions} block).
#'
#' @export
digest_definitions <- function(config = NULL) {
    cfg <- config %||% digest_config()
    cfg$definitions
}


#' Convenience accessor for the \code{runtime:} namespace
#'
#' @param config Optional pre-loaded config (from
#'     \code{\link{digest_config}}); \code{NULL} → load fresh.
#'
#' @return Named list (the \code{runtime} block).
#'
#' @export
digest_runtime <- function(config = NULL) {
    cfg <- config %||% digest_config()
    cfg$runtime
}


#' Build the \code{stats.json -> metadata.definitions} block
#'
#' Combines the YAML's \code{definitions:} defaults with the
#' chosen-at-runtime branches (which pathway / reaction definition
#' fired, which annotation slot 3 carries) so the returned structure
#' validates against
#' \code{contracts/panel-a-stats.schema.json} § \code{metadata.definitions}.
#'
#' Only the chosen branch strings end up in the result; the
#' pathway / reaction full \code{{preferred, fallback}} objects from
#' the YAML are NOT mirrored.
#'
#' @param definitions The \code{definitions:} namespace from
#'     \code{\link{digest_config}}.
#' @param branches Named list with one or more of:
#'     \code{pathway} (\code{"Pathway:OM:0014"} or \code{"annotation"}),
#'     \code{reaction} (\code{"Reaction:OM:0015"} or
#'     \code{"predicate-classification"}),
#'     \code{annotation_slot3} (\code{"diseases"} | \code{"phenotypes"}
#'     | \code{"tissues"} | \code{"go_biological_process"}). Missing
#'     keys fall back to the YAML default for that key (\code{preferred}
#'     for pathway / reaction; the YAML's \code{annotation_slot3}
#'     literal).
#'
#' @return Named list shaped to the schema's
#'     \code{metadata.definitions} block.
#'
#' @export
mirror_to_metadata <- function(definitions, branches = list()) {

    pathway_choice <-
        branches$pathway %||% definitions$pathway$preferred
    reaction_choice <-
        branches$reaction %||% definitions$reaction$preferred
    slot3_choice <-
        branches$annotation_slot3 %||% definitions$annotation_slot3

    list(
        transporter      = definitions$transporter,
        receptor         = definitions$receptor,
        pathway          = pathway_choice,
        reaction         = reaction_choice,
        annotation_slot3 = slot3_choice
    )
}


#' Validate the YAML's \code{definitions:} block has every required key
#'
#' Throws a clear error if a required key is missing rather than
#' letting a downstream lookup return \code{NULL}. The validator is
#' aligned with \code{contracts/panel-a-stats.schema.json}.
#'
#' @param defs Named list (the \code{definitions:} block).
#' @param path Character: file path (for the error message).
#'
#' @return Invisibly \code{NULL}; called for side effect of throwing.
#'
#' @importFrom rlang abort
#' @keywords internal
#' @noRd
validate_definitions_block <- function(defs, path) {

    required <- list(
        list(
            keys  = c("transporter", "resources"),
            label = "definitions.transporter.resources"
        ),
        list(
            keys  = c("transporter", "uniprot_keywords"),
            label = "definitions.transporter.uniprot_keywords"
        ),
        list(
            keys  = c("receptor", "resources"),
            label = "definitions.receptor.resources"
        ),
        list(
            keys  = c("receptor", "uniprot_keywords"),
            label = "definitions.receptor.uniprot_keywords"
        ),
        list(
            keys  = c("pathway", "preferred"),
            label = "definitions.pathway.preferred"
        ),
        list(
            keys  = c("reaction", "preferred"),
            label = "definitions.reaction.preferred"
        ),
        list(
            keys  = "annotation_slot3",
            label = "definitions.annotation_slot3"
        )
    )

    purrr::walk(required, function(entry) {
        if (is.null(purrr::pluck(defs, !!!entry$keys))) {
            rlang::abort(sprintf(
                "digest-config.yaml at %s is missing required key '%s'",
                path,
                entry$label
            ))
        }
    })

    invisible(NULL)
}


#' Validate the YAML's \code{runtime:} block has every required key
#'
#' @param rt Named list (the \code{runtime:} block).
#' @param path Character: file path (for the error message).
#'
#' @return Invisibly \code{NULL}.
#'
#' @importFrom rlang abort
#' @keywords internal
#' @noRd
validate_runtime_block <- function(rt, path) {

    required <- list(
        list(
            keys  = c("structures", "inchikey_column"),
            label = "runtime.structures.inchikey_column"
        ),
        list(
            keys  = c("annotation", "disease_min_rows"),
            label = "runtime.annotation.disease_min_rows"
        ),
        list(
            keys  = c("annotation", "slot3_fallback_order"),
            label = "runtime.annotation.slot3_fallback_order"
        )
    )

    purrr::walk(required, function(entry) {
        if (is.null(purrr::pluck(rt, !!!entry$keys))) {
            rlang::abort(sprintf(
                "digest-config.yaml at %s is missing required key '%s'",
                path,
                entry$label
            ))
        }
    })

    invisible(NULL)
}
