#' Known OmniPath deployments on the beauty workstation
#'
#' Static table of the five deployments documented in
#' \code{saezverse/human/dev-deployments.md}. Used by
#' \code{\link{load_connection}} to resolve a deployment name into
#' connection parameters.
#'
#' @return A tibble with columns \code{name}, \code{url}, \code{db_port},
#'     \code{db_name}, \code{usable}.
#'
#' @importFrom tibble tribble
#' @keywords internal
#' @noRd
known_deployments <- function() {
    tibble::tribble(
        ~name,    ~url,                       ~db_port, ~db_name,    ~usable,
        "prod",   "dev.omnipathdb.org",       5485L,    "omnipath",  TRUE,
        "dev2",   "dev2.omnipathdb.org",      5402L,    "omnipath",  TRUE,
        "dev3",   "dev3.omnipathdb.org",      5403L,    "omnipath",  TRUE,
        "dev4",   "dev4.omnipathdb.org",      5404L,    "omnipath",  TRUE,
        "dev5",   "dev5.omnipathdb.org",      5405L,    "omnipath",  FALSE
    )
}


#' Read and validate the pipeline's deployment connection config
#'
#' Reads \code{~/.config/metabo-figures/connection.yaml} (or a
#' project-local \code{connection.yaml}, or env-only) and resolves it
#' into a deployment record. Refuses \code{dev5} as long as it is
#' marked unpopulated in \code{\link{known_deployments}} (spec edge
#' case: unreachable deployment).
#'
#' Credential fields are NOT returned in the deployment record — they
#' are resolved at query time from the environment, a credentials
#' file, or \code{~/.pgpass}, depending on \code{credentials_source}.
#'
#' @param config_path Character: explicit path to a connection.yaml
#'     file. \code{NULL} (default) → check project-local
#'     \code{connection.yaml}, then user config.
#'
#' @return A list with elements \code{name}, \code{url},
#'     \code{db_host}, \code{db_port}, \code{db_name},
#'     \code{credentials_source}, \code{credentials_path}.
#'
#' @examples
#' \dontrun{
#' conn <- load_connection()
#' conn$name        # e.g. "dev3"
#' conn$db_port     # e.g. 5403
#' }
#'
#' @importFrom yaml read_yaml
#' @importFrom logger log_info log_warn
#' @importFrom dplyr filter
#' @importFrom rlang abort
#' @export
load_connection <- function(config_path = NULL) {

    # NSE vs. R CMD check workaround
    name <- usable <- NULL

    cfg <- resolve_connection_config(config_path)

    deployment <-
        known_deployments() %>%
        dplyr::filter(name == cfg$deployment)

    if (nrow(deployment) == 0L) {
        rlang::abort(sprintf(
            "Unknown deployment '%s'. Known: %s",
            cfg$deployment,
            paste(known_deployments()$name, collapse = ", ")
        ))
    }

    if (isFALSE(deployment$usable)) {
        rlang::abort(sprintf(
            paste0(
                "Deployment '%s' is currently unusable per ",
                "dev-deployments.md (e.g. dev5 is empty). ",
                "Pick a different deployment in connection.yaml."
            ),
            cfg$deployment
        ))
    }

    logger::log_info(
        "Resolved deployment {cfg$deployment} (port {deployment$db_port})"
    )

    list(
        name                = deployment$name,
        url                 = deployment$url,
        db_host             = Sys.getenv("PGHOST", unset = "localhost"),
        db_port             = deployment$db_port,
        db_name             = deployment$db_name,
        credentials_source  = cfg$credentials_source,
        credentials_path    = cfg$credentials_path
    )
}


#' Locate and parse the connection.yaml the user has configured
#'
#' @param config_path Optional explicit override path.
#' @return Named list with \code{deployment}, \code{credentials_source},
#'     \code{credentials_path}.
#'
#' @importFrom yaml read_yaml
#' @importFrom rlang abort
#' @keywords internal
#' @noRd
resolve_connection_config <- function(config_path = NULL) {

    candidates <- if (!is.null(config_path)) {
        config_path
    } else {
        c(
            "connection.yaml",
            file.path("~", ".config", "metabo-figures", "connection.yaml")
        )
    }
    candidates <- path.expand(candidates)
    found <- candidates[file.exists(candidates)]

    if (length(found) == 0L) {
        rlang::abort(paste0(
            "No connection.yaml found. ",
            "Copy inst/extdata/connection.template.yaml ",
            "to ~/.config/metabo-figures/connection.yaml."
        ))
    }

    cfg <- yaml::read_yaml(found[[1L]])

    if (is.null(cfg$deployment)) {
        rlang::abort("connection.yaml is missing 'deployment'")
    }
    if (is.null(cfg$credentials_source)) {
        cfg$credentials_source <- "env"
    }

    cfg
}


#' Quick sanity-check of the configured deployment
#'
#' Loads the connection, opens a Postgres handle, prints the resolved
#' deployment and (if implemented) the per-build snapshot identifiers
#' and partial-build flags. Intended for use during configuration of a
#' fresh checkout (see CONFIGURATION.md and quickstart.md).
#'
#' Phase 2A note: the Postgres connection step and snapshot derivation
#' land in Phases 2C (T025–T026); this function logs the resolved
#' deployment and exits successfully without attempting the query.
#'
#' @return Invisibly the resolved deployment list.
#'
#' @importFrom logger log_info
#' @export
check_deployment <- function() {

    conn <- load_connection()
    logger::log_info(
        "Deployment OK: {conn$name} at {conn$db_host}:{conn$db_port}"
    )

    invisible(conn)
}
