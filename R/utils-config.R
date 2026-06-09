#' Known OmniPath deployments on the beauty workstation
#'
#' Static table of the five deployments documented in
#' \code{saezverse/human/dev-deployments.md}. Used by
#' \code{\link{load_connection}} to resolve a deployment name into
#' connection parameters. Per the 2026-06-09 cycle-001 + 002 handover,
#' the figure pipeline routes panels to \code{dev3} (gene-centric +
#' labels, default) and \code{dev4} (protein-centric + RDKit
#' structural specificity + RaMP-conflict tables); \code{dev5} is
#' reserved for the ongoing integrated build and is refused; \code{prod}
#' is excluded from the active rotation and requires an explicit
#' opt-in.
#'
#' @return A tibble with columns \code{name}, \code{url}, \code{db_port},
#'     \code{db_name}, \code{role}. \code{role} classifies the
#'     deployment as \code{"figures-default"},
#'     \code{"figures-structures"}, \code{"reserved"}, or
#'     \code{"opt-in"}.
#'
#' @importFrom tibble tribble
#' @keywords internal
#' @noRd
known_deployments <- function() {
    tibble::tribble(
        ~name,  ~url,                  ~db_port, ~db_name,   ~role,
        "prod", "dev.omnipathdb.org",  5485L,    "omnipath", "opt-in",
        "dev2", "dev2.omnipathdb.org", 5402L,    "omnipath", "opt-in",
        "dev3", "dev3.omnipathdb.org", 5403L,    "omnipath", "figures-default",
        "dev4", "dev4.omnipathdb.org", 5404L,    "omnipath", "figures-structures",
        "dev5", "dev5.omnipathdb.org", 5405L,    "omnipath", "reserved"
    )
}


#' Read and validate the pipeline's deployment connection config
#'
#' Reads \code{~/.config/metabo-figures/connection.yaml} (or a
#' project-local \code{connection.yaml}, or env-only) and resolves a
#' deployment name into a connection record. The 2026-06-09 cycle-001
#' + 002 handover put the figure pipeline on a per-panel deployment
#' matrix: \code{dev3} (gene-centric + labels) is the default;
#' \code{dev4} (protein-centric + RDKit structural specificity +
#' RaMP-conflict tables) is the required override for the panels
#' driven by FR-007a Structures facet, FR-007e, FR-007f, FR-015.
#' \code{dev5} is reserved for the ongoing integrated build and is
#' refused. \code{prod} / \code{dev2} are excluded from the active
#' rotation and require \code{allow_optin = TRUE} to load.
#'
#' Per-deployment port overrides are read from the env vars
#' \code{PGPORT_DEV3}, \code{PGPORT_DEV4}, \code{PGPORT_DEV5},
#' \code{PGPORT_PROD}, \code{PGPORT_DEV2} — useful when the pipeline
#' runs against tunnelled deployments. \code{PGPORT} (the libpq
#' default) is honoured only when no per-deployment override is set.
#'
#' Credential fields are NOT returned in the deployment record — they
#' are resolved at query time from the environment, a credentials
#' file, or \code{~/.pgpass}, depending on \code{credentials_source}.
#'
#' @param deployment Character: explicit deployment name. \code{NULL}
#'     (default) → read \code{default_deployment} (or legacy
#'     \code{deployment}) from \code{connection.yaml}.
#' @param config_path Character: explicit path to a connection.yaml
#'     file. \code{NULL} (default) → check project-local
#'     \code{connection.yaml}, then user config.
#' @param allow_optin Logical: when \code{TRUE}, permit loading an
#'     opt-in deployment (\code{prod} / \code{dev2}). Default
#'     \code{FALSE}: opt-in deployments fail-fast with a clear
#'     diagnostic.
#'
#' @return A list with elements \code{name}, \code{url},
#'     \code{db_host}, \code{db_port}, \code{db_name},
#'     \code{role}, \code{credentials_source},
#'     \code{credentials_path}.
#'
#' @examples
#' \dontrun{
#' conn <- load_connection()           # → dev3 by default
#' conn <- load_connection("dev4")     # → structural-specificity build
#' }
#'
#' @importFrom yaml read_yaml
#' @importFrom logger log_info
#' @importFrom dplyr filter
#' @importFrom rlang abort
#' @export
load_connection <- function(
    deployment = NULL,
    config_path = NULL,
    allow_optin = FALSE
) {

    # NSE vs. R CMD check workaround
    name <- NULL

    cfg <- resolve_connection_config(config_path)
    requested <- deployment %||% cfg$default_deployment

    record <-
        known_deployments() %>%
        dplyr::filter(name == requested)

    if (nrow(record) == 0L) {
        rlang::abort(sprintf(
            "Unknown deployment '%s'. Known: %s",
            requested,
            paste(known_deployments()$name, collapse = ", ")
        ))
    }

    if (identical(record$role, "reserved")) {
        rlang::abort(sprintf(
            paste0(
                "Deployment '%s' is reserved for the ongoing integrated ",
                "build and MUST NOT be used by the figure pipeline. ",
                "Use dev3 (default) or dev4 (structural specificity / ",
                "RaMP) instead."
            ),
            requested
        ))
    }

    if (identical(record$role, "opt-in") && !isTRUE(allow_optin)) {
        rlang::abort(sprintf(
            paste0(
                "Deployment '%s' is excluded from the active figure ",
                "rotation (cycles 001 + 002 promoted dev3 + dev4 as ",
                "the targets). Pass allow_optin = TRUE to override."
            ),
            requested
        ))
    }

    list(
        name                = record$name,
        url                 = record$url,
        db_host             = Sys.getenv("PGHOST", unset = "localhost"),
        db_port             = resolve_port(record$name, record$db_port),
        db_name             = record$db_name,
        role                = record$role,
        credentials_source  = cfg$credentials_source,
        credentials_path    = cfg$credentials_path
    )
}


#' Apply per-deployment PGPORT env override, falling back to the table
#'
#' @param name Character: deployment name (e.g. \code{"dev3"}).
#' @param default_port Integer: port from \code{known_deployments()}.
#'
#' @return Integer port.
#'
#' @keywords internal
#' @noRd
resolve_port <- function(name, default_port) {

    env_key <- paste0("PGPORT_", toupper(name))
    override <- Sys.getenv(env_key, unset = "")

    if (nzchar(override)) {
        as.integer(override)
    } else {
        as.integer(default_port)
    }
}


#' Default null-coalescing helper used by \code{load_connection}
#'
#' Kept tiny and local — magrittr's pipe is already imported, but
#' \code{\%||\%} is not exported from anywhere convenient.
#'
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x


#' Locate and parse the connection.yaml the user has configured
#'
#' Supports the cycle-001+002 schema (\code{default_deployment} +
#' \code{overrides}) and the legacy single-deployment shape
#' (\code{deployment}). The legacy field is silently promoted to
#' \code{default_deployment} so older configs keep working.
#'
#' @param config_path Optional explicit override path.
#'
#' @return Named list with \code{default_deployment},
#'     \code{credentials_source}, \code{credentials_path}, and
#'     optionally \code{overrides} (a nested list keyed by panel id).
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

    # Legacy single-deployment compatibility.
    if (is.null(cfg$default_deployment) && !is.null(cfg$deployment)) {
        cfg$default_deployment <- cfg$deployment
        cfg$deployment <- NULL
    }
    if (is.null(cfg$default_deployment)) {
        rlang::abort(
            "connection.yaml is missing 'default_deployment'"
        )
    }
    if (is.null(cfg$credentials_source)) {
        cfg$credentials_source <- "env"
    }
    if (is.null(cfg$overrides)) {
        cfg$overrides <- list()
    }

    cfg
}


#' Quick sanity-check of the figure-pipeline deployments
#'
#' Resolves the per-deployment connection records the figure pipeline
#' will touch — by default \code{dev3} (gene-centric + labels, the
#' default per FR-030's matrix) and \code{dev4} (protein-centric
#' build, required for FR-007a Structures facet, FR-007e, FR-007f,
#' FR-015) — and logs a one-line summary per deployment. The full
#' \code{build_manifest} probe lives in
#' \code{\link{build_manifest_for}}; this helper does not open a
#' Postgres handle.
#'
#' @param deployments Character vector of deployment names to probe.
#'     Defaults to the figure-pipeline pair \code{c("dev3", "dev4")}.
#'
#' @return Invisibly a list keyed by deployment name with the
#'     resolved connection record.
#'
#' @importFrom logger log_info log_warn
#' @importFrom purrr map
#' @export
check_deployment <- function(deployments = c("dev3", "dev4")) {

    records <- purrr::map(deployments, function(name) {
        record <- load_connection(deployment = name)
        logger::log_info(
            "Deployment OK: {record$name} at ",
            "{record$db_host}:{record$db_port} (role={record$role})"
        )
        record
    })

    names(records) <- deployments
    invisible(records)
}
