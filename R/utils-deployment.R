#' Per-panel deployment resolver
#'
#' Maps a (panel, facet) pair to the OmniPath deployment label that
#' SHOULD serve its data. Encodes the FR-030 deployment matrix from
#' the 2026-06-09 cycle-001 + 002 handover: \code{dev3} is the
#' default; \code{dev4} is required for FR-007a Structures facet,
#' FR-007e (structural specificity × chemical category), FR-007f
#' (RaMP comparison), and FR-015 (Methods table — RaMP comparison).
#'
#' User overrides from \code{connection.yaml::overrides} win over the
#' built-in registry — useful when a fresh integrated build lands on
#' \code{dev5} or a future port and the team wants to collapse the
#' per-panel split without a code change.
#'
#' @param panel_id Character: panel identifier
#'     (e.g. \code{"fig01-overview"} or
#'     \code{"tab02-ramp-comparison"}).
#' @param facet Character or \code{NULL}: facet identifier
#'     (e.g. \code{"structures"}, \code{"ramp_conflict"},
#'     \code{"panel_e"}). Optional; the registry first looks for an
#'     exact (panel, facet) match, then falls back to the panel-level
#'     default, then to the global default.
#' @param config Optional connection config to consult for user
#'     overrides. \code{NULL} (default) → load via
#'     \code{\link{resolve_connection_config}} from the standard
#'     candidates.
#'
#' @return Character: deployment label
#'     (e.g. \code{"dev3"}, \code{"dev4"}).
#'
#' @examples
#' \dontrun{
#' panel_deployment("fig01-overview")                          # "dev3"
#' panel_deployment("fig01-overview", "structures")            # "dev4"
#' panel_deployment("fig01-overview", "ramp_conflict")         # "dev4"
#' panel_deployment("tab02-ramp-comparison")                   # "dev4"
#' }
#'
#' @importFrom logger log_trace
#' @export
panel_deployment <- function(panel_id, facet = NULL, config = NULL) {

    cfg <- config %||% resolve_connection_config()
    registry <- deployment_registry()

    panel_overrides <- cfg$overrides[[panel_id]]
    panel_registry  <- registry[[panel_id]]

    chosen <- NULL

    # 1. User-config (panel, facet) override.
    if (!is.null(facet) && !is.null(panel_overrides[[facet]])) {
        chosen <- panel_overrides[[facet]]
    }

    # 2. Built-in (panel, facet) entry.
    if (is.null(chosen) && !is.null(facet) &&
        !is.null(panel_registry[[facet]])) {
        chosen <- panel_registry[[facet]]
    }

    # 3. User-config panel-level default.
    if (is.null(chosen) && !is.null(panel_overrides$default)) {
        chosen <- panel_overrides$default
    }

    # 4. Built-in panel-level default.
    if (is.null(chosen) && !is.null(panel_registry$default)) {
        chosen <- panel_registry$default
    }

    # 5. Global default from connection.yaml.
    if (is.null(chosen)) {
        chosen <- cfg$default_deployment
    }

    logger::log_trace(
        "panel_deployment(panel_id='{panel_id}', facet='{facet %||% \"\"}') ",
        "→ {chosen}"
    )

    chosen
}


#' Built-in (panel, facet) → deployment registry
#'
#' Encodes the FR-030 deployment matrix as a nested named list. Each
#' panel maps to a list whose \code{default} (optional) is the
#' panel-level default and whose other keys are facet-level overrides.
#'
#' @return Named list of named character vectors.
#'
#' @keywords internal
#' @noRd
deployment_registry <- function() {
    list(
        `fig01-overview` = list(
            structures     = "dev4",
            panel_e        = "dev4",
            ramp_conflict  = "dev4"
        ),
        `tab02-ramp-comparison` = list(
            default = "dev4"
        ),
        `fig03-metalinks-versions` = list(
            default = "dev4"
        )
    )
}


#' Open and cache a Postgres handle per deployment label
#'
#' Caches connections inside the package environment so the same
#' rebuild reuses the handle across many panels. Callers SHOULD NOT
#' close these handles directly; \code{\link{pg_close_all_panel}}
#' tears the pool down at the end of the rebuild.
#'
#' @param deployment Character: deployment label
#'     (\code{"dev3"} | \code{"dev4"} | …).
#'
#' @return A live DBI handle.
#'
#' @importFrom DBI dbConnect dbIsValid
#' @importFrom RPostgres Postgres
#' @importFrom logger log_info
#' @export
pg_connect_panel <- function(deployment) {

    pool <- panel_connection_pool()
    handle <- pool[[deployment]]

    if (is.null(handle) || isFALSE(DBI::dbIsValid(handle))) {
        record <- load_connection(deployment = deployment)
        creds  <- resolve_credentials(
            record$credentials_source,
            record$credentials_path
        )
        # bigint = "numeric" forces RPostgres to return BIGINT columns
        # as R double instead of the default bit64::integer64. The
        # latter silently misbehaves in base R (`max`, `cut`, ggplot
        # scale transforms etc. reinterpret the bit pattern as ~1e-317
        # / ~2e-23) and is the root cause of multiple FR-007a band /
        # bar-width bugs.
        handle <- DBI::dbConnect(
            RPostgres::Postgres(),
            host     = record$db_host,
            port     = record$db_port,
            dbname   = record$db_name,
            user     = creds$user,
            password = creds$password,
            bigint   = "numeric"
        )
        pool[[deployment]] <- handle
        logger::log_info(
            "Connected to {record$name} ",
            "({record$db_host}:{record$db_port}, role={record$role})"
        )
    }

    handle
}


#' Run a query against the deployment a (panel, facet) routes to
#'
#' Thin shim around \code{\link{pg_query}} that resolves the
#' deployment first, opens / re-uses a handle via
#' \code{\link{pg_connect_panel}}, runs the query, and tags the
#' returned tibble's attribute \code{"deployment"} with the
#' resolved label so the provenance sidecar can collect every
#' deployment the artifact touched (FR-031, sidecar
#' \code{deployments: array}).
#'
#' @param panel_id Character: panel identifier (see
#'     \code{\link{panel_deployment}}).
#' @param sql Character: SQL text.
#' @param ... Parameters bound via \code{DBI::dbGetQuery(params = ...)}.
#' @param facet Character or \code{NULL}: facet identifier.
#'
#' @return A tibble of query results carrying \code{"sql"},
#'     \code{"result_hash"}, and \code{"deployment"} attributes.
#'
#' @examples
#' \dontrun{
#' rows <- pg_query_panel(
#'     "fig01-overview",
#'     "SELECT COUNT(*) FROM entity",
#'     facet = "structures"      # routes to dev4
#' )
#' attr(rows, "deployment")
#' }
#'
#' @importFrom logger log_trace
#' @export
pg_query_panel <- function(panel_id, sql, ..., facet = NULL) {

    deployment <- panel_deployment(panel_id, facet)
    con <- pg_connect_panel(deployment)
    rows <- pg_query(con, sql, ...)

    attr(rows, "deployment") <- deployment
    rows
}


#' Close every cached panel-deployment handle
#'
#' Called by \code{rebuild.R} at the end of a run to release
#' connection slots. No-op if the pool is empty.
#'
#' @return Invisibly \code{NULL}.
#'
#' @importFrom DBI dbDisconnect dbIsValid
#' @importFrom logger log_info
#' @importFrom purrr walk
#' @export
pg_close_all_panel <- function() {

    pool <- panel_connection_pool()
    names_open <- names(pool)

    purrr::walk(names_open, function(name) {
        handle <- pool[[name]]
        if (!is.null(handle) && isTRUE(DBI::dbIsValid(handle))) {
            DBI::dbDisconnect(handle)
            logger::log_info("Disconnected from {name}")
        }
        rm(list = name, envir = pool)
    })

    invisible(NULL)
}


#' Internal panel-connection pool (package-private environment)
#'
#' @return The pool environment (created on first use).
#'
#' @keywords internal
#' @noRd
panel_connection_pool <- function() {
    if (!exists(".panel_pool", envir = .metabo_figures_state)) {
        assign(
            ".panel_pool",
            new.env(parent = emptyenv()),
            envir = .metabo_figures_state
        )
    }
    get(".panel_pool", envir = .metabo_figures_state)
}


#' Package-local state environment
#'
#' Holds the panel-connection pool (and any other per-session state
#' that future tasks may add). Created on package load via
#' \code{zzz.R}; this fallback constructor keeps the file safe to
#' source standalone during development.
#'
#' @keywords internal
#' @noRd
.metabo_figures_state <- new.env(parent = emptyenv())
