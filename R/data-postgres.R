#' Open a Postgres handle to the configured OmniPath deployment
#'
#' Reads connection parameters via \code{\link{load_connection}} plus
#' credentials per \code{credentials_source}: \code{env},
#' \code{file}, or \code{pgpass}. Returns a live DBI handle the caller
#' is responsible for closing (\code{DBI::dbDisconnect}).
#'
#' @return A \code{DBIConnection} to the OmniPath Postgres.
#'
#' @examples
#' \dontrun{
#' con <- pg_connect()
#' on.exit(DBI::dbDisconnect(con), add = TRUE)
#' }
#'
#' @importFrom DBI dbConnect
#' @importFrom RPostgres Postgres
#' @importFrom logger log_info
#' @importFrom rlang abort
#' @export
pg_connect <- function() {

    conn <- load_connection()
    creds <- resolve_credentials(conn$credentials_source, conn$credentials_path)

    # bigint = "numeric" returns BIGINT as R double, avoiding the
    # bit64::integer64 default which silently misbehaves in base R
    # and ggplot scale transforms (see pg_connect_panel docstring).
    handle <- DBI::dbConnect(
        RPostgres::Postgres(),
        host     = conn$db_host,
        port     = conn$db_port,
        dbname   = conn$db_name,
        user     = creds$user,
        password = creds$password,
        bigint   = "numeric"
    )

    logger::log_info(
        "Connected to {conn$name} ({conn$db_host}:{conn$db_port})"
    )

    handle
}


#' Resolve Postgres credentials per the configured source
#'
#' @param source One of \code{env}, \code{file}, \code{pgpass}.
#' @param path Optional path required when \code{source == "file"}.
#'
#' @return A list with \code{user}, \code{password} (\code{NA} when
#'     pgpass is used and libpq resolves credentials itself).
#'
#' @importFrom yaml read_yaml
#' @importFrom rlang abort
#' @keywords internal
#' @noRd
resolve_credentials <- function(source, path = NULL) {

    if (identical(source, "env")) {
        list(
            user     = Sys.getenv("PGUSER",     unset = ""),
            password = Sys.getenv("PGPASSWORD", unset = "")
        )
    } else if (identical(source, "file")) {
        if (is.null(path) || !file.exists(path)) {
            rlang::abort(sprintf(
                "credentials_source = file requires credentials_path; ",
                "got '%s'", path
            ))
        }
        creds <- yaml::read_yaml(path)
        list(user = creds$user, password = creds$password)
    } else if (identical(source, "pgpass")) {
        list(user = NA_character_, password = NA_character_)
    } else {
        rlang::abort(sprintf("Unknown credentials_source: '%s'", source))
    }
}


#' Run a query, return a tibble, and record the result for provenance
#'
#' Wraps \code{DBI::dbGetQuery} so every query the pipeline executes is
#' logged and fingerprinted. The returned tibble carries two
#' attributes: \code{"sql"} (the query text) and \code{"result_hash"}
#' (SHA-256 truncated to 12 hex chars over the canonical result
#' representation). The provenance sidecar writer picks these up.
#'
#' @param con A DBI Postgres connection from \code{\link{pg_connect}}.
#' @param sql Character: parameterized SQL.
#' @param ... Parameters bound via \code{DBI::dbGetQuery(params = ...)}.
#'
#' @return A tibble of query results with \code{sql} and
#'     \code{result_hash} attributes.
#'
#' @examples
#' \dontrun{
#' con <- pg_connect()
#' rows <- pg_query(con, "SELECT COUNT(*) FROM data_source")
#' attr(rows, "result_hash")
#' }
#'
#' @importFrom DBI dbGetQuery
#' @importFrom tibble as_tibble
#' @importFrom digest digest
#' @importFrom logger log_info
#' @export
pg_query <- function(con, sql, ...) {

    params <- list(...)
    rows <- tibble::as_tibble(
        if (length(params) == 0L) {
            DBI::dbGetQuery(con, sql)
        } else {
            DBI::dbGetQuery(con, sql, params = params)
        }
    )

    hash <- substr(
        digest::digest(rows, algo = "sha256", serialize = TRUE),
        1L, 12L
    )

    attr(rows, "sql")         <- trimws(sql)
    attr(rows, "result_hash") <- hash

    logger::log_info(
        "query rows={nrow(rows)} hash={hash}"
    )

    rows
}
