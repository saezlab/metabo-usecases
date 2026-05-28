#' The lead categorical palette
#'
#' Returns the 9-colour OmniPath Metabo lead palette in canonical order
#' per spec FR-020. The last entry, \code{"#BEBEBE"}, is reserved for
#' "None / Unknown / NA" and is named accordingly.
#'
#' @return Named character vector of length 9 — hex codes in canonical
#'     order.
#'
#' @examples
#' palette_lead()
#'
#' @export
palette_lead <- function() {
    c(
        teal    = "#006384",
        magenta = "#9F0162",
        amber   = "#FEAF16",
        lime    = "#BBCC33",
        coral   = "#EA6572",
        green   = "#009E73",
        sky     = "#99DDFF",
        pink    = "#D03293",
        unknown = "#BEBEBE"
    )
}


#' The rwth extended palette
#'
#' Parses \code{inst/extdata/palettes/rwth.gpl} into a named hex
#' vector. Used for categorical needs that exceed the lead palette
#' (FR-021).
#'
#' @return Named character vector of hex codes.
#'
#' @examples
#' head(palette_rwth())
#'
#' @importFrom readr read_lines
#' @importFrom stringr str_match str_squish
#' @importFrom purrr discard map_chr map
#' @export
palette_rwth <- function() {
    gpl_path <- system.file(
        "extdata", "palettes", "rwth.gpl",
        package = "metabo.figures",
        mustWork = TRUE
    )
    parse_gpl_palette(gpl_path)
}


#' Pick N colours from the lead palette, honouring the unknown rule
#'
#' Returns the first \code{n} entries of \code{\link{palette_lead}}.
#' When \code{unknown = TRUE}, the last position is replaced by
#' \code{"#BEBEBE"} so the palette always carries a clear None/Unknown
#' slot (FR-020). When \code{n == 1L}, returns \code{"#BEBEBE"} alone
#' (single-colour-defaults-to-grey rule, FR-020).
#'
#' @param n Integer: number of categories to colour.
#' @param unknown Logical: whether one of the N categories is null-like
#'     (None / Unknown / NA) — if so, the last slot becomes
#'     \code{"#BEBEBE"}.
#'
#' @return Unnamed character vector of \code{n} hex codes.
#'
#' @examples
#' palette_n(4L)
#' palette_n(4L, unknown = TRUE)
#' palette_n(1L)
#'
#' @importFrom rlang abort
#' @export
palette_n <- function(n, unknown = FALSE) {

    if (!is.integer(n) || length(n) != 1L || n < 1L) {
        rlang::abort("`n` must be a positive integer scalar")
    }

    lead <- palette_lead()

    if (n == 1L) {
        unname(lead["unknown"])
    } else if (isTRUE(unknown)) {
        if (n > length(lead)) {
            rlang::abort(sprintf(
                "n = %d exceeds lead-palette size (%d)",
                n, length(lead)
            ))
        }
        unname(c(lead[seq_len(n - 1L)], lead["unknown"]))
    } else {
        if (n >= length(lead)) {
            rlang::abort(sprintf(
                "n = %d would reach the reserved unknown slot; ",
                "set unknown = TRUE if that's intended"
            ))
        }
        unname(lead[seq_len(n)])
    }
}


#' Parse a GIMP .gpl palette into a named hex vector
#'
#' @param path Character: path to the .gpl file.
#' @return Named character vector of hex codes.
#'
#' @importFrom readr read_lines
#' @importFrom stringr str_squish str_match
#' @importFrom purrr discard map_chr keep
#' @keywords internal
#' @noRd
parse_gpl_palette <- function(path) {

    body <-
        readr::read_lines(path) %>%
        stringr::str_squish() %>%
        purrr::discard(function(x) {
            x == "" || startsWith(x, "#") ||
                startsWith(x, "GIMP Palette") ||
                startsWith(x, "Name:") ||
                startsWith(x, "Columns:")
        })

    matches <- stringr::str_match(
        body, "^(\\d+)\\s+(\\d+)\\s+(\\d+)\\s*(.*)$"
    )
    rows <- matches[!is.na(matches[, 1L]), , drop = FALSE]

    hex <- purrr::map_chr(seq_len(nrow(rows)), function(i) {
        sprintf(
            "#%02X%02X%02X",
            as.integer(rows[i, 2L]),
            as.integer(rows[i, 3L]),
            as.integer(rows[i, 4L])
        )
    })
    names(hex) <- ifelse(
        rows[, 5L] == "",
        sprintf("c%02d", seq_along(hex)),
        gsub("[^a-zA-Z0-9_]+", "_", tolower(rows[, 5L]))
    )

    hex
}
