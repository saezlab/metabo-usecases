#' Category-colour registry environment
#'
#' Holds the mapping from recurring categorical values (resource names,
#' interaction types, ontologies, etc.) to specific hex colours from the
#' lead or rwth palettes. Per FR-022 / SC-004: the same category always
#' renders in the same colour across every figure.
#'
#' Stored in an environment so the registry can be extended at
#' package-load time (see \code{\link{register_category_colours}}) while
#' still being shared across every figure script that imports it.
#'
#' @keywords internal
#' @noRd
.metabo_category_colours <- new.env(parent = emptyenv())


#' Initialize the registry with the canonical mappings
#'
#' Called from \code{.onLoad()}. Idempotent.
#'
#' @keywords internal
#' @noRd
init_category_colours <- function() {

    # All categories start empty. Figure scripts register the values
    # they observe in the live snapshot via
    # \code{\link{register_category_colours}}; the assertive accessor
    # \code{\link{assert_category_known}} fails fast on unregistered
    # values so missing registrations surface as build errors rather
    # than silently autocoloring (spec Edge Case).
    #
    # Once Figure 4's COSMOS PKN comparison lands (Milestone B), we
    # will hard-code the curated interaction_types palette so its
    # categories render the same colour across every figure.
    .metabo_category_colours$interaction_types <- character(0)
    .metabo_category_colours$resources         <- character(0)
    .metabo_category_colours$ontologies        <- character(0)

    invisible()
}


#' List the values currently registered in a category
#'
#' Lets figure scripts compute the "new" subset they need to
#' \code{\link{register_category_colours}} for, without poking at the
#' internal registry environment directly.
#'
#' @param category Character scalar.
#' @return Character vector of registered value names (empty when the
#'     category has no entries yet).
#'
#' @examples
#' registered_category_values("interaction_types")
#'
#' @export
registered_category_values <- function(category) {
    if (!exists(category, envir = .metabo_category_colours)) {
        return(character(0))
    }
    names(.metabo_category_colours[[category]])
}


#' Get the colour for a (category, value) pair
#'
#' @param category Character scalar: one of the registered category
#'     names (e.g. \code{"interaction_types"}, \code{"resources"}).
#' @param value Character: one or more category values to look up.
#'
#' @return Character vector of hex codes, same length as \code{value}.
#'
#' @examples
#' register_category_colours(
#'     "interaction_types",
#'     c(signaling = "#006384", transport = "#D03293")
#' )
#' category_colour("interaction_types", c("signaling", "transport"))
#'
#' @importFrom rlang abort
#' @export
category_colour <- function(category, value) {

    assert_category_known(category, value)
    mapping <- .metabo_category_colours[[category]]
    unname(mapping[value])
}


#' Assert that every value is registered in the named category
#'
#' Fails fast per the spec Edge Case "new category appears in the data
#' that has no entry in the colour registry": pipeline MUST fail with a
#' clear error so the registry is updated deliberately.
#'
#' @param category Character scalar.
#' @param values Character vector.
#'
#' @return Invisibly \code{TRUE} on success; \code{rlang::abort()} on
#'     missing entries.
#'
#' @importFrom rlang abort
#' @export
assert_category_known <- function(category, values) {

    if (!exists(category, envir = .metabo_category_colours)) {
        rlang::abort(sprintf(
            "Category '%s' is not registered. Known: %s",
            category,
            paste(ls(.metabo_category_colours), collapse = ", ")
        ))
    }

    mapping <- .metabo_category_colours[[category]]
    missing <- setdiff(unique(values), names(mapping))

    if (length(missing) > 0L) {
        rlang::abort(sprintf(
            paste0(
                "Category '%s' is missing colour entries for: %s. ",
                "Register them in R/style/category_colours.R ",
                "(spec Edge Case: never silently autocolor)."
            ),
            category, paste(missing, collapse = ", ")
        ))
    }

    invisible(TRUE)
}


#' Extend the registry with additional category-value → colour entries
#'
#' Lets figure scripts (e.g. R/plots/db_content.R) add resource-colour
#' mappings discovered from the live snapshot before they render. The
#' added entries are merged into the existing category mapping; an
#' attempt to overwrite an existing entry triggers
#' \code{rlang::abort()} to surface accidental palette drift.
#'
#' @param category Character scalar.
#' @param mapping Named character vector: name = category value,
#'     value = hex code.
#'
#' @return Invisibly the updated mapping for the category.
#'
#' @importFrom rlang abort
#' @export
register_category_colours <- function(category, mapping) {

    current <- if (exists(category, envir = .metabo_category_colours)) {
        .metabo_category_colours[[category]]
    } else {
        character(0)
    }

    overlap <- intersect(names(current), names(mapping))
    conflict <- overlap[current[overlap] != mapping[overlap]]
    if (length(conflict) > 0L) {
        rlang::abort(sprintf(
            "Refusing to overwrite category-colour entries: %s",
            paste(conflict, collapse = ", ")
        ))
    }

    merged <- c(current, mapping[!names(mapping) %in% names(current)])
    .metabo_category_colours[[category]] <- merged

    invisible(merged)
}
