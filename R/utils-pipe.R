#' Pipe operator
#'
#' Re-export of \code{magrittr::\link[magrittr]{\%>\%}} so internal
#' uses of the pipe in this package's source resolve at install /
#' load time.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or magrittr placeholder.
#' @param rhs A function call using magrittr semantics.
#' @return The result of \code{rhs(lhs)}.
NULL


#' Compound assignment pipe operator
#'
#' Re-export of \code{magrittr::\link[magrittr]{\%<>\%}} for in-place
#' updates.
#'
#' @name %<>%
#' @rdname compound-pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %<>%
#' @usage lhs \%<>\% rhs
#' @param lhs A value or magrittr placeholder.
#' @param rhs A function call.
#' @return Used for side effects on \code{lhs}.
NULL
