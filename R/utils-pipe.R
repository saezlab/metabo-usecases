#' Pipe operator
#'
#' Re-export of magrittr's `\%>\%` so internal callers in this package
#' resolve it at install / load time.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @importFrom magrittr %>%
#' @export
`%>%` <- magrittr::`%>%`


#' Compound assignment pipe
#'
#' Re-export of magrittr's `\%<>\%` for in-place updates.
#'
#' @name %<>%
#' @rdname compound-pipe
#' @keywords internal
#' @importFrom magrittr %<>%
#' @export
`%<>%` <- magrittr::`%<>%`
