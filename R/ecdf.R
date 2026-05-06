#' Empirical cumulative distribution function
#'
#' Calculates the ECDF using the Weibull plotting position formula.
#'
#' @param x Numeric vector of values.
#'
#' @return A \code{data.table} with columns \code{p} (plotting positions)
#'   and \code{value} (sorted unique values).
#'
#' @keywords internal
#'
ECDF <- function(x) {

  st <- sort(x)

  aux <- data.table(
    p     = rank(x = st, ties.method = "first") / (length(st) + 1),
    value = st,
    key   = "value"
  )

  J <- value <- NULL  # data.table global variable declaration

  aux[J(unique(value)), mult = "last"]
}
