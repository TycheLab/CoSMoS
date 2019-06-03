#' Empirical cummulative distrubution function
#'
#' Calculates ecdf based on weibull plotting position
#'
#' @param x vector of values
#'
#' @keywords internal
#' @export
#'
#' @examples
#'
#' ECDF(round(rnorm(100)))
#'
ECDF <- function(x) {

  ## sort the values: x1 < x2 < x3 < ... < xN
  st <- sort(x)

  ## weibull plotting position
  aux <- data.table(p = rank(
    x = st,
    ties.method = 'first') / (length(st) + 1),
                    value = st,
                    key = 'value')

  J <- value <- NULL

  out <- aux[J(unique(value)),
             mult = 'last']

  structure(.Data = out)
}
