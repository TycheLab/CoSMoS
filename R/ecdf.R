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

  st <- sort(x) ## sort the values: x1 < x2 < x3 < ... < xN
  p <- rank(st)/(length(st) + 1) ## weibull plotting position

  out <- unique(data.table(p = p,
                           value = st))

  structure(.Data = out)
}
