#' Ratio mean square error
#'
#' @param x vector of observed values
#' @param y vector of simulated values
#'
#' @keywords internal
#' @rdname error
#' @export
#'
#' @examples
#'
#' rMSE(rnorm(10), rnorm(10))
#'
#' MSE(rnorm(10), rnorm(10))
#'
rMSE <- function(x, y) {

  sum((x/y - 1)^2)/length(y) ## ratio MSE
}


#' @keywords internal
#' @rdname error
#' @export

MSE <- function(x, y) {

  sum((x - y)^2)/length(y) ## ordinary MSE
}
