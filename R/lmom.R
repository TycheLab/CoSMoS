#' L-Moments calculation
#'
#' @param x vector of values
#'
#' @keywords internal
#' @export
#'
#' @examples
#'
#' lmom(rnorm(100))
#'
lmom <- function (x) { ## calculate l moments

  x <- sort(x)

  n <- length(x)
  nn <- rep(n - 1, n)

  pp <- seq(0, n - 1)
  p1 <- pp/nn
  p2 <- p1*(pp - 1)/(nn - 1)
  p3 <- p2*(pp - 2)/(nn - 2)

  b0 <- sum(x)/n
  b1 <- sum(p1*x)/n
  b2 <- sum(p2*x)/n
  b3 <- sum(p3*x)/n

  l1 <- b0
  l2 <- 2*b1 - b0
  l3 <- 2*(3*b2 - b0)/(2*b1 - b0) - 3
  l4 <- 5*(2*(2*b3 - 3*b2) + b0)/(2*b1 - b0) + 6

  unlist(mget(paste0('l', 1:4)))
}
