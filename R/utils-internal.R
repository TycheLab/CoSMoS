# Internal utility functions for CoSMoS
# These functions are not exported and not intended for direct use.

#' Compute spectral values for DHM circulant-embedding simulation
#'
#' Computes the one-sided spectral values \eqn{S_j} from an autocovariance
#' function (ACVF) vector using the Dietrich-Newsam-Haskard circulant-embedding
#' approach.
#'
#' @param acvf numeric vector of autocovariance values starting at lag 0
#'
#' @return numeric vector of \eqn{S_j} values (length = \code{length(acvf)})
#'
#' @keywords internal
DHMgenSj <- function(acvf) {
  MM <- length(acvf) - 1
  Sj <- Re(fft(c(acvf, rev(acvf[2:MM]))))[1:(MM + 1)]
  Sj <- abs(Sj)
  if (!all(Sj >= 0))
    stop("some of the S_j's are negative")
  Sj
}

#' Simulate one realisation using DHM circulant-embedding
#'
#' Simulates one realisation of a stationary Gaussian process from precomputed
#' spectral values \eqn{S_j} via the Dietrich-Newsam-Haskard algorithm.
#'
#' @param Sj numeric vector of spectral values from \code{DHMgenSj}
#' @param rn optional vector of standard normal random numbers of length
#'   \code{2 * length(Sj) - 2}; generated internally if \code{NULL}
#'
#' @return numeric vector of length \code{length(Sj) - 1}
#'
#' @keywords internal
DHMgenSim <- function(Sj, rn = NULL) {
  if (is.null(rn))
    rn <- rnorm(2 * length(Sj) - 2)
  M  <- 2 * length(Sj) - 2
  N  <- M / 2
  Yj <- 0:N
  Yj[1]     <- sqrt(Sj[1]) * rn[1]
  Yj[N + 1] <- sqrt(Sj[N + 1]) * rn[M]
  js         <- 2:N
  Yj[js]    <- sqrt(0.5 * Sj[js]) *
    complex(real      = rn[2 * (1:(N - 1))],
            imaginary = rn[2 * (1:(N - 1)) + 1])
  Re(fft(c(Yj, Conj(rev(Yj[js])))))[1:N] / sqrt(M)
}


#' L-moments
#'
#' Computes the first four L-moments of a sample using probability-weighted
#' moments.
#'
#' @param x numeric vector of values
#'
#' @return named numeric vector of length 4: \code{l1}, \code{l2},
#'   \code{l3} (L-skewness ratio), \code{l4} (L-kurtosis ratio)
#'
#' @seealso \code{\link{reportTS}}
#' @keywords internal
lmom <- function(x) {

  x  <- sort(x)
  n  <- length(x)
  nn <- rep(n - 1, n)
  pp <- seq(0, n - 1)

  p1 <- pp / nn
  p2 <- p1 * (pp - 1) / (nn - 1)
  p3 <- p2 * (pp - 2) / (nn - 2)

  b0 <- sum(x) / n
  b1 <- sum(p1 * x) / n
  b2 <- sum(p2 * x) / n
  b3 <- sum(p3 * x) / n

  l1 <- b0
  l2 <- 2 * b1 - b0
  l3 <- 2 * (3 * b2 - b0) / (2 * b1 - b0) - 3
  l4 <- 5 * (2 * (2 * b3 - 3 * b2) + b0) / (2 * b1 - b0) + 6

  unlist(mget(paste0("l", 1:4)))
}


#' Complementary error function
#'
#' @param x numeric vector
#'
#' @return numeric vector
#'
#' @rdname erfc
#' @keywords internal
erfc <- function(x) {
  2 * pnorm(-sqrt(2) * x)
}

#' Inverse complementary error function
#'
#' @rdname erfc
#' @keywords internal
inv.erfc <- function(x) {
  x[x < 0 | x > 2] <- NA
  -qnorm(x / 2) / sqrt(2)
}
