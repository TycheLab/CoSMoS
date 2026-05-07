#' Numerical estimation of moments
#'
#' Uses numerical integration to compute the theoretical raw or central moments
#' of the specified distribution.
#'
#' @param dist character; distribution name (e.g. \code{"norm"},
#'   \code{"paretoII"})
#' @param distarg list of distribution arguments
#' @param distbounds numeric vector of length 2; distribution bounds (default
#'   \code{c(-Inf, Inf)})
#' @param p0 numeric; probability zero (default 0)
#' @param raw logical; compute raw moments?
#' @param central logical; compute central moments?
#' @param coef logical; compute standardised coefficients (CV, skewness,
#'   kurtosis)?
#' @param order integer vector; raw moment orders (default \code{1:4})
#'
#' @return a named list with zero or more of:
#'   \describe{
#'     \item{\code{m}}{raw moments}
#'     \item{\code{mu}}{central moments}
#'     \item{\code{coefficients}}{CV, skewness, kurtosis}
#'   }
#'
#' @seealso \code{\link{sample.moments}}, \code{\link{populationstat}}
#'
#' @import stats
#' @export
#'
#' @keywords moments
#'
#' @examples
#'
#' library(CoSMoS)
#'
#' ## Normal distribution
#' moments("norm", list(mean = 2, sd = 1))
#'
#' ## Pareto type II
#' moments(dist    = "paretoII",
#'         distarg = list(shape = 0.2, scale = 1))
#'
moments <- function(dist, distarg, p0 = 0, raw = TRUE, central = TRUE,
                    coef = TRUE, distbounds = c(-Inf, Inf), order = 1:4) {

  if ((central | coef) & (max(order) < 4)) order <- 1:4

  if (existsFunction(paste0("m", dist), where = asNamespace("CoSMoS"))) {
    m <- lapply(order, function(i) do.call(paste0("m", dist),
                                           args = c(list(r = i), distarg)))
  } else {
    f <- function(x, dist, q, ...) {
      .arg <- c(list(x), ...)
      x^q * do.call(paste0("d", dist), args = .arg)
    }
    m <- vector("list", max(order))
    for (i in order) {
      m[[i]] <- integrate(f = f, distbounds[1], distbounds[2],
                          q = i, distarg, dist = dist,
                          stop.on.error = FALSE)$value
    }
  }

  m        <- lapply(m, function(x) (1 - p0) * x)
  names(m) <- paste0("m", order)

  out <- list()

  if (raw) out$m <- unlist(m)

  if (central | coef) {
    mu <- data.frame(
      mu1 = m$m1,
      mu2 = m$m2 - m$m1^2,
      mu3 = 2 * m$m1^3 - 3 * m$m2 * m$m1 + m$m3,
      mu4 = -3 * m$m1^4 + 6 * m$m2 * m$m1^2 - 4 * m$m3 * m$m1 + m$m4
    )
    if (central) out$mu <- unlist(mu)
  }

  if (coef) {
    coef_df <- data.frame(
      cv = sqrt(mu$mu2) / mu$mu1,
      cs = mu$mu3 / sqrt(mu$mu2)^3,
      ck = mu$mu4 / mu$mu2^2
    )
    out$coefficients <- unlist(coef_df)
  }

  out
}


#' Sample moments
#'
#' Computes raw moments, central moments, and standardised coefficients (CV,
#' skewness, kurtosis) from a numeric sample.
#'
#' @param x numeric vector of values
#' @param na.rm logical; strip \code{NA} values before computation?
#' @inheritParams moments
#'
#' @return a named list with zero or more of:
#'   \describe{
#'     \item{\code{m}}{raw moments}
#'     \item{\code{mu}}{central moments}
#'     \item{\code{coefficients}}{CV, skewness, kurtosis}
#'   }
#'
#' @seealso \code{\link{moments}}, \code{\link{checkTS}}
#'
#' @import stats
#' @export
#'
#' @keywords moments
#'
#' @examples
#'
#' library(CoSMoS)
#'
#' x <- rnorm(1000)
#' sample.moments(x)
#'
#' y <- rparetoII(1000, 10, .1)
#' sample.moments(y)
#'
sample.moments <- function(x, na.rm = FALSE, raw = TRUE, central = TRUE,
                           coef = TRUE, order = 1:4) {

  if ((central | coef) & (max(order) < 4)) order <- 1:4

  n    <- if (na.rm) length(which(!is.na(x))) else length(x)
  xbar <- mean(x, na.rm = na.rm)
  out  <- list()

  if (raw) {
    m        <- lapply(order, function(i) 1 / n * sum(x^i))
    names(m) <- paste0("m", order)
    out$m    <- unlist(m)
  }

  if (central | coef) {
    mu        <- lapply(order, function(i) 1 / n * sum((x - xbar)^i))
    names(mu) <- paste0("mu", order)
    if (central) out$mu <- unlist(mu)
  }

  if (coef) {
    coef_df <- data.frame(
      cv = sqrt(mu$mu2) / mu$mu1,
      cs = mu$mu3 / sqrt(mu$mu2)^3,
      ck = mu$mu4 / mu$mu2^2
    )
    out$coefficients <- unlist(coef_df)
  }

  out
}
