## distributions.R
## Merged from: ggamma.R, paretoII.R, burrXII.R, burrIII.R, gev.R
## Group C refactoring.  No algorithmic changes; style-only pass.
## norm.R belongs in fitting.R (Group F) — not included here.


# ===========================================================================
# Generalized Gamma (GGamma)
# ===========================================================================

#' Generalized Gamma distribution
#'
#' Provides density, distribution function, quantile function, random value
#' generation, and raw moments of order \emph{r} for the generalized gamma
#' distribution.
#'
#' @param x,q     vector of quantiles.
#' @param p       vector of probabilities.
#' @param n       number of observations. If \code{length(n) > 1}, the length
#'   is taken to be the number required.
#' @param r       raw moment order.
#' @param scale,shape1,shape2   scale and shape parameters; the shape
#'   arguments cannot be vectors (must have length one).
#' @param log,log.p   logical; if \code{TRUE}, probabilities \code{p} are
#'   given as \code{log(p)}.
#' @param lower.tail  logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return \code{dggamma} returns a numeric vector of density values.
#'   \code{pggamma} returns a numeric vector of cumulative probabilities.
#'   \code{qggamma} returns a numeric vector of quantiles.
#'   \code{rggamma} returns a numeric vector of random deviates.
#'   \code{mggamma} returns the raw moment of order \code{r}.
#'
#' @name GGamma
#' @aliases GGamma
#' @aliases dggamma
#'
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @import stats ggplot2
#' @export
#'
#' @references Papalexiou, S.M., Koutsoyiannis, D. (2012). Entropy based
#'   derivation of probability distributions: A case study to daily rainfall.
#'   Advances in Water Resources, 45, 51-57,
#'   \doi{10.1016/j.advwatres.2011.11.007}
#'
#' @seealso \code{\link{fitDist}}, \code{\link{moments}}
#'
#' @examples
#'
#' ## plot the density
#'
#' ggplot(data.frame(x = c(0, 20)),
#'        aes(x)) +
#'   stat_function(fun = dggamma,
#'                 args = list(scale = 5,
#'                             shape1 = .25,
#'                             shape2 = .75),
#'                 colour = "royalblue4") +
#'   labs(x = "",
#'        y = "Density") +
#'   theme_classic()

dggamma <- function(x, scale, shape1, shape2, log = FALSE) {

  if ((shape1 <= 0) | (shape2 <= 0)) {
    NaN
  } else {

    d <- (shape2 * (x / scale) ^ (-1 + shape1)) /
         (exp(x / scale) ^ shape2 * scale * gamma(shape1 / shape2))

    if (log) d <- log(d)

    d
  }
}

#' @rdname GGamma
#' @export

pggamma <- function(q, scale, shape1, shape2,
                    lower.tail = TRUE, log.p = FALSE) {

  if ((shape1 <= 0) | (shape2 <= 0)) {
    NaN
  } else {

    p <- pgamma((q / scale) ^ shape2, scale = 1, shape = shape1 / shape2)

    if (!lower.tail) p <- 1 - p
    if (log.p)       p <- log(p)

    p
  }
}

#' @rdname GGamma
#' @export

qggamma <- function(p, scale, shape1, shape2,
                    lower.tail = TRUE, log.p = FALSE) {

  if ((shape1 <= 0) | (shape2 <= 0)) {
    NaN
  } else {

    if (!lower.tail) p <- 1 - p
    if (log.p)       p <- log(p)

    scale * qgamma(p, scale = 1, shape = shape1 / shape2) ^ (1 / shape2)
  }
}

#' @rdname GGamma
#' @export

rggamma <- function(n, scale, shape1, shape2) {

  if ((shape1 <= 0) | (shape2 <= 0)) {
    NaN
  } else {
    qggamma(pnorm(rnorm(n)), scale, shape1, shape2)
  }
}

#' @rdname GGamma
#' @export

mggamma <- function(r, scale, shape1, shape2) {

  if ((shape1 <= 0) | (shape2 <= 0)) {
    NaN
  } else {
    (scale ^ r * gamma(r / shape2 + shape1 / shape2)) / gamma(shape1 / shape2)
  }
}


# ===========================================================================
# Pareto Type II (ParetoII)
# ===========================================================================

#' Pareto Type II distribution
#'
#' Provides density, distribution function, quantile function, random value
#' generation, and raw moments of order \emph{r} for the Pareto type II
#' distribution.
#'
#' @param x,q     vector of quantiles.
#' @param p       vector of probabilities.
#' @param n       number of observations. If \code{length(n) > 1}, the length
#'   is taken to be the number required.
#' @param r       raw moment order.
#' @param scale,shape   scale and shape parameters; the shape argument cannot
#'   be a vector (must have length one).
#' @param log,log.p   logical; if \code{TRUE}, probabilities \code{p} are
#'   given as \code{log(p)}.
#' @param lower.tail  logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return \code{dparetoII} returns a numeric vector of density values.
#'   \code{pparetoII} returns a numeric vector of cumulative probabilities.
#'   \code{qparetoII} returns a numeric vector of quantiles.
#'   \code{rparetoII} returns a numeric vector of random deviates.
#'   \code{mparetoII} returns the raw moment of order \code{r}.
#'
#' @name ParetoII
#' @aliases ParetoII
#' @aliases dparetoII
#'
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @import stats ggplot2
#' @export
#'
#' @seealso \code{\link{fitDist}}, \code{\link{moments}}
#'
#' @examples
#'
#' ## plot the density
#'
#' ggplot(data.frame(x = c(0, 20)),
#'        aes(x)) +
#'   stat_function(fun = dparetoII,
#'                 args = list(scale = 1,
#'                             shape = .3),
#'                 colour = "royalblue4") +
#'   labs(x = "",
#'        y = "Density") +
#'   theme_classic()

dparetoII <- function(x, scale, shape, log = FALSE) {

  if ((scale <= 0) | (shape <= 0)) {
    NaN
  } else {

    d <- (1 + (shape * x) / scale) ^ (-1 - 1 / shape) / scale

    if (log) d <- log(d)

    d
  }
}

#' @rdname ParetoII
#' @export

pparetoII <- function(q, scale, shape, lower.tail = TRUE, log.p = FALSE) {

  if ((scale <= 0) | (shape <= 0)) {
    NaN
  } else {

    p <- 1 - (1 + (shape * q) / scale) ^ (-1 / shape)

    if (!lower.tail) p <- 1 - p
    if (log.p)       p <- log(p)

    p
  }
}

#' @rdname ParetoII
#' @export

qparetoII <- function(p, scale, shape, lower.tail = TRUE, log.p = FALSE) {

  if ((scale <= 0) | (shape <= 0)) {
    NaN
  } else {

    if (!lower.tail) p <- 1 - p
    if (log.p)       p <- log(p)

    (scale * (-1 + (1 - p) ^ (-shape))) / shape
  }
}

#' @rdname ParetoII
#' @export

rparetoII <- function(n, scale, shape) {

  if ((scale <= 0) | (shape <= 0)) {
    NaN
  } else {
    qparetoII(runif(n), scale, shape)
  }
}

#' @rdname ParetoII
#' @export

mparetoII <- function(r, scale, shape) {

  if ((scale <= 0) | (shape <= 0)) {
    NaN
  } else {
    ((scale / shape) ^ r * gamma(1 / shape - r) * gamma(1 + r)) /
      gamma(1 / shape)
  }
}


# ===========================================================================
# Burr Type XII (BurrXII)
# ===========================================================================

#' Burr Type XII distribution
#'
#' Provides density, distribution function, quantile function, random value
#' generation, and raw moments of order \emph{r} for the Burr Type XII
#' distribution.
#'
#' @param x,q     vector of quantiles.
#' @param p       vector of probabilities.
#' @param n       number of observations. If \code{length(n) > 1}, the length
#'   is taken to be the number required.
#' @param r       raw moment order.
#' @param scale,shape1,shape2   scale and shape parameters; the shape
#'   arguments cannot be vectors (must have length one).
#' @param log,log.p   logical; if \code{TRUE}, probabilities \code{p} are
#'   given as \code{log(p)}.
#' @param lower.tail  logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return \code{dburrXII} returns a numeric vector of density values.
#'   \code{pburrXII} returns a numeric vector of cumulative probabilities.
#'   \code{qburrXII} returns a numeric vector of quantiles.
#'   \code{rburrXII} returns a numeric vector of random deviates.
#'   \code{mburrXII} returns the raw moment of order \code{r}.
#'
#' @name BurrXII
#' @aliases BurrXII
#' @aliases dburrXII
#'
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @import stats ggplot2
#' @export
#'
#' @references Papalexiou, S.M. (2018). Unified theory for stochastic
#'   modelling of hydroclimatic processes: Preserving marginal distributions,
#'   correlation structures, and intermittency. Advances in Water Resources,
#'   115, 234-252, \doi{10.1016/j.advwatres.2018.02.013}
#'
#' @seealso \code{\link{fitDist}}, \code{\link{moments}}
#'
#' @examples
#'
#' ## plot the density
#'
#' ggplot(data.frame(x = c(0, 10)),
#'        aes(x)) +
#'   stat_function(fun = dburrXII,
#'                 args = list(scale = 5,
#'                             shape1 = .25,
#'                             shape2 = .75),
#'                 colour = "royalblue4") +
#'   labs(x = "",
#'        y = "Density") +
#'   theme_classic()

dburrXII <- function(x, scale, shape1, shape2, log = FALSE) {

  if ((shape1 <= 0) | (shape2 <= 0)) {
    NaN
  } else {

    d <- ((x / scale) ^ (-1 + shape1) *
            (1 + shape2 * (x / scale) ^ shape1) ^ (-1 - 1 / (shape1 * shape2))) /
         scale

    if (log) d <- log(d)

    d
  }
}

#' @rdname BurrXII
#' @export

pburrXII <- function(q, scale, shape1, shape2,
                     lower.tail = TRUE, log.p = FALSE) {

  if ((shape1 <= 0) | (shape2 <= 0)) {
    NaN
  } else {

    p <- 1 - (1 + shape2 * (q / scale) ^ shape1) ^ (-1 / (shape1 * shape2))

    if (!lower.tail) p <- 1 - p
    if (log.p)       p <- log(p)

    p
  }
}

#' @rdname BurrXII
#' @export

qburrXII <- function(p, scale, shape1, shape2,
                     lower.tail = TRUE, log.p = FALSE) {

  if ((shape1 <= 0) | (shape2 <= 0)) {
    NaN
  } else {

    if (!lower.tail) p <- 1 - p
    if (log.p)       p <- log(p)

    ## Q(u) = beta * (((1-u)^(-g1*g2) - 1) / g2)^(1/g1)
    scale * (((1 - p) ^ (-(shape1 * shape2)) - 1) / shape2) ^ (1 / shape1)
  }
}

#' @rdname BurrXII
#' @export

rburrXII <- function(n, scale, shape1, shape2) {

  if ((shape1 <= 0) | (shape2 <= 0)) {
    NaN
  } else {
    qburrXII(runif(n), scale, shape1, shape2)
  }
}

#' @rdname BurrXII
#' @export

mburrXII <- function(r, scale, shape1, shape2) {

  if ((shape1 <= 0) | (shape2 <= 0)) {
    NaN
  } else {
    (scale ^ r *
       shape2 ^ (-1 - r / shape1) *
       beta((r + shape1) / shape1,
            (1 - r * shape2) / (shape1 * shape2))) /
      shape1
  }
}


# ===========================================================================
# Burr Type III (BurrIII)
# ===========================================================================

#' Burr Type III distribution
#'
#' Provides density, distribution function, quantile function, random value
#' generation, and raw moments of order \emph{r} for the Burr Type III
#' distribution.
#'
#' @param x,q     vector of quantiles.
#' @param p       vector of probabilities.
#' @param n       number of observations. If \code{length(n) > 1}, the length
#'   is taken to be the number required.
#' @param r       raw moment order.
#' @param scale,shape1,shape2   scale and shape parameters; the shape
#'   arguments cannot be vectors (must have length one).
#' @param log,log.p   logical; if \code{TRUE}, probabilities \code{p} are
#'   given as \code{log(p)}.
#' @param lower.tail  logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return \code{dburrIII} returns a numeric vector of density values.
#'   \code{pburrIII} returns a numeric vector of cumulative probabilities.
#'   \code{qburrIII} returns a numeric vector of quantiles.
#'   \code{rburrIII} returns a numeric vector of random deviates.
#'   \code{mburrIII} returns the raw moment of order \code{r}.
#'
#' @name BurrIII
#' @aliases BurrIII
#' @aliases dburrIII
#'
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @import stats ggplot2
#' @export
#'
#' @references Papalexiou, S.M. (2018). Unified theory for stochastic
#'   modelling of hydroclimatic processes: Preserving marginal distributions,
#'   correlation structures, and intermittency. Advances in Water Resources,
#'   115, 234-252, \doi{10.1016/j.advwatres.2018.02.013}
#'
#' @seealso \code{\link{fitDist}}, \code{\link{moments}}
#'
#' @examples
#'
#' ## plot the density
#'
#' ggplot(data.frame(x = c(1, 15)),
#'        aes(x)) +
#'   stat_function(fun = dburrIII,
#'                 args = list(scale = 5,
#'                             shape1 = .25,
#'                             shape2 = .75),
#'                 colour = "royalblue4") +
#'   labs(x = "",
#'        y = "Density") +
#'   theme_classic()

dburrIII <- function(x, scale, shape1, shape2, log = FALSE) {

  if ((shape1 <= 0) | (shape2 <= 0)) {
    NaN
  } else {

    d <- ((x / scale) ^ (-1 - 1 / shape2) *
            (1 + 1 / (shape1 * (x / scale) ^ (1 / shape2))) ^ (-1 - shape1 * shape2)) /
         scale

    if (log) d <- log(d)

    d
  }
}

#' @rdname BurrIII
#' @export

pburrIII <- function(q, scale, shape1, shape2,
                     lower.tail = TRUE, log.p = FALSE) {

  if ((shape1 <= 0) | (shape2 <= 0)) {
    NaN
  } else {

    p <- (1 + 1 / (shape1 * (q / scale) ^ (1 / shape2))) ^ (-(shape1 * shape2))

    if (!lower.tail) p <- 1 - p
    if (log.p)       p <- log(p)

    p
  }
}

#' @rdname BurrIII
#' @export

qburrIII <- function(p, scale, shape1, shape2,
                     lower.tail = TRUE, log.p = FALSE) {

  if ((shape1 <= 0) | (shape2 <= 0)) {
    NaN
  } else {

    if (!lower.tail) p <- 1 - p
    if (log.p)       p <- log(p)

    scale * (shape1 * (p ^ (-1 / (shape1 * shape2)) - 1)) ^ (-shape2)
  }
}

#' @rdname BurrIII
#' @export

rburrIII <- function(n, scale, shape1, shape2) {

  if ((shape1 <= 0) | (shape2 <= 0)) {
    NaN
  } else {
    qburrIII(runif(n), scale, shape1, shape2)
  }
}

#' @rdname BurrIII
#' @export

mburrIII <- function(r, scale, shape1, shape2) {

  if ((shape1 <= 0) | (shape2 <= 0)) {
    NaN
  } else {
    (scale ^ r *
       gamma((r + shape1) * shape2) *
       gamma(1 - r * shape2)) /
      (shape1 ^ (r * shape2) * gamma(shape1 * shape2))
  }
}


# ===========================================================================
# Generalized Extreme Value (GEV)
# ===========================================================================

#' Generalized Extreme Value distribution
#'
#' Provides density, distribution function, quantile function, random value
#' generation, and raw moments of order \emph{r} for the generalized extreme
#' value distribution.
#'
#' @param x,q     vector of quantiles.
#' @param p       vector of probabilities.
#' @param n       number of observations. If \code{length(n) > 1}, the length
#'   is taken to be the number required.
#' @param r       raw moment order.
#' @param loc,scale,shape   location, scale, and shape parameters.
#' @param log,log.p   logical; if \code{TRUE}, probabilities \code{p} are
#'   given as \code{log(p)}.
#' @param lower.tail  logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return \code{dgev} returns a numeric vector of density values.
#'   \code{pgev} returns a numeric vector of cumulative probabilities.
#'   \code{qgev} returns a numeric vector of quantiles.
#'   \code{rgev} returns a numeric vector of random deviates.
#'   \code{mgev} returns the raw moment of order \code{r} (via numerical
#'   integration).
#'
#' @name GEV
#' @aliases GEV
#' @aliases dgev
#'
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @import stats ggplot2
#' @export
#'
#' @seealso \code{\link{fitDist}}, \code{\link{moments}}
#'
#' @examples
#'
#' ## plot the density
#'
#' ggplot(data.frame(x = c(0, 20)),
#'        aes(x)) +
#'   stat_function(fun = dgev,
#'                 args = list(loc = 1,
#'                             scale = .5,
#'                             shape = .15),
#'                 colour = "royalblue4") +
#'   labs(x = "",
#'        y = "Density") +
#'   theme_classic()

dgev <- function(x, loc, scale, shape, log = FALSE) {

  x <- (x - loc) / scale

  if (shape == 0) {

    d <- log(1 / scale) - x - exp(-x)
  } else {

    n     <- length(x)
    x     <- 1 + shape * x
    i     <- x > 0
    scale <- rep(scale, length.out = n)[i]

    d     <- c()
    d[i]  <- log(1 / scale) - x[i] ^ (-1 / shape) - (1 / shape + 1) * log(x[i])
    d[!i] <- -Inf
  }

  if (!log) d <- exp(d)

  d
}

#' @rdname GEV
#' @export

pgev <- function(q, loc, scale, shape, lower.tail = TRUE, log.p = FALSE) {

  if (shape == 0) {
    p <- (q - loc) / scale
  } else {
    p <- -1 / shape * log(pmax(0, 1 - shape * (q - loc) / scale))
  }

  p <- exp(-exp(-p))

  if (!lower.tail) p <- 1 - p
  if (log.p)       p <- log(p)

  p
}

#' @rdname GEV
#' @export

qgev <- function(p, loc, scale, shape, lower.tail = TRUE, log.p = FALSE) {

  if (!lower.tail) p <- 1 - p
  if (log.p)       p <- log(p)

  if (shape == 0) {
    loc - scale * log(-log(p))
  } else {
    loc + scale / shape * (1 - (-log(p)) ^ shape)
  }
}

#' @rdname GEV
#' @export

rgev <- function(n, loc, scale, shape) {

  qgev(p     = runif(n),
       loc   = loc,
       scale = scale,
       shape = shape)
}

#' @rdname GEV
#' @export

mgev <- function(r, loc, scale, shape) {

  if (scale == 0) {
    NaN
  } else {

    integrand <- function(x, q, ...) {
      x ^ q * dgev(x, loc = loc, scale = scale, shape = shape)
    }

    ## integration bounds vary with the sign of loc (GEV support boundary)
    if (loc == 0) {
      lwr <- -Inf
      upr <- Inf
    } else if (loc > 0) {
      lwr <- loc - scale / shape
      upr <- Inf
    } else {
      lwr <- -Inf
      upr <- loc - scale / shape
    }

    integrate(f             = integrand,
              lower         = lwr,
              upper         = upr,
              q             = r,
              loc           = loc,
              scale         = scale,
              shape         = shape,
              stop.on.error = FALSE)$value
  }
}
