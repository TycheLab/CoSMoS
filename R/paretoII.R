#' Pareto type II distribution
#'
#' Provides density, distribution function, quantile function, random value generation
#' and raw moments of order \emph{r} for the Pareto type II distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param r               raw moment order
#' @param scale,shape     scale and shape parameters; the shape argument cannot be a vector (must have length one).
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.
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
#' @examples
#'
#' ## plot the density
#'
#' ggplot(data.frame(x = c(0, 20)),
#'        aes(x)) +
#'   stat_function(fun = dparetoII,
#'                 args = list(scale = 1,
#'                             shape = .3),
#'                 colour = 'royalblue4') +
#'   labs(x = '',
#'        y = 'Density') +
#'   theme_classic()
#'

dparetoII <- function(x, scale, shape, log = FALSE) { ## density

  if((scale <= 0) | (shape <= 0)) {

    return(NaN)
  } else {

    d <- (1 + (shape*x)/scale)^(-1 - 1/shape)/scale

    if (log) {d <- log(d)}

    return(d)
  }
}

#' @rdname ParetoII
#' @export

pparetoII <- function(q, scale, shape, lower.tail = TRUE, log.p = FALSE) { ## cdf

  if((scale <= 0) | (shape <= 0)) {

    return(NaN)
  } else {

    p <- 1 - (1 + (shape*q)/scale)^(-1/shape)

    if (!lower.tail) {p <- 1 - p}
    if (log.p) {p <- log(p)}

    return(p)
  }
}

#' @rdname ParetoII
#' @export

qparetoII <- function(p, scale, shape, lower.tail = TRUE, log.p = FALSE) { ## quantile function

  if((scale <= 0) | (shape <= 0)) {

    return(NaN)
  } else {

    if (!lower.tail) {p <- 1 - p}
    if (log.p) {p <- log(p)}

    q <- (scale*(-1 + (1 - p)^(-shape)))/shape

    return(q)
  }
}

#' @rdname ParetoII
#' @export

rparetoII <- function(n, scale, shape) { ## random numbers generator

  if((scale <= 0) | (shape <= 0)) {

    return(NaN)
  } else {

    return(qparetoII(pnorm(rnorm(n)), scale, shape))
  }
}

#' @rdname ParetoII
#' @export

mparetoII <- function(r, scale, shape) {

  if((scale <= 0) | (shape <= 0)) {

    return(NaN)
  } else {

    return(((scale/shape)^r*gamma(1/shape - r)*gamma(1 + r))/gamma(1/shape))
  }
}

