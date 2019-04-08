#' Generalized gamma distribution
#'
#' Provides density, distribution function, quantile function, random value generation,
#' and raw moments of order \emph{r} for the generalized gamma distribution.
#'
#' @param x,q	                  vector of quantiles.
#' @param p	                    vector of probabilities.
#' @param n	                    number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param r                     raw moment order
#' @param scale,shape1,shape2   scale and shape parameters; the shape arguments cannot be a vectors (must have length one).
#' @param log,log.p	            logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	          logical; if TRUE (default), probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.
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
#'                 colour = 'royalblue4') +
#'   labs(x = '',
#'        y = 'Density') +
#'   theme_classic()


dggamma <- function(x, scale, shape1, shape2, log = FALSE) { ## density

  if((shape1 <= 0) | (shape2 <= 0)) {

    return(NaN)
  } else {

    d <- (shape2*(x/scale)^(-1 + shape1))/(exp(x/scale)^shape2*scale*gamma(shape1/shape2))

    if (log) {d <- log(d)}

    return(d)
  }
}

#' @rdname GGamma
#' @export

pggamma <- function(q, scale, shape1, shape2, lower.tail = TRUE, log.p = FALSE) { ## cdf

  if((shape1 <= 0) | (shape2 <= 0)) {

    return(NaN)
  } else {

    p <- pgamma((q/scale)^shape2, scale = 1, shape = shape1/shape2)

    if (!lower.tail) {p <- 1 - p}
    if (log.p) {p <- log(p)}

    return(p)
  }
}

#' @rdname GGamma
#' @export

qggamma <- function(p, scale, shape1, shape2, lower.tail = TRUE, log.p = FALSE) { ## cdf

  if((shape1 <= 0) | (shape2 <= 0)) {

    return(NaN)
  } else {

    if (!lower.tail) {p <- 1 - p}
    if (log.p) {p <- log(p)}

    q <- scale*qgamma(p, scale = 1, shape = shape1/shape2)^(1/shape2)

    return(q)
  }
}

#' @rdname GGamma
#' @export

rggamma <- function(n, scale, shape1, shape2) {

  if((shape1 <= 0) | (shape2 <= 0)) {

    return(NaN)
  } else {

    return(qggamma(pnorm(rnorm(n)), scale, shape1, shape2))
  }
}

#' @rdname GGamma
#' @export

mggamma <- function(r, scale, shape1, shape2) {

  if((shape1 <= 0) | (shape2 <= 0)) {

    return(NaN)
  } else {

    return((scale^r*gamma(r/shape2 + shape1/shape2))/gamma(shape1/shape2))
  }
}
