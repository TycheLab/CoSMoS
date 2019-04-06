#' Burr Type XII distribution
#'
#' Density, distribution function, quantile function, random value generation,
#' and raw moments of order \emph{r} for the Burr Type XII distribution.
#'
#' @param x,q	                  vector of quantiles.
#' @param p	                    vector of probabilities.
#' @param n	                    number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param r                     raw moment order
#' @param scale,shape1,shape2   scale and shape parameters; the shape arguments cannot be a vectors (must have length one).
#' @param log,log.p	            logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	          logical; if TRUE (default), probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.
#'
#' @name BurrXII
#' @aliases BurrXII
#' @aliases dburrXII
#'
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @import stats
#' @export

dburrXII <- function(x, scale, shape1, shape2, log = FALSE) { ## density

  if((shape1 <= 0) | (shape2 <= 0)) {

    return(NaN)
  } else {

    d <- ((x/scale)^(-1 + shape1)*(1 + shape2*(x/scale)^shape1)^(-1 - 1/(shape1*shape2)))/scale

    if (log) {d <- log(d)}

    return(d)
  }
}

#' @rdname BurrXII
#' @export

pburrXII <- function(q, scale, shape1, shape2, lower.tail = TRUE, log.p = FALSE) { ## cdf

  if((shape1 <= 0) | (shape2 <= 0)) {

    return(NaN)
  } else {

    p <- 1 - (1 + shape2*(q/scale)^shape1)^(-(1/(shape1*shape2)))

    if (!lower.tail) {p <- 1 - p}
    if (log.p) {p <- log(p)}

    return(p)
  }
}

#' @rdname BurrXII
#' @export

qburrXII <- function(p, scale, shape1, shape2, lower.tail = TRUE, log.p = FALSE) { ## cdf

  if((shape1 <= 0) | (shape2 <= 0)) {

    return(NaN)
  } else {

    if (!lower.tail) {p <- 1 - p}
    if (log.p) {p <- log(p)}

    q <- scale*(-((1 - (1 - p)^(-(shape1*shape2)))/shape2))^shape1^(-1)

    return(q)
  }
}

#' @rdname BurrXII
#' @export

rburrXII <- function(n, scale, shape1, shape2) {

  if((shape1 <= 0) | (shape2 <= 0)) {

    return(NaN)
  } else {

    return(qburrXII(pnorm(rnorm(n)), scale, shape1, shape2))
  }
}

#' @rdname BurrXII
#' @export

mburrXII <- function(r, scale, shape1, shape2) {

  if((shape1 <= 0) | (shape2 <= 0)) {

    return(NaN)
  } else {

    return((scale^r*shape2^(-1 - r/shape1)*beta((r + shape1)/shape1, (1 - r*shape2)/(shape1*shape2)))/shape1)
  }
}

