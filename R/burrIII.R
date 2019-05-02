#' Burr Type III distribution
#'
#' Provides density, distribution function, quantile function, random value generation,
#' and raw moments of order \emph{r} for the Burr Type III distribution.
#'
#' @param x,q	                  vector of quantiles.
#' @param p	                    vector of probabilities.
#' @param n	                    number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param r                     raw moment order
#' @param scale,shape1,shape2   scale and shape parameters; the shape arguments cannot be a vectors (must have length one).
#' @param log,log.p	            logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	          logical; if TRUE (default), probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.
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
#'                 colour = 'royalblue4') +
#'   labs(x = '',
#'        y = 'Density') +
#'   theme_classic()

dburrIII <- function(x, scale, shape1, shape2, log = FALSE) { ## density

  if((shape1 <= 0) | (shape2 <= 0)) {

    return(NaN)
  } else {

    d <- ((x/scale)^(-1 - shape2^(-1))*(1 + 1/(shape1*(x/scale)^shape2^(-1)))^(-1 - shape1*shape2))/scale

    if (log) {d <- log(d)}

    return(d)
  }
}

#' @rdname BurrIII
#' @export

pburrIII <- function(q, scale, shape1, shape2, lower.tail = TRUE, log.p = FALSE) { ## cdf

  if((shape1 <= 0) | (shape2 <= 0)) {

    return(NaN)
  } else {

    p <- (1 + 1/(shape1*(q/scale)^shape2^(-1)))^(-(shape1*shape2))

    if (!lower.tail) {p <- 1 - p}
    if (log.p) {p <- log(p)}

    return(p)
  }
}

#' @rdname BurrIII
#' @export

qburrIII <- function(p, scale, shape1, shape2, lower.tail = TRUE, log.p = FALSE) { ## quntile

  if((shape1 <= 0) | (shape2 <= 0)) {

    return(NaN)
  } else {

    if (!lower.tail) {p <- 1 - p}
    if (log.p) {p <- log(p)}

    q <- scale*(-((1 - (1 - p)^(-(shape1*shape2)))/shape2))^shape1^(-1)

    return(q)
  }
}

#' @rdname BurrIII
#' @export

rburrIII <- function(n, scale, shape1, shape2) {

  if((shape1 <= 0) | (shape2 <= 0)) {

    return(NaN)
  } else {

    return(qburrIII(pnorm(rnorm(n)), scale, shape1, shape2))
  }
}

#' @rdname BurrIII
#' @export

mburrIII <- function(r, scale, shape1, shape2) {

  if((shape1 <= 0) | (shape2 <= 0)) {

    return(NaN)
  } else {

    return((scale^r*gamma((r + shape1)*shape2)*gamma(1 - r*shape2))/(shape1^(r*shape2)*gamma(shape1*shape2)))
  }
}
