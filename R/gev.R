#' Generalized extreme value distribution
#'
#' Provides density, distribution function, quantile function,
#' and random value generation, for the generalized extreme value distribution.
#'
#' @param x,q	                  vector of quantiles.
#' @param p	                    vector of probabilities.
#' @param n	                    number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param r                     raw moment order
#' @param loc,scale,shape       location, scale and shape parameters.
#' @param log,log.p	            logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	          logical; if TRUE (default), probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.
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
#'                 colour = 'royalblue4') +
#'   labs(x = '',
#'        y = 'Density') +
#'   theme_classic()

dgev <- function (x, loc, scale, shape, log = FALSE) {

  x <- (x - loc)/scale

  if (shape == 0) {

    d <- log(x = 1 / scale) - x - exp(-x)
  } else {

    n <- length(x = x)
    x <- 1 + shape * x
    i <- x > 0
    scale <- rep(x = scale,
                 length.out = n)[i]

    d <- c()
    d[i] <- log(1 / scale) - x[i] ^ (-1 / shape) - (1 / shape + 1) * log(x[i])
    d[!i] <- -Inf
  }

  if (!log) {
    d <- exp(d)
  }

  return(d)
}

#' @rdname GEV
#' @export

pgev <- function(q, loc, scale, shape, lower.tail = TRUE, log.p = FALSE) {

  if (shape == 0) {

    p <- (q - loc) / scale
  } else {

    p <- -1 / shape * log(pmax(0, 1 - shape*(q - loc) / scale))
  }

  p <- exp(-exp(-p))

  if (!lower.tail) {

    p <- 1 - p
  }
  if (log.p) {

    p <- log(p)
  }

  return(p)
}

#' @rdname GEV
#' @export

qgev <- function(p, loc, scale, shape, lower.tail = TRUE, log.p = FALSE) {

  if (!lower.tail) {

    p <- 1 - p
  }
  if (log.p) {

    p <- log(p)
  }

  if (shape == 0) {

    q <- loc - scale * log(-log(p))
  } else {

    q <- loc + scale/shape*(1 - (-log(p)) ^ shape)
  }

  return(q)
}

#' @rdname GEV
#' @export

rgev <- function(n, loc, scale, shape) {

  qgev(p = pnorm(rnorm(n)),
       loc = loc,
       scale = scale,
       shape =shape)
}

#' @rdname GEV
#' @export

mgev <- function(r, loc, scale, shape) {

  if(scale == 0) {

    return(NaN)
  } else {

    f <- function(x, q, ...){

      x^q*dgev(x, loc = loc, scale = scale, shape = shape)
    }

    if (loc == 0) {

      upr <- Inf
      lwr <- -Inf
    }

    if (loc > 0) {

      upr <- Inf
      lwr <- loc - scale / shape
    }

    if (loc < 0) {

      upr <- loc - scale / shape
      lwr <- -Inf
    }

    m <- integrate(f = f,
                   lower = lwr,
                   upper = upr,
                   q = r,
                   loc = loc,
                   scale = scale, shape = shape,
                   stop.on.error = F)$value

    return(m)
  }
}
