#' Numerical estimation of moments
#'
#' Uses numerical integration to caclulate the theoretical raw or central moments of the specified distribution
#'
#' @param dist distribution
#' @param distarg list of distribution arguments
#' @param distbounds distribution bounds (default set to c(-Inf, Inf))
#' @param p0 probability zero
#' @param raw logical - calculate raw moments?
#' @param central logical - calculate central moments?
#' @param coef logical - calculate coeffitients (coeffitient of variation, skweness and kurtosis)?
#' @param order vector of integers - raw moment orders
#'
#' @import stats
#' @export
#'
#' @keywords moments, moment
#'
#' @examples
#'
#' library(CoSMoS)
#'
#' ## Normal Distribution
#' moments('norm', list(mean = 2, sd = 1))
#'
#' ## Pareto type II
#' scale <- 1
#' shape <- .2
#'
#' moments(dist = 'paretoII',
#'         distarg = list(shape = shape,
#'                        scale = scale))
#'
moments <- function(dist, distarg, p0 = 0, raw = T, central = T, coef = T, distbounds = c(-Inf, Inf), order = 1:4) {

  if((central | coef) & (max(order) < 4)) {

    order <- 1:4
  }

  # tryCatch({ ## error handeling
    if(exists(paste0('m', dist))) { ## condition - if raw moments function is specified, use it

      m <- lapply(order, function(i) do.call(paste0('m', dist), args = c(list(r = i), distarg)))

    } else { ## - else integrate the PDF

      f <- function(x, dist, q, ...){

        .arg <- c(list(x), ...) ## pull the distribution arguments

        x^q * do.call(paste0('d', dist), args = .arg) ## function to integrate
      }

      m <- list()

      for (i in order) {

        m[[i]] <- integrate(f = f, distbounds[1], distbounds[2], q = i, distarg, dist = dist, stop.on.error = F)$value
      }
    }

    m <- lapply(m, function(x) (1 - p0)*x)

    names(m) <- paste0('m', order)

    out <- list()

    if (raw) {

      out$m <- unlist(m)
    }

    if (central | coef) { ## central moments calculation

      mu <- data.frame(mu1 = m$m1,
                       mu2 = m$m2 - m$m1^2,
                       mu3 = 2*m$m1^3 - 3*m$m2*m$m1 + m$m3,
                       mu4 = -3*m$m1^4 + 6*m$m2*m$m1^2 - 4*m$m3*m$m1 + m$m4)

      if (central) {

        out$mu <- unlist(mu)
      }
    }

    if (coef) { ## coeffitiens calculation

      coef <- data.frame(cv = sqrt(mu$mu2)/mu$mu1,
                         cs = mu$mu3/sqrt(mu$mu2)^3,
                         ck = mu$mu4/mu$mu2^2)

      out$coefficients <- unlist(coef)
    }

    return(out)
  # },
  # error = function(e) {
  #   message('Please input a correct distribution name with valid parameters \n(ditributions have to be defined in form dxxx, pxxx, qxxx, etc...)')
  # }
  # )
}

#' Estimation of sample moments
#'
#' @param x a numeric vector of values
#' @param na.rm a logical value indicating whether NA values should be stripped before the computation proceeds
#' @inheritParams moments
#'
#' @import stats
#' @export
#'
#' @keywords moments, moment
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
sample.moments <- function(x, na.rm = FALSE, raw = T, central = T, coef = T, order = 1:4) {

  if((central | coef) & (max(order) < 4)) {

    order <- 1:4
  }

  n <- ifelse(na.rm, length(which(!is.na(x))), length(x))

  mean <- mean(x, na.rm = na.rm)

  out <- list()

  if (raw) {

    m <- lapply(order, function(i) (1/n * sum((x)^i)))
    names(m) <- paste0('m', order)

    out$m <- unlist(m)
  }

  if (central | coef) { ## central moments calculation

    mu <- lapply(order, function(i) 1/n * sum((x - mean)^i))
    names(mu) <- paste0('mu', order)

    if (central) {

      out$mu <- unlist(mu)
    }
  }

  if (coef) { ## coeffitiens calculation

    coef <- data.frame(cv = sqrt(mu$mu2)/mu$mu1,
                       cs = mu$mu3/sqrt(mu$mu2)^3,
                       ck = mu$mu4/mu$mu2^2)

    out$coefficients <- unlist(coef)
  }

  return(out)
}

