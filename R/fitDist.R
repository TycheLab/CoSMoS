#' Distribution fitting
#'
#' Uses Nelder-Mead simplex algorithm to minimize fitting norms.
#'
#' @inheritParams N
#' @param data value to be fitted
#' @param constrain logical - constrain shape2 parametes for finite tails
#'
#' @export
#' @import nloptr
#'
#' @examples
#'
#' x <- fitDist(rnorm(1000), 'norm', 30, 'N1', FALSE)
#' x
#'
fitDist <- function(data, dist, n.points, norm, constrain) {

  a <- getDistArg(dist)

  lwr <- NULL

  if (constrain) {

    if (any(a == 'shape2')) {

      start <- c(rep(1, length(a) - 1),
                 .1)
      upr <- c(rep(Inf, length(a) - 1),
               .499)
    } else {

      start <- rep(1, length(a))
      upr <- rep(Inf, length(a))
    }

  } else {

    if (any(a == 'shape2')) {

      start <- c(rep(1, length(a) - 1),
                 .1)
      upr <- rep(Inf, length(a))
    } else {

      start <- rep(1, length(a))
      upr <- rep(Inf, length(a))
    }
  }

  fit <-  neldermead(x0 = start, ## parameter starting position
                     fn = N, ## function to minimize
                     val = data, ## input data
                     dist = dist, ## distribution to be fitted
                     n.points = n.points, ## number of points for fitting
                     norm = norm, ## norm used to fit
                     lower = lwr, ## lower bound
                     upper = upr) ## upper bound

  out <- as.list(setNames(fit$par, ## return list of dist arguments
                          a))

  structure(.Data = out,
            dist = dist,
            edf = ECDF(data),
            err = fit$value,
            class = 'fitDist')
}

#' Plot method for fitDist
#'
#' @param x fitDist object
#' @param ... other args
#'
#' @keywords internal
#' @export
#' @import ggplot2
#' @method plot fitDist
#'
#' @examples
#'
#' x <- fitDist(rnorm(1000), 'norm', 30, 'N1', FALSE)
#'
#' plot(x)
#'
plot.fitDist <- function(x, ...) {

  edf <- attr(x, 'edf') ## get ecdf
  dist <- attr(x, 'dist') ## get dist name
  err <- attr(x, 'err')

  mm <- do.call(paste0('p', dist), ## find the min / max prob for the plot
                c(list(q = range(edf$value)),
                  x))
  p <- exp(seq(log(mm[1]), ## make seq of propabilities
               log(mm[2]),
               length.out = 10000))
  cdf <- data.frame(p = p, ## calculate theoretical cdf
                    value = do.call(paste0('q', dist),
                                    c(list(p = p), x)))

  p <- ggplot() + ## plot the values
    geom_line(data = cdf, ## theoretical
              aes(x = cdf$value,
                  y = log(1 - cdf$p)),
              colour = 'grey25',
              lwd = 1,
              alpha = .75) +
    geom_point(data = edf, ## empirical
               aes(x = edf$value,
                   y = log(1 - edf$p)),
               colour = 'red4',
               alpha = .5) +
    labs(x = 'Nonzero values',
         y = 'Exceedence probability',
         title = paste('Fitting norm (error) value =', round(err, 5))) +
    scale_y_continuous(breaks = seq(-10, ## auxiliary axis ticks
                                    0,
                                    length.out = 5),
                       labels = format(exp(seq(log(.0001),
                                               log(1),
                                               length.out = 5)),
                                       scientific = TRUE)) +
    theme_grey()

  return(p)
}
