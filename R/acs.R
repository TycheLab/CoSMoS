#' AutoCorrelation Structure
#'
#' Provides a parametric function that describes the values of the linear autocorrelation up to desired lags. For more details on the parametric autocorrelation structures see section 3.2 in Papalexiou (2018).
#'
#' @param id autocorrelation structure id
#' @param ... other arguments (t as lag and acs parameters)
#'
#' @name acs
#'
#' @import ggplot2
#' @export
#'
#' @references Papalexiou, S.M. (2018). Unified theory for stochastic modelling of hydroclimatic
#' processes: Preserving marginal distributions, correlation structures,
#' and intermittency. Advances in Water Resources, 115, 234-252, \doi{10.1016/j.advwatres.2018.02.013}
#'
#' @examples
#'
#' library(CoSMoS)
#'
#' ## specify lag
#' t <- 0:10
#'
#' ## get the ACS
#' f <- acs('fgn', t = t, H = .75)
#' b <- acs('burrXII', t = t, scale = 1, shape1 = .6, shape2 = .4)
#' w <- acs('weibull', t = t, scale = 2, shape = 0.8)
#' p <- acs('paretoII', t = t, scale = 3, shape = 0.3)
#'
#' ## visualize the ACS
#' dta <- data.table(t, f, b, w, p)
#'
#' m.dta <- melt(dta, id.vars = 't')
#'
#' ggplot(m.dta,
#'        aes(x = t,
#'            y = value,
#'            group = variable,
#'            colour = variable)) +
#'   geom_point(size = 2.5) +
#'   geom_line(lwd = 1) +
#'   scale_color_manual(values = c('steelblue4', 'red4', 'green4', 'darkorange'),
#'                      labels = c('FGN', 'Burr XII', 'Weibull', 'Pareto II'),
#'                      name = '') +
#'   labs(x = bquote(lag ~ tau),
#'        y = 'Acf') +
#'   scale_x_continuous(breaks = t) +
#'   theme_classic()
#'
acs <- function(id, ...) {

  .args <- list(...)

  do.call(paste0('acf', id), args = .args)
}

#' Autocorrelation structure functions
#'
#' @keywords internal
#' @name ACSfunctions
#' @export
#'
#' @examples
#'
#' t <- 1
#'
#' H <- .75
#' scale <- .2
#' shape <- .3
#' shape1 <- .5
#' shape2 <- .2
#'
#' acfburrXII(t, scale, shape1, shape2)
#'
#' acfparetoII(t, scale, shape)
#'
#' acffgn(t, H)
#'
#' acfweibull(t, scale, shape)
#'

acfburrXII <- function(t, scale, shape1, shape2) {

  if((shape1 <= 0) | (shape2 <= 0)) {

    return(NaN)
  } else {

    return((1 + shape2 * (t / scale) ^ shape1) ^ -(1 / shape1 * shape2))
  }
}

#' @keywords internal
#' @rdname ACSfunctions
#' @export

acfparetoII <- function(t, scale, shape) {

  if ((scale <= 0) | (shape <= 0)) {

    return(NaN)
  } else {

    return((1 + (shape * t) / scale) ^ (-shape ^ (-1)))
  }
}

#' @keywords internal
#' @rdname ACSfunctions
#' @export


acffgn <- function(t, H) {

  (abs(-1 + t) ^ (2 * H) - 2 * abs(t) ^ (2 * H) + abs(1 + t) ^ (2 * H)) / 2
}

#' @keywords internal
#' @rdname ACSfunctions
#' @export

acfweibull <- function(t, scale, shape) {
  exp(-(t / scale) ^ shape)
}
