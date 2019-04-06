#' Generate timeseries
#'
#' Generates timeseries with given properties, just provide (1) the target marginal distribution and its parameters, (2) the target autocorrelation structure or individual autocorrelation values up to a desired lag, and (3) the probablility zero if you wish to simulate an intermittent process.
#'
#' @inheritParams actpnts
#' @inheritParams ARp
#' @param TSn number of timeseries to be generated
#'
#' @import reshape2 ggplot2
#' @export
#'
#' @examples
#'
#' library(CoSMoS)
#'
#' ## define target and marginal distributions and arguments
#' x <- generateTS(margdist = 'burrXII',
#'                 margarg = list(scale = 1,
#'                                shape1 = .75,
#'                                shape2 = .25),
#'                 acsvalue = acs(id = 'weibull',
#'                                t = 0:30,
#'                                scale = 10,
#'                                shape = .75),
#'                 n = 1000, p = 30, p0 = .5, TSn = 3)
#'
#' ## see the results
#' dta <- data.frame(n = seq_along(x[[1]]), do.call(cbind, x))
#' names(dta)[2:(length(x) + 1)] <- paste('timeseries', 1:length(x), sep = '_')
#'
#' m.dta <- melt(dta, id.vars = 'n')
#'
#' ggplot(m.dta, aes(x = n, y = value)) +
#'   geom_line() +
#'   labs(x = '') +
#'   facet_wrap(~variable, ncol = 1) +
#'   theme_classic()

#'
generateTS <- function(margdist, margarg, n, p, p0 = 0, TSn = 1, distbounds = c(-Inf, Inf), limitorder = TRUE, acsvalue = NULL) {

  pnts <- actpnts(margdist = margdist, ## estimate act points
                  margarg = margarg,
                  p0 = p0,
                  distbounds = distbounds)

  par <- fitactf(pnts) ## fit actf to the points

  out <- lapply(1:TSn, function(i) { ## get timeseries using ARp

    ARp(margdist = margdist,
        margarg = margarg,
        n = n,
        p = p,
        p0 = p0,
        actfpara = par,
        limitorder = limitorder,
        acsvalue = acsvalue)
  })

  return(out)
}

#' Bulk Timeseries generation
#'
#' You have used the generateTS function and you wish to generate more time series. Instead of re-running the generateTS you can use the regenerateTS Timeseries that used all the parameters calculated by the generateTS function and thus it is faster.
#'
#' @param ts generated timeseries using ARp
#' @param TSn number of timeseries to be (re)generated
#' @export
#'
#' @examples
#'
#' library(CoSMoS)
#'
#' ## define marginal distribution and arguments with target autocorrelation structure
#' x <- generateTS(margdist = 'burrXII',
#'                 margarg = list(scale = 1,
#'                                shape1 = .75,
#'                                shape2 = .25),
#'                 acsvalue = acs(id = 'weibull',
#'                                t = 0:30,
#'                                scale = 10,
#'                                shape = .75),
#'                 n = 1000, p = 30, p0 = .5, TSn = 3)
#'
#' ## generate new values with same parameters
#' regenerateTS(x)
#'
regenerateTS <- function(ts, TSn = 1) {

  X <- ts ## argument rename

  ## getting all necesary info from attributes

  margdist <- attr(X[[1]], 'margdist')
  margarg <- attr(X[[1]], 'margarg')

  n <- length(X[[1]])
  p <- length(attr(X[[1]], 'a'))
  a <- attr(X[[1]], 'a')
  p0 <- attr(X[[1]], 'p0')
  esd <- attr(X[[1]], 'esd')

  out <- lapply(1:TSn, function(i) {

    gn <- c(rep(0, (2*p)), rnorm(n, mean = 0, sd = esd)) ## generate gaussian noise

    val <- c(AR1((2*p), alpha = a[1]), rep(0, n)) ## generate vector of values

    for (i in (p + 1):(n + (2*p))) { ## AR
      val[i] <- sum(val[(i - p):(i - 1)]*a) + gn[i]
    }

    uval <- (pnorm(val[-1:-(2*p)]) - p0)/(1 - p0) ## p0 + gaussian probabilities calculation
    uval[uval < 0] <- 0

    do.call(paste0('q', margdist), args = c(list(p = uval), margarg)) ## quantile mapping
  })

  return(out)
}
