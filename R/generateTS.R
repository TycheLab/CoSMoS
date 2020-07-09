#' Generate timeseries
#'
#' Generates timeseries with given properties, just provide (1) the target marginal
#' distribution and its parameters, (2) the target autocorrelation structure or
#' individual autocorrelation values up to a desired lag, and (3) the probablility
#' zero if you wish to simulate an intermittent process.
#'
#' @details
#'
#' A step-by-step guide:
#' * First define the target marginal (margdist), that is, the probability distribution
#' of the generated data. For example set margdist = 'ggamma' if you wish to generate
#' data following the Generalized Gamma distribution, margidst = 'burrXII' for Burr
#' type XII distribution etc. For a full list of the distributions we support see the
#' help vignette: \code{vignette('vignette', package = 'CoSMoS')}. In general, the package
#' supports all build-in distribution functions of R and of other packages.
#'
#' * Define the parameters’ values (margarg) of the distribution you selected. For example
#' the Generalized Gamma has one scale and two shape parameters so set the desired value,
#' e.g., margarg = list(scale = 2, shape1 = 0.9, shape2 = 0.8). Note distributions might
#' have different number of parameters and different type of parameters (location, scale, shape).
#' To see the parameters of each distribution we support, see the help vignette:
#' \code{vignette('vignette', package = 'CoSMoS')}.
#'
#' * If you wish your time series to be intermittent (e.g., precipitation), then define the
#' probability zero. For example, set p0 = 0.9, if you wish your generated data to have
#' 90\% of zero values (dry days).
#'
#' * Define your linear autocorrelations.
#'     + You can supply specific lag autocorrelations starting from lag 0
#'     and up to a desired lag, e.g., acs = c(1, 0.9, 0.8, 0.7); this will generate
#'     a process with lag1, 2 and 3 autocorrelations equal with 0.9, 0.8 and 0.7.
#'
#'     + Alternatively, you can use a parametric autocorrelation structure (see section 3.2 in
#'     \href{https://doi.org/10.1016/j.advwatres.2018.02.013}{Papalexiou 2018}).
#'     We support the following autocorrelation structures (acs) weibull, paretoII,
#'     fgn and burrXII. See also \link[CoSMoS]{acs} examples.
#'
#' * Define the order to the autoregressive model p. For example if you aim to preserve
#' the first 10 lag autocorrelations then just set p = 10. Otherwise set it p = NULL and
#' the model will decide the value of p in order to preserve the whole autocorrelation
#' structure.
#'
#' * Lastly just define the time series length, e.g., n = 1000 and number of time series
#' you wish to generate, e.g., TSn = 10.
#'
#' Play around with the following given examples which will make the whole
#' process a piece of cake.
#'
#' @inheritParams actpnts
#' @inheritParams ARp
#' @param TSn number of timeseries to be generated
#'
#' @import data.table ggplot2
#' @export
#'
#' @examples
#'
#' library(CoSMoS)
#'
#' ## Case1:
#' ## You wish to generate 3 time series of size 1000 each
#' ## that follow the Generalized Gamma distribution with parameters
#' ## scale = 1, shape1 = 0.8, shape2 = 0.8
#' ## and autocorrelation structure the ParetoII
#' ## with parameters scale = 1 and shape = .75
#' x <- generateTS(margdist = 'ggamma',
#'                 margarg = list(scale = 1,
#'                                shape1 = .8,
#'                                shape2 = .8),
#'                 acsvalue = acs(id = 'paretoII',
#'                                t = 0:30,
#'                                scale = 1,
#'                                shape = .75),
#'                 n = 1000,
#'                 p = 30,
#'                 TSn = 3)
#'
#' ## see the results
#' plot(x)
#'
#' \donttest{
#'
#' ## Case2:
#' ## You wish to generate time series the same distribution
#' ## and autocorrelations as is Case1 but intermittent
#' ## with probability zero equal to 90%
#' y <- generateTS(margdist = 'ggamma',
#'                 margarg = list(scale = 1,
#'                                shape1 = .8,
#'                                shape2 = .8),
#'                 acsvalue = acs(id = 'paretoII',
#'                                t = 0:30,
#'                                scale = 1,
#'                                shape = .75),
#'                 p0 = .9,
#'                 n = 1000,
#'                 p = 30,
#'                 TSn = 3)
#'
#' ## see the results
#' plot(y)
#'
#' ## Case3:
#' ## You wish to generate a time series of size 1000
#' ## that follows the Beta distribution
#' ## (e.g., relative humidity ranging from 0 to 1)
#' ## with parameters shape1 = 0.8, shape2 = 0.8, is defined from 0 to 1
#' ## and autocorrelation structure the ParetoII
#' ## with parameters scale = 1 and shape = .75
#' z <- generateTS(margdist = 'beta',
#'                 margarg = list(shape1 = .6,
#'                                shape2 = .8),
#'                 distbounds = c(0, 1),
#'                 acsvalue = acs(id = 'paretoII',
#'                                t = 0:30,
#'                                scale = 1,
#'                                shape = .75),
#'                 n = 1000,
#'                 p = 20)
#'
#' ## see the results
#' plot(z)
#'
#' ## Case4:
#' ## Same in previous case but now you provide specific
#' ## autocorrelation values for the first three lags,
#' ## ie.., lag 1 to 3 equal to 0.9, 0.8 and 0.7
#'
#' z <- generateTS(margdist = 'beta',
#'                 margarg = list(shape1 = .6,
#'                                shape2 = .8),
#'                 distbounds = c(0, 1),
#'                 acsvalue = c(1, .9, .8, .7),
#'                 n = 1000,
#'                 p = TRUE)
#'
#' ## see the results
#' plot(z)
#'
#'}
#'
generateTS <- function(n, margdist, margarg, p = NULL, p0 = 0, TSn = 1, distbounds = c(-Inf, Inf), acsvalue = NULL) {

  pnts <- actpnts(margdist = margdist, ## estimate act points
                  margarg = margarg,
                  p0 = p0,
                  distbounds = distbounds)

  par <- fitactf(pnts) ## fit actf to the points

  out <- lapply(1:TSn, function(i) { ## get timeseries using ARp

    return(
      value = ARp(margdist = margdist,
                  margarg = margarg,
                  n = n,
                  p = p,
                  p0 = p0,
                  actfpara = par,
                  acsvalue = acsvalue)
    )


  })

  structure(.Data = out,
            class = 'cosmosts')
}

#' Bulk Timeseries generation
#'
#' Resamples given Timeseries
#'
#' @details
#'
#' You have used the generateTS function and you wish to generate more time series. Instead of re-running generateTS you can use regenerateTS, which generates timeseries using the parameters previously calculated by the generateTS function,  and thus it is faster.
#'
#' @param ts generated timeseries using ARp
#' @param TSn number of timeseries to be (re)generated
#' @export
#'
#' @examples
#'
#' library(CoSMoS)
#'
#' ## define marginal distribution and arguments with target
#' ## autocorrelation structure
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
#' r <- regenerateTS(x)
#'
#' plot(r)
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

  structure(.Data = out,
            class = 'cosmosts')
}
