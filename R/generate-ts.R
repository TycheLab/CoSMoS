#' Generate time series
#'
#' Generates time series with given properties. Provide (1) the target marginal
#' distribution and its parameters, (2) the target autocorrelation structure or
#' individual autocorrelation values up to a desired lag, and (3) the probability
#' zero if you wish to simulate an intermittent process.
#'
#' @details
#'
#' A step-by-step guide:
#' * First define the target marginal (\code{margdist}), that is, the probability
#' distribution of the generated data. For example set \code{margdist = 'ggamma'}
#' for the Generalised Gamma distribution, \code{margdist = 'burrXII'} for Burr
#' type XII etc. For a full list of supported distributions see the help
#' \href{https://CRAN.R-project.org/package=CoSMoS/vignettes/vignette.html}{vignette}.
#' In general, the package supports all built-in distribution functions of R and
#' of other packages.
#'
#' * Define the parameters (\code{margarg}) of the selected distribution.
#' For example the Generalised Gamma has one scale and two shape parameters,
#' e.g. \code{margarg = list(scale = 2, shape1 = 0.9, shape2 = 0.8)}.
#' See the help vignette for details on each distribution's parameters.
#'
#' * If you wish your time series to be intermittent (e.g. precipitation), define
#' the probability zero. For example \code{p0 = 0.9} produces 90\% zero values.
#'
#' * Define your linear autocorrelations.
#'   + Supply specific lag autocorrelations starting from lag 0 up to a desired
#'     lag, e.g. \code{acsvalue = c(1, 0.9, 0.8, 0.7)}; this preserves lag-1,
#'     lag-2 and lag-3 autocorrelations equal to 0.9, 0.8 and 0.7.
#'
#'   + Alternatively, use a parametric autocorrelation structure (see section 3.2
#'     in Papalexiou (2018)). Supported structures: \code{weibull}, \code{paretoII},
#'     \code{fgn} and \code{burrXII}. See also \code{\link[CoSMoS]{acs}}.
#'
#' * Define the AR model order \code{p}. For example if you aim to preserve the
#' first 10 lag autocorrelations then set \code{p = 10}. Set \code{p = NULL} and
#' the model will choose \code{p} to preserve the whole autocorrelation structure.
#'
#' * Set the time series length, e.g. \code{n = 1000}, and the number of time
#' series to generate, e.g. \code{TSn = 10}.
#'
#' @inheritParams actpnts
#' @inheritParams ARp
#' @param TSn number of time series to generate
#'
#' @return An object of class \code{'cosmosts'}: a list of \code{TSn} numeric
#'   vectors, each of length \code{n}, with per-series attributes recording the
#'   fitted model parameters.
#'
#' @seealso \code{\link{regenerateTS}}, \code{\link{ARp}}, \code{\link{actpnts}}
#'
#' @import data.table ggplot2
#' @export
#'
#' @references Papalexiou, S.M. (2018). Unified theory for stochastic modelling of
#' hydroclimatic processes: Preserving marginal distributions, correlation structures,
#' and intermittency. Advances in Water Resources, 115, 234-252,
#' \doi{10.1016/j.advwatres.2018.02.013}
#'
#'
#' @examples
#'
#' library(CoSMoS)
#'
#' ## Case 1:
#' ## Generate 3 time series of length 1000 following the Generalised Gamma
#' ## distribution with scale = 1, shape1 = 0.8, shape2 = 0.8 and ParetoII
#' ## autocorrelation structure with scale = 1 and shape = 0.75.
#' x <- generateTS(margdist = "ggamma",
#'                 margarg = list(scale = 1,
#'                                shape1 = .8,
#'                                shape2 = .8),
#'                 acsvalue = acs(id = "paretoII",
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
#' ## Case 2:
#' ## Same as Case 1 but intermittent with probability zero equal to 90%.
#' y <- generateTS(margdist = "ggamma",
#'                 margarg = list(scale = 1,
#'                                shape1 = .8,
#'                                shape2 = .8),
#'                 acsvalue = acs(id = "paretoII",
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
#' ## Case 3:
#' ## Generate a time series of length 1000 following the Beta distribution
#' ## (e.g. relative humidity in [0, 1]) with shape1 = 0.6, shape2 = 0.8
#' ## and ParetoII autocorrelation structure.
#' z <- generateTS(margdist = "beta",
#'                 margarg = list(shape1 = .6,
#'                                shape2 = .8),
#'                 distbounds = c(0, 1),
#'                 acsvalue = acs(id = "paretoII",
#'                                t = 0:30,
#'                                scale = 1,
#'                                shape = .75),
#'                 n = 1000,
#'                 p = 20)
#'
#' ## see the results
#' plot(z)
#'
#' ## Case 4:
#' ## Same as Case 3 but providing specific autocorrelation values for the
#' ## first three lags (lag 1 to 3 equal to 0.9, 0.8, 0.7).
#' z <- generateTS(margdist = "beta",
#'                 margarg = list(shape1 = .6,
#'                                shape2 = .8),
#'                 distbounds = c(0, 1),
#'                 acsvalue = c(1, .9, .8, .7),
#'                 n = 1000,
#'                 p = NULL)
#'
#' ## see the results
#' plot(z)
#'
#' }
#'
generateTS <- function(n, margdist, margarg, p = NULL, p0 = 0, TSn = 1,
                       distbounds = c(-Inf, Inf), acsvalue = NULL) {

  pnts <- actpnts(margdist = margdist,
                  margarg = margarg,
                  p0 = p0,
                  distbounds = distbounds)

  par <- fitactf(pnts)

  out <- lapply(1:TSn, function(i) {
    ARp(margdist = margdist,
        margarg = margarg,
        n = n,
        p = p,
        p0 = p0,
        actfpara = par,
        acsvalue = acsvalue)
  })

  structure(.Data = out, class = "cosmosts")
}

#' Bulk time series generation
#'
#' Generates additional time series using parameters already fitted by
#' \code{\link{generateTS}}, avoiding recomputation of the ACTF.
#'
#' @details
#' After calling \code{\link{generateTS}}, use \code{regenerateTS} to generate
#' more time series with the same fitted parameters. This is faster than
#' re-running \code{\link{generateTS}} because the ACTF fitting step is skipped.
#'
#' @param ts a \code{cosmosts} object returned by \code{\link{generateTS}}
#' @param TSn number of time series to generate
#'
#' @return An object of class \code{'cosmosts'}: a list of \code{TSn} numeric
#'   vectors of the same length as those in \code{ts}.
#'
#' @seealso \code{\link{generateTS}}, \code{\link{ARp}}
#'
#' @export
#'
#' @examples
#'
#' library(CoSMoS)
#'
#' ## Fit once
#' x <- generateTS(margdist = "burrXII",
#'                 margarg = list(scale = 1,
#'                                shape1 = .75,
#'                                shape2 = .25),
#'                 acsvalue = acs(id = "weibull",
#'                                t = 0:30,
#'                                scale = 10,
#'                                shape = .75),
#'                 n = 1000, p = 30, p0 = .5, TSn = 3)
#'
#' ## Generate more realisations with the same parameters
#' r <- regenerateTS(x)
#'
#' plot(r)
#'
regenerateTS <- function(ts, TSn = 1) {

  X <- ts

  margdist <- attr(X[[1]], "margdist")
  margarg  <- attr(X[[1]], "margarg")

  n   <- length(X[[1]])
  p   <- length(attr(X[[1]], "ar_coef"))
  a   <- attr(X[[1]], "ar_coef")
  p0  <- attr(X[[1]], "p0")
  esd <- attr(X[[1]], "noise_sd")

  out <- lapply(1:TSn, function(i) {

    gn  <- c(rep(0, (2 * p)), rnorm(n, mean = 0, sd = esd))
    val <- c(AR1((2 * p), alpha = a[1]), rep(0, n))

    for (j in (p + 1):(n + (2 * p))) {
      val[j] <- sum(rev(a) * val[(j - p):(j - 1)]) + gn[j]
    }

    uval <- (pnorm(val[-1:-(2 * p)]) - p0) / (1 - p0)
    uval[uval < 0] <- 0

    do.call(paste0("q", margdist), args = c(list(p = uval), margarg))
  })

  structure(.Data = out, class = "cosmosts")
}
