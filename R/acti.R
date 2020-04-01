#' ACTI - autocorrelation transformation integral function
#'
#' Expression supplied to double integral
#'
#' @param x x-plain value
#' @param y y-plain value
#' @param dist distribution
#' @param distarg a list of distribution arguments
#' @param rhoz gaussian correlation
#' @param p0 probability od zero values
#'
#' @export
#' @import stats
#' @keywords internal
#'
#' @examples
#'
#' acti(1, -1, 'norm', list(), .3, 0)
#'
acti <- function(x, y, dist, distarg, rhoz, p0) {

  do.call(what = paste0('q', dist),
          args = c(list(p = (erfc(-x / sqrt(2)) / 2 - p0) / (1 - p0)),
                   distarg)) *
    do.call(what = paste0('q', dist),
            args = c(list(p = (erfc(-y / sqrt(2)) / 2 - p0) / (1 - p0)),
                     distarg)) *
    exp((x ^ 2 + y ^ 2 - 2 * x * y * rhoz) /
          (2 * (-1 + rhoz ^ 2))) / (2 * pi * sqrt(1 - rhoz ^ 2))
}

#' AutoCorrelation Transformed Points
#'
#' Transforms a gaussian process in order to match a target marginal lowers its autocorrelation values. The actpnts evaluates the corresponding autocorrelations for the given target marginal for a set of gaussian correlations, i.e., it returns  (\eqn{\rho_x , \rho_z}) points where \eqn{\rho_x and \rho_z} represent, respectively, the autocorrelations of the target and gaussian process.
#'
#' @param margdist target marginal distribution
#' @param margarg list of marginal distribution arguments
#' @param p0 probability zero
#' @inheritParams moments
#'
#' @export
#' @import stats ggplot2
#'
#' @examples
#'
#' library(CoSMoS)
#'
#' ## here we target to a process that has the Pareto type II marginal distribution
#' ## with scale parameter 1 and shape parameter 0.3
#' ## (note that all parameters have to be named)
#' dist <- 'paretoII'
#' distarg <- list(scale = 1, shape = .3)
#'
#' x <- actpnts(margdist = dist, margarg = distarg, p0 = 0)
#' x
#'
#' ## you can see the points by using
#' ggplot(x,
#'        aes(x = rhox,
#'            y = rhoz)) +
#'   geom_point(colour = 'royalblue4', size = 2.5) +
#'   geom_abline(lty = 5) +
#'   labs(x = bquote(Autocorrelation ~ rho[x]),
#'        y = bquote(Gaussian ~ rho[z])) +
#'   scale_x_continuous(limits = c(0, 1)) +
#'   scale_y_continuous(limits = c(0, 1)) +
#'   theme_classic()
#'
actpnts <- function(margdist, margarg, p0 = 0, distbounds = c(-Inf, Inf)) {

  rho <- data.frame(rhoz = c(seq(from = 0.1,
                                 to = .9,
                                 by = .1),
                             .95),
                    rhox = 0) ## create data frame of marginal ACS values

  .min <- ifelse(test = p0 == 0,
                 yes = -8,
                 no = -sqrt(2) * inv.erfc(2 * p0)) ## double integral lower bound
  .max <- 8 ## double integral upper bound

  m <- moments(dist = margdist, ## moment calculation
               distarg = margarg,
               raw = F,
               central = T,
               coef = F,
               distbounds = distbounds,
               p0 = p0,
               order = 1:2)

  for (i in 1:dim(rho)[1]) {

    # temp <- integral2(acti, ## ACTI calculation using pracma
    #                   ymin = .min,
    #                   ymax = .max,
    #                   xmin = .min,
    #                   xmax = .max,
    #                   rhoz = rho[i, 'rhoz'],
    #                   p0 = p0,
    #                   dist = margdist,
    #                   distarg = margarg)$Q

    temp <- integrate( ## ACTI using base
      f = function(y) {
        sapply(y, function(y) {
          integrate(
            f = function(x) {
              acti(x = x,
                   y = y,
                   rhoz = rho[i, 'rhoz'],
                   p0 = p0,
                   dist = margdist,
                   distarg = margarg)
            },
            lower = .min,
            upper = .max,
            subdivisions = 1.0e3,
            rel.tol = 1.0e-5#,
            # stop.on.error = FALSE
          )$value
        })
      },
      lower = .min,
      upper = .max,
      subdivisions = 1.0e3,
      rel.tol = 1.0e-5#,
      # stop.on.error = FALSE
    )$value

    rho[i, 'rhox'] <- (temp - m[[1]]['mu1'] ^ 2) / (m[[1]]['mu2'])
  }

  structure(.Data = rho)
}

#' Fit the AutoCorrelation Transformation Function
#'
#' Fits the ACTF (Autocorrelation Transformation Function) to the estimated points (\eqn{\rho_x, \rho_z}) using nls
#'
#' @param actpnts estimated ACT points
#' @param discrete logical - is the marginal distribution discrete?
#'
#' @export
#'
#' @examples
#'
#' library(CoSMoS)
#'
#' ## choose the marginal distribution as Pareto type II with corresponding parameters
#' dist <- 'paretoII'
#' distarg <- list(scale = 1, shape = .3)
#'
#' ## estimate rho 'x' and 'z' points using ACTI
#' p <- actpnts(margdist = dist, margarg = distarg, p0 = 0)
#'
#' ## fit ACTF
#' fit <- fitactf(p)
#'
#' ## plot the result
#' plot(fit)
#'
fitactf <- function(actpnts, discrete = FALSE) {

  suppressWarnings(fit <- nls(rhoz ~ do.call(what = ifelse(test = discrete, ## actf fit using nls
                                                           yes = 'actfdiscrete',
                                                           no = 'actf'),
                                             args = list(rhox, b, c)),
                              data = list(rhoz = actpnts$rhoz,
                                          rhox = actpnts$rhox),
                              start = c(b = 1, c = 0),
                              lower = c(b = .001, c = 0),
                              algorithm = 'port',
                              control = nls.control(maxiter = 50000,
                                                    warnOnly = TRUE)))

  structure(.Data = list(actfcoef = coefficients(fit),
                         actfpoints = actpnts),
            class = 'acti')
}

