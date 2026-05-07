# actf.R
# Merged: actf, actfInv, acti, actpnts, actpntsB6, fitactf
# Group C — style pass applied; C++ fast path added to actpnts()

#' ACTF auto-correlation transformation function
#'
#' Provides transformation for continuous distributions, based on two parameters.
#'
#' @param rhox marginal correlation value
#' @param b 1st parameter
#' @param c 2nd parameter
#'
#' @return numeric vector of transformed Gaussian correlations
#' @seealso \code{\link{actpnts}}, \code{\link{fitactf}}
#' @export
#' @keywords internal
#'
#' @examples
#' actf(.4, 1, 0)
#'
actf <- function(rhox, b, c) {
  ((1 + b * rhox) ^ (1 - c) - 1) / ((1 + b) ^ (1 - c) - 1)
}

#' ACTF auto-correlation transformation function for discrete distributions
#'
#' Provides transformation for discrete distributions, based on two parameters.
#'
#' @inheritParams actf
#'
#' @return numeric vector of transformed Gaussian correlations
#' @seealso \code{\link{actpnts}}, \code{\link{fitactf}}
#' @export
#' @keywords internal
#'
#' @examples
#' actfdiscrete(.4, .2, 1)
#'
actfdiscrete <- function(rhox, b, c) {
  1 - (1 - rhox ^ b) ^ c
}

#' Inverse ACTF
#'
#' Inverse of the autocorrelation transformation function.
#'
#' @inheritParams actf
#' @param rhoz Gaussian correlation value
#'
#' @return numeric vector of back-transformed marginal correlations
#' @seealso \code{\link{actf}}
#' @keywords internal
#' @export
#'
#' @examples
#' actfInv(.4, 1, 0)
#'
actfInv <- function(rhoz, b, c) {
  ((rhoz * ((1 + b) ^ (1 - c) - 1) + 1) ^ (1 / (1 - c)) - 1) / b
}

#' ACTI — autocorrelation transformation integral function
#'
#' Expression supplied to the double integral inside \code{\link{actpnts}}.
#'
#' @param x x-plane value
#' @param y y-plane value
#' @param dist distribution name
#' @param distarg a list of distribution arguments
#' @param rhoz Gaussian correlation
#' @param p0 probability of zero values
#'
#' @return numeric scalar — integrand value at (x, y)
#' @seealso \code{\link{actpnts}}
#' @export
#' @import stats
#' @keywords internal
#'
#' @examples
#' acti(1, -1, "norm", list(mean = 0, sd = 1), .3, 0)
#'
acti <- function(x, y, dist, distarg, rhoz, p0) {

  do.call(what = paste0("q", dist),
          args = c(list(p = (erfc(-x / sqrt(2)) / 2 - p0) / (1 - p0)),
                   distarg)) *
    do.call(what = paste0("q", dist),
            args = c(list(p = (erfc(-y / sqrt(2)) / 2 - p0) / (1 - p0)),
                     distarg)) *
    exp((x ^ 2 + y ^ 2 - 2 * x * y * rhoz) /
          (2 * (-1 + rhoz ^ 2))) / (2 * pi * sqrt(1 - rhoz ^ 2))
}

#' AutoCorrelation Transformed Points
#'
#' Evaluates the (rho_x, rho_z) mapping between the target marginal
#' autocorrelation and the underlying Gaussian autocorrelation, using a
#' double numerical integral.
#'
#' @details
#' When the package is compiled with Rcpp support (i.e., \code{actpnts_cpp}
#' is available), the double integral is evaluated in C++ via the Cubature
#' algorithm, which is substantially faster than the nested base-R
#' \code{integrate()} fallback. The C++ path supports the following
#' distributions natively: \code{ggamma}, \code{paretoII}, \code{burrXII},
#' \code{burrIII}, \code{gev}, \code{norm}, \code{beta}, \code{gamma},
#' \code{exp}, \code{weibull}, \code{lnorm}, \code{unif}. Any other
#' distribution falls back to the R quantile function automatically, so
#' correctness is always preserved.
#'
#' @param margdist target marginal distribution
#' @param margarg list of marginal distribution arguments
#' @param p0 probability zero
#' @inheritParams moments
#'
#' @return A data frame with columns \code{rhoz} (Gaussian correlations) and
#'   \code{rhox} (corresponding target marginal correlations).
#'
#' @seealso \code{\link{fitactf}}, \code{\link{acti}}, \code{\link{generateTS}}
#' @references Papalexiou, S.M. (2018). Unified theory for stochastic modelling
#'   of hydroclimatic processes: Preserving marginal distributions, correlation
#'   structures, and intermittency. Advances in Water Resources, 115, 234-252,
#'   \doi{10.1016/j.advwatres.2018.02.013}
#'
#' @export
#' @import stats ggplot2
#'
#' @examples
#'
#' library(CoSMoS)
#'
#' ## Pareto type II marginal
#' x <- actpnts(margdist = "paretoII",
#'              margarg  = list(scale = 1, shape = .3),
#'              p0 = 0)
#' x
#'
actpnts <- function(margdist, margarg, p0 = 0, distbounds = c(-Inf, Inf)) {

  rhoz_vals <- c(seq(from = 0.1, to = 0.9, by = 0.1), 0.95)

  .min <- ifelse(test = p0 == 0,
                 yes  = -7.5,
                 no   = -sqrt(2) * inv.erfc(2 * p0))
  .max <- 7.5

  m <- moments(dist      = margdist,
               distarg   = margarg,
               raw       = FALSE,
               central   = TRUE,
               coef      = FALSE,
               distbounds = distbounds,
               p0        = p0,
               order     = 1:2)

  mu1 <- m$mu["mu1"]
  mu2 <- m$mu["mu2"]

  ## --- Fast C++ path (requires compiled actpnts_cpp) ---
  cpp_available <- exists("actpnts_cpp",
                          envir  = asNamespace("CoSMoS"),
                          mode   = "function",
                          inherits = FALSE)

  if (cpp_available) {

    ## Convert margarg list to named numeric vector for C++
    margarg_vec <- unlist(margarg)

    rho <- actpnts_cpp(
      margdist  = margdist,
      margarg   = margarg_vec,
      p0        = p0,
      lower     = .min,
      upper     = .max,
      mu1       = mu1,
      mu2       = mu2,
      rhoz_vals = rhoz_vals
    )

    return(structure(.Data = rho))
  }

  ## --- Pure R fallback (original implementation) ---
  rho <- data.frame(rhoz = rhoz_vals, rhox = 0)

  for (i in seq_len(nrow(rho))) {

    temp <- integrate(
      f = function(y) {
        sapply(y, function(y) {
          integrate(
            f           = function(x) {
              acti(x        = x,
                   y        = y,
                   rhoz     = rho[i, "rhoz"],
                   p0       = p0,
                   dist     = margdist,
                   distarg  = margarg)
            },
            lower        = .min,
            upper        = .max,
            subdivisions = 1.0e4,
            rel.tol      = 1.0e-5
          )$value
        })
      },
      lower        = .min,
      upper        = .max,
      subdivisions = 1.0e4,
      rel.tol      = 1.0e-5
    )$value

    rho[i, "rhox"] <- (temp - mu1 ^ 2) / mu2
  }

  structure(.Data = rho)
}

#' AutoCorrelation Transformed Points for Bardossy dependence structure
#'
#' Evaluates the (rho_x, rho_z) mapping using Monte Carlo integration for the
#' Bardossy copula dependence structure.
#'
#' @param margdist target marginal distribution
#' @param margarg list of marginal distribution arguments
#' @param m mean of the parent Gaussian processes controlling asymmetry
#' @param p0 probability of zero values
#'
#' @return A data frame with columns \code{rhoz} and \code{rhox}.
#'
#' @seealso \code{\link{actpnts}}, \code{\link{generateMTSFast}}
#' @references Bardossy, A. (2006). Copula-based geostatistical models for
#'   groundwater quality parameters. Water Resources Research, 42(11),
#'   \doi{10.1029/2005WR004754}
#'
#' @export
#' @import stats mvtnorm ggplot2
#'
#' @examples
#'
#' library(CoSMoS)
#'
#' x <- actpntsB6(margdist = "paretoII",
#'                margarg  = list(scale = 1, shape = .3),
#'                m = 1,
#'                p0 = 0)
#' x
#'
actpntsB6 <- function(margdist, margarg, m, p0 = 0) {

  rho <- data.frame(rhoz = c(0.05, seq(from = 0.1, to = 0.9, by = 0.1), 0.95),
                    rhox = 0)

  set.seed(666)

  for (i in seq_len(nrow(rho))) {

    rhoz     <- rho[i, "rhoz"]
    mc.size  <- 2e5
    rnd      <- rmvnorm(mc.size,
                        mean  = rep(m, 2),
                        sigma = matrix(c(1, rhoz, rhoz, 1), 2, 2))

    u <- apply(rnd ^ 2, 2,
               function(x) jitter(rank(x) / (mc.size + 1),
                                  amount = 1 / (mc.size + 1)))

    x        <- do.call(what = paste0("q", margdist),
                        args = c(list(p = (u - p0) / (1 - p0)), margarg))
    x[u <= p0] <- 0

    rho[i, "rhox"] <- (mean(x[, 1] * x[, 2]) -
                         mean(x[, 1]) * mean(x[, 2])) /
      (sd(x[, 1]) * sd(x[, 2]))
  }

  structure(.Data = rho)
}

#' Fit the AutoCorrelation Transformation Function
#'
#' Fits the ACTF to the estimated (rho_x, rho_z) points using \code{nls}.
#'
#' @param actpnts estimated ACT points (output of \code{\link{actpnts}})
#' @param discrete logical — is the marginal distribution discrete?
#'
#' @return An object of class \code{"acti"} with components:
#'   \describe{
#'     \item{actfcoef}{fitted ACTF coefficients \code{b} and \code{c}}
#'     \item{actfpoints}{the input ACT points data frame}
#'   }
#'
#' @seealso \code{\link{actpnts}}, \code{\link{actf}}
#' @export
#'
#' @examples
#'
#' library(CoSMoS)
#'
#' p   <- actpnts(margdist = "paretoII",
#'                margarg  = list(scale = 1, shape = .3),
#'                p0 = 0)
#' fit <- fitactf(p)
#' plot(fit)
#'
fitactf <- function(actpnts, discrete = FALSE) {

  suppressWarnings(
    fit <- nls(
      rhoz ~ do.call(
        what = ifelse(discrete, "actfdiscrete", "actf"),
        args = list(rhox, b, c)
      ),
      data      = list(rhoz = actpnts$rhoz, rhox = actpnts$rhox),
      start     = c(b = 1, c = 0),
      lower     = c(b = .001, c = 0),
      algorithm = "port",
      control   = nls.control(maxiter = 50000, warnOnly = TRUE)
    )
  )

  structure(
    .Data = list(actfcoef   = coefficients(fit),
                 actfpoints = actpnts),
    class = "acti"
  )
}
