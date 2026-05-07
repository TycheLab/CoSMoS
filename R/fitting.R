# --- Error norms (internal) ---------------------------------------------------

#' Ratio mean squared error
#'
#' @param x numeric vector of observed values
#' @param y numeric vector of reference values
#'
#' @return scalar numeric
#'
#' @seealso \code{\link{MSE}}
#' @keywords internal
rMSE <- function(x, y) {
  sum((x / y - 1)^2) / length(y)
}

#' Mean squared error
#'
#' @param x numeric vector of observed values
#' @param y numeric vector of reference values
#'
#' @return scalar numeric
#'
#' @seealso \code{\link{rMSE}}
#' @keywords internal
MSE <- function(x, y) {
  sum((x - y)^2) / length(y)
}


# --- Empirical CDF (internal) -------------------------------------------------

#' Empirical cumulative distribution function
#'
#' Computes the ECDF using the Weibull plotting position.
#'
#' @param x numeric vector of values
#'
#' @return a \code{data.table} with columns \code{p} (plotting positions) and
#'   \code{value} (sorted unique values), keyed on \code{value}
#'
#' @seealso \code{\link{fitDist}}
#' @keywords internal
ECDF <- function(x) {

  st <- sort(x)

  aux <- data.table(p     = rank(x = st, ties.method = "first") / (length(st) + 1),
                    value = st,
                    key   = "value")

  J <- value <- NULL

  structure(.Data = aux[J(unique(value)), mult = "last"])
}


# --- Fitting norm (internal) --------------------------------------------------

#' Norm minimised during distribution fitting
#'
#' Computes one of four fitting norms used as the objective function in
#' \code{\link{fitDist}}.
#'
#' @param par numeric vector of distribution parameter values
#' @param val numeric vector of empirical data
#' @param dist character; name of the distribution to fit (e.g. \code{"norm"},
#'   \code{"ggamma"})
#' @param norm character; norm identifier — one of \code{"N1"}, \code{"N2"},
#'   \code{"N3"}, \code{"N4"}
#' @param n.points integer; number of ECDF points used in the norm computation
#'
#' @return scalar numeric; the norm value (returns \code{10e9} if non-finite)
#'
#' @seealso \code{\link{fitDist}}
#' @keywords internal
N <- function(par, val, dist, norm, n.points) {

  edf <- ECDF(val)
  aux <- nrow(edf)

  if (aux != n.points) {
    if (aux < n.points) n.points <- aux
    edf <- edf[seq(from = 1, to = aux, length.out = n.points), ]
  }

  a   <- getDistArg(dist)
  arg <- as.list(setNames(par, a))

  if (norm %in% c("N1", "N2")) {
    FN <- edf$value
    F  <- do.call(paste0("q", dist), args = c(list(edf$p), arg))
    out <- if (norm == "N1") rMSE(F, FN) else MSE(F, FN)
  }

  if (norm %in% c("N3", "N4")) {
    Xi <- edf$p
    Xu <- do.call(paste0("p", dist), args = c(list(edf$value), arg))
    out <- if (norm == "N3") rMSE(Xu, Xi) else MSE(Xu, Xi)
  }

  if (!is.finite(out)) out <- 10e9

  out
}


# --- ACS fitting objective (internal) -----------------------------------------

#' Auxiliary objective function for \code{fitACS}
#'
#' Computes the mean squared error between a theoretical autocorrelation
#' structure and the empirical ACS; passed as \code{fn} to \code{optim} inside
#' \code{\link{fitACS}}.
#'
#' @param par numeric vector of ACS parameter values
#' @param id character; ACS identifier (e.g. \code{"weibull"})
#' @param eACS numeric vector of empirical ACS values
#' @param error character; error criterion — currently only \code{"MSE"}
#'
#' @return scalar numeric; the MSE between theoretical and empirical ACS
#'
#' @seealso \code{\link{fitACS}}
#' @keywords internal
optimACS <- function(par, id, eACS, error = "MSE") {

  arg <- as.list(setNames(par, getACSArg(id)))
  ACS <- do.call(acs, c(list(id = id, t = 0:(length(eACS) - 1)), arg))
  do.call(error, list(x = ACS, y = eACS))
}


# --- Public fitting functions -------------------------------------------------

#' Distribution fitting
#'
#' Fits a parametric distribution to data using the Nelder-Mead simplex
#' algorithm to minimise one of four fitting norms.
#'
#' @param data numeric vector of values to fit
#' @param dist character; distribution name (e.g. \code{"norm"}, \code{"ggamma"})
#' @param n.points integer; number of ECDF points used in norm computation
#' @param norm character; norm identifier — one of \code{"N1"} (ratio RMSE on
#'   quantiles), \code{"N2"} (MSE on quantiles), \code{"N3"} (ratio RMSE on
#'   probabilities), \code{"N4"} (MSE on probabilities)
#' @param constrain logical; if \code{TRUE}, constrains \code{shape2} parameters
#'   to \code{(0, 0.48)} to enforce finite upper tails
#' @param opts list of \code{nloptr} minimisation options
#'
#' @return An object of class \code{"fitDist"}: a named list of fitted
#'   distribution parameters with attributes \code{dist}, \code{edf} (empirical
#'   CDF), and \code{nfo} (full \code{nloptr} output).
#'
#' @seealso \code{\link{fitACS}}, \code{\link{plot.fitDist}}
#'
#' @export
#' @import nloptr
#'
#' @examples
#'
#' x <- fitDist(rnorm(1000), "norm", 30, "N1", FALSE)
#' x
#'
fitDist <- function(data,
                    dist,
                    n.points,
                    norm,
                    constrain,
                    opts = list("algorithm" = "NLOPT_LN_NELDERMEAD",
                                "xtol_rel"  = 1.0e-8,
                                "maxeval"   = 1.0e4)) {

  a   <- getDistArg(dist)
  lwr <- NULL

  if (constrain) {
    if (any(a == "shape2")) {
      start <- c(rep(1, length(a) - 1), .1)
      upr   <- c(rep(Inf, length(a) - 1), .48)
    } else {
      start <- rep(1, length(a))
      upr   <- rep(Inf, length(a))
    }
  } else {
    if (any(a == "shape2")) {
      start <- c(rep(1, length(a) - 1), .1)
      upr   <- rep(Inf, length(a))
    } else {
      start <- rep(1, length(a))
      upr   <- rep(Inf, length(a))
    }
  }

  fit <- nloptr(x0     = start,
                eval_f = N,
                lb     = lwr,
                ub     = upr,
                val    = data,
                dist   = dist,
                n.points = n.points,
                norm   = norm,
                opts   = opts)

  out <- as.list(setNames(fit$solution, a))

  structure(.Data = out,
            dist  = dist,
            edf   = ECDF(data),
            nfo   = fit,
            class = "fitDist")
}


#' Autocorrelation structure fitting
#'
#' Fits a parametric autocorrelation structure (ACS) to empirical ACF values
#' using Nelder-Mead optimisation with MSE criterion.
#'
#' @param acf numeric vector of autocorrelation function values from lag 0
#' @param ID character; ACS identifier (e.g. \code{"weibull"}, \code{"paretoII"})
#' @param start numeric vector of starting parameter values; if \code{NULL},
#'   all parameters start at 1
#' @param lag integer; number of lags to use; if \code{NULL}, lags up to the
#'   first value \eqn{\le 0.01} are used (or all lags if none drops below 0.01)
#'
#' @return An object of class \code{"fitACS"}: a named list of fitted ACS
#'   parameters with attributes \code{ID} and \code{eACS} (empirical ACS used
#'   for fitting).
#'
#' @seealso \code{\link{fitDist}}, \code{\link{plot.fitACS}}, \code{\link{acs}}
#'
#' @export
#' @import stats
#'
#' @examples
#'
#' x <- arima.sim(model = list(ar = 0.8), n = 1000)
#'
#' acsfit <- fitACS(acf(x, plot = FALSE)$acf, "weibull", c(1, 1))
#'
fitACS <- function(acf, ID, start = NULL, lag = NULL) {

  if (is.null(lag)) {
    if (length(which(acf <= .01)) == 0) {
      eACS <- acf
    } else {
      eACS <- acf[1:which(acf <= .01)[1]]
    }
  } else {
    eACS <- acf[1:(lag + 1)]
  }

  if (is.null(start)) start <- rep(1, length(getACSArg(ID)))

  par <- optim(par = start,
               fn  = optimACS,
               id  = ID,
               eACS = eACS,
               error = "MSE")$par

  out <- as.list(setNames(par, getACSArg(ID)))

  structure(.Data = out,
            ID   = ID,
            eACS = eACS,
            class = "fitACS")
}
