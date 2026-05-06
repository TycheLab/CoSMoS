#' Autoregressive model of order 1
#'
#' Generates a Gaussian time series from a first-order autoregressive (AR(1))
#' model with specified lag-1 autocorrelation, mean, and standard deviation.
#'
#' @param n      Positive integer. Number of values to generate.
#' @param alpha  Numeric in \eqn{(-1, 1)}. Lag-1 autocorrelation coefficient.
#' @param mean   Numeric. Process mean. Default \code{0}.
#' @param sd     Positive numeric. Process standard deviation. Default \code{1}.
#'
#' @return A numeric vector of length \code{n}.
#'
#' @keywords internal
#'
AR1 <- function(n, alpha, mean = 0, sd = 1) {

  ## --- input validation -------------------------------------------------------
  n <- suppressWarnings(as.integer(n))
  if (is.na(n) || n < 1L) {
    stop("'n' must be a positive integer.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || abs(alpha) >= 1) {
    stop("'alpha' must be a single numeric value in (-1, 1).")
  }
  if (!is.numeric(sd) || length(sd) != 1L || sd <= 0) {
    stop("'sd' must be a single positive numeric value.")
  }

  ## --- innovation parameters --------------------------------------------------
  noise_mean <- (1 - alpha) * mean
  noise_sd   <- sqrt(1 - alpha^2) * sd

  ## --- initialise output ------------------------------------------------------
  val <- rnorm(n = 1L, mean = mean, sd = sd)

  ## --- AR(1) recursion --------------------------------------------------------
  if (n > 1L) {
    noise <- rnorm(n = n, mean = noise_mean, sd = noise_sd)
    for (i in seq(2L, n)) {
      val[i] <- alpha * val[i - 1L] + noise[i]
    }
  }

  val
}


#' Yule-Walker solver
#'
#' Solves the Yule-Walker equations to compute AR(\eqn{p}) coefficients from
#' an autocorrelation structure (ACS) vector. Used internally by
#' \code{\link{ARp}} and \code{\link{seasonalAR}}.
#'
#' @param ACS Numeric vector. ACS values starting from lag 0
#'   (\code{ACS[1] = 1}). The AR order is \code{length(ACS) - 1}.
#'
#' @return Numeric matrix of AR coefficients, shape \code{(p, 1)}.
#'
#' @keywords internal
#'
#' @seealso \code{\link{ARp}}, \code{\link{seasonalAR}}
#'
YW <- function(ACS) {

  p   <- length(ACS) - 1L
  idx <- seq_len(p)

  ## symmetric Toeplitz covariance matrix: P[i,j] = ACS[|i-j| + 1]
  P   <- outer(idx, idx, function(i, j) ACS[abs(i - j) + 1L])
  rho <- matrix(ACS[seq(2L, p + 1L)], nrow = p, ncol = 1L)

  solve(P, rho)
}


#' Autoregressive model of order \emph{p}
#'
#' Generates a time series from an AR(\eqn{p}) model. The parent Gaussian
#' process is constructed to match a target autocorrelation structure via the
#' Yule-Walker equations, then transformed to a target marginal distribution
#' and optional intermittency (probability of zeros).
#'
#' @param margdist  Character. Name of the target marginal distribution
#'   (e.g. \code{"ggamma"}, \code{"paretoII"}). Must have a quantile function
#'   \code{q<margdist>} available.
#' @param margarg   Named list. Parameters of the marginal distribution passed
#'   to \code{q<margdist>}.
#' @param acsvalue  Numeric vector. Target autocorrelation structure starting
#'   from lag 0 (i.e. \code{acsvalue[1] = 1}).
#' @param actfpara  List returned by \code{\link{fitactf}}, containing the
#'   fitted ACTF coefficients in \code{actfpara$actfcoef}.
#' @param n         Positive integer. Length of the generated time series.
#' @param p         Positive integer or \code{NULL}. AR model order. When
#'   \code{NULL} (default), the order is chosen automatically as the number of
#'   lags where the transformed ACS exceeds 0.01, capped at 1000.
#' @param p0        Numeric in \eqn{[0, 1)}. Probability of zero values
#'   (intermittency). Default \code{0}.
#'
#' @return A numeric vector of length \code{n} with the following attributes:
#' \describe{
#'   \item{margdist}{Target marginal distribution name.}
#'   \item{margarg}{Target marginal distribution parameters.}
#'   \item{acsvalue}{Original (untransformed) target ACS.}
#'   \item{p0}{Probability of zeros.}
#'   \item{ar_coef}{AR(\eqn{p}) coefficients (Yule-Walker solution).}
#'   \item{noise_sd}{Standard deviation of the Gaussian innovations.}
#'   \item{gaussian}{The underlying Gaussian process (length \code{n}).}
#'   \item{transformed_acs}{ACTF-transformed ACS used for the AR model.}
#' }
#'
#' @keywords internal
#'
#' @importFrom matrixcalc is.positive.definite
#'
#' @seealso \code{\link{generateTS}}, \code{\link{actpnts}},
#'   \code{\link{fitactf}}, \code{\link{acs}}
#'
ARp <- function(margdist, margarg, acsvalue, actfpara, n, p = NULL, p0 = 0) {

  ## --- input validation -------------------------------------------------------
  if (!is.character(margdist) || length(margdist) != 1L) {
    stop("'margdist' must be a single character string.")
  }
  if (!is.list(margarg)) {
    stop("'margarg' must be a named list of distribution parameters.")
  }
  if (!is.numeric(acsvalue) || acsvalue[1] != 1) {
    stop("'acsvalue' must be a numeric vector starting with 1 (lag-0 value).")
  }
  n <- suppressWarnings(as.integer(n))
  if (is.na(n) || n < 1L) {
    stop("'n' must be a positive integer.")
  }
  if (!is.numeric(p0) || length(p0) != 1L || p0 < 0 || p0 >= 1) {
    stop("'p0' must be a single numeric value in [0, 1).")
  }

  ## --- transform ACS via ACTF -------------------------------------------------
  b <- actfpara$actfcoef[1]
  c <- actfpara$actfcoef[2]
  transformed_acs <- actf(acsvalue, b = b, c = c)

  ## --- determine AR order p ---------------------------------------------------
  max_order <- 1000L

  if (is.null(p) || p > max_order) {
    n_sig <- length(transformed_acs[transformed_acs > 0.01]) - 1L
    p     <- min(n_sig, max_order)
    message(sprintf('AR order "p" automatically set to %d.', p))
  }

  if (length(transformed_acs) - 1L < p) {
    stop(sprintf(
      "ACS is too short for order p = %d. Supply a longer 'acsvalue' or set p = NULL.",
      p
    ))
  }

  ## --- solve Yule-Walker via shared YW() helper -------------------------------
  rho     <- matrix(transformed_acs[seq(2L, p + 1L)], nrow = p, ncol = 1L)

  ## build Toeplitz covariance matrix for positive-definiteness check
  idx <- seq_len(p)
  P   <- outer(idx, idx, function(i, j) transformed_acs[abs(i - j) + 1L])

  if (!is.positive.definite(P, tol = 1e-8)) {
    stop(paste0(
      "Covariance matrix is not positive definite. ",
      "Try a different ACS (e.g. Pareto II, Weibull, or fGn)."
    ))
  }

  ar_coef  <- YW(transformed_acs[seq_len(p + 1L)])  # uses shared helper
  noise_sd <- sqrt(1 - sum(rho * ar_coef))

  ## --- generate parent Gaussian process ---------------------------------------
  gauss <- AR1(n = p, alpha = rho[1])

  innovations <- rnorm(n = p + n, mean = 0, sd = noise_sd)
  ar_coef_rev <- rev(ar_coef)

  for (i in seq(p + 1L, n + p)) {
    gauss[i] <- sum(gauss[(i - p):(i - 1L)] * ar_coef_rev) + innovations[i]
  }

  ## --- transform to target marginal -------------------------------------------
  gauss_out      <- gauss[-seq_len(p)]
  uval           <- (pnorm(gauss_out) - p0) / (1 - p0)
  uval[uval < 0] <- 0

  out <- do.call(paste0("q", margdist), args = c(list(p = uval), margarg))

  ## --- attach metadata as attributes ------------------------------------------
  structure(
    .Data           = out,
    margdist        = margdist,
    margarg         = margarg,
    acsvalue        = acsvalue,
    p0              = p0,
    ar_coef         = ar_coef,
    noise_sd        = noise_sd,
    gaussian        = gauss_out,
    transformed_acs = transformed_acs
  )
}
