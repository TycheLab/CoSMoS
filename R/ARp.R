#' Autoregressive model of order \emph{p}
#'
#' Generates time series from an Autoregressive model of order \emph{p}
#'
#' @inheritParams actpnts
#' @param n number of values
#' @param p integer - model order (if NULL - limits maximum model order according to auto-correlation structure values)
#' @param actfpara auto-correlation structure transformation parameters
#' @param acsvalue target auto-correlation structure (from lag 0)
#' @param p0 probability zero
#'
#' @keywords internal
#' @import stats ggplot2
#' @export
#'
#' @examples
#'
#' library(CoSMoS)
#'
#' ## choose the marginal distribution as Pareto type II with corresponding parameters
#' dist <- 'paretoII'
#' distarg <- list(scale = 1, shape = .3)
#' p0 <- .5
#'
#' ## estimate rho 'x' and 'z' points using ACTI
#' pnts <- actpnts(margdist = dist, margarg = distarg, p0 = p0)
#'
#' ## fit ACTF
#' fit <- fitactf(pnts)
#'
#' ## define target auto-correlation structure and model order
#' order <- 1000
#' acsvalue <- acs(id = 'weibull', t = 0:order, scale = 10, shape = .75)
#'
#' ## limit ACS lag (recomended)
#' system.time(val <- ARp(margdist = dist,
#'                        margarg = distarg,
#'                        acsvalue = acsvalue,
#'                        actfpara = fit,
#'                        n = 5000,
#'                        p0 = p0))
#' \donttest{
#' ## order w/o limit
#' system.time(val <- ARp(margdist = dist,
#'                        margarg = distarg,
#'                        acsvalue = acsvalue,
#'                        actfpara = fit,
#'                        n = 5000,
#'                        p = order,
#'                        p0 = p0))
#' }
#'
#' ## see the result
#' ggplot() +
#'   geom_col(aes(x = seq_along(val),
#'                y = val)) +
#'   labs(x = '',
#'        y = 'value') +
#'   theme_classic()
#'
ARp <- function(margdist, margarg, acsvalue, actfpara, n, p = NULL, p0 = 0) {

  transacsvalue <- actf(acsvalue,
                        b = actfpara$actfcoef[1],
                        c = actfpara$actfcoef[2])

  if (any(c(is.null(p), p > 1000))) { ## limit ARp order by the ACS values or by hard cap 1000

    temp <- length(transacsvalue[transacsvalue > .01]) - 1

    p <- ifelse(temp > 1000, 1000, temp)
    message(paste('Order "p" limited to ', p))
  }

  if(length(transacsvalue) - 1 < p) {

    stop(paste0('Please supply ACS of length up to order p = ', p, ', or leave argument p = NULL'))
  }

  P <- matrix(NA, p, p) ## cov matrix generation

  for (i in 1:p) {
    P[i, i:p] <- transacsvalue[1:(p - i + 1)]
    P[i, 1:i] <- transacsvalue[i:1]
  }

  rho <- matrix(transacsvalue[2:(p + 1)], p, 1)

  a <- solve(P, rho) ## Yule-Walker

  esd <- sqrt(1 - sum(rho*a)) ## gaussian noise standard deviation

  val <- AR1(n = p,
             alpha = rho[1]) ## values vector (first values are generated using AR1 to ensure ACS)

  gn <- rnorm(n = (p + n),
              mean = 0,
              sd = esd) ## Gaussian noise generation

  a.rev <- rev(a)

  for (i in (p + 1):(n + p)) { ## AR
    val[i] <- sum(val[(i - p):(i - 1)]*a.rev) + gn[i]
  }

  uval <- (pnorm(val[-1:-p]) - p0)/(1 - p0) ## p0 + gaussian probabilities calculation
  uval[uval < 0] <- 0

  out <- do.call(paste0('q', margdist), args = c(list(p = uval), margarg))

  structure(.Data = out,
            margdist = margdist,
            margarg = margarg,
            acsvalue = acsvalue,
            p0 = p0,
            a = a,
            esd = esd,
            gaussian = val[-1:-p],
            trACS = transacsvalue)
}
