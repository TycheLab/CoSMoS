#' Autoregressive model of first order
#'
#' Generates time series from an AR1 model
#'
#' @param n number of values
#' @param alpha lag-1 autocorrelation
#' @param mean mean
#' @param sd standard deviation
#'
#' @keywords internal
#' @export
#' @import stats ggplot2
#'
#' @examples
#'
#' library(CoSMoS)
#'
#' ## generate 500 values from an AR1 having lag-1 autocorrelation 0.8,
#' ## mean value equal to 0, and standard deviation equal to 1.
#' n <- 500
#'
#' ## generate white noise for comparsion
#' x <- rnorm(n)
#' ggplot() +
#'  geom_line(aes(x = 1:n,
#'                y = x)) +
#'    labs(x = '',
#'         y = 'value') +
#'    theme_classic()
#'
#' ## generete values using AR1
#' y <- AR1(n, .8)
#' ggplot() +
#'   geom_line(aes(x = 1:n,
#'                 y = y)) +
#'   labs(x = '',
#'        y = 'value') +
#'   theme_classic()
#'
AR1 <- function(n, alpha, mean = 0, sd = 1) {

  emean <- (1 - alpha) * mean ## Gaussian noise mean
  esd <- sqrt(1 - alpha ^ 2) * sd ## Gaussian noise sd

  val <- c(rnorm(n = 1,
                 mean = mean,
                 sd = sd)) ## values vector

  if(n != 1) {

    gn <- rnorm(n = n,
                mean = emean,
                sd = esd) ## Gaussain noise vector

    for (i in 2:n) { ## AR
      val[i] <- val[(i - 1)] * alpha + gn[i]
    }
  }

  return(val)
}
