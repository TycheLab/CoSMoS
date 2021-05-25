#' AutoCorrelation Transformation Function visualisation
#'
#' Visualizes the autocorrelation tranformation integral
#' (there are two possible methods for plotting - base graphics and ggplot2 package).
#'
#' @param x \code{\link{fitactf}} result object
#' @param ... other arguments
#'
#' @export
#' @method plot acti
#' @import ggplot2 graphics utils methods
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
#' ## plot the results
#' plot(fit)
#' plot(fit, main = 'Pareto type II distribution \nautocorrelation tranformation')
#'
plot.acti <- function(x, ...) {

  args <- list(...)

  main <- ifelse(!is.null(args[['main']]), args[['main']], '')

  acs <- x

  temp <- seq(0, 1, .001)

  dta <- data.frame(x = temp,
                    y = actf(temp, acs$actfcoef[1], acs$actfcoef[2]))

  p <- ggplot() +
    geom_line(aes(x = dta$x,
                  y = dta$y),
              colour = 'steelblue4',
              lwd = 1.5) +
    geom_point(aes(x = acs$actfpoints$rhox, y = acs$actfpoints$rhoz),
               colour = 'grey35',
               size = 3.5) +
    geom_abline(lty = 5) +
    scale_x_continuous(limits = c(0, 1),
                       expand = c(0.01, 0),
                       breaks = seq(0, 1, .2)) +
    scale_y_continuous(limits = c(0, 1),
                       expand = c(0.01, 0),
                       breaks = seq(0, 1, .2)) +
    labs(x = bquote(Autocorrelation ~ rho[x]),
         y = bquote(Gaussian ~ rho[z]),
         title = main) +
    theme_gray() +
    theme(legend.position = 'bottom',
          strip.background = element_rect(fill = 'grey5'),
          strip.text = element_text(colour = 'grey95'),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 15, face = 'bold'))

  return(p)

}
