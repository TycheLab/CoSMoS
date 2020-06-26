#' Quick visualization of basic timeseries properties
#'
#' Return timeseries diagram, empirical density function, and empirical autocorrelation function
#'
#' @param TS timeseries to plot
#' @param xlimdf x-axis limit for the plot the empirical density function (default set to NULL). It can be specified for better visualization of highly skewed distributions)
#' @param ci confidence interval around the zero autocorrelation value (default set to 0.95, i.e. 95% CI)
#'
#' @name quickTSPlot
#'
#' @import cowplot
#'
#' @export
#'
#' @examples
#' no <- 1000
#' ggamma_sim <- rggamma(n = no, scale = 1, shape1 = 1, shape2 = .5)
#' quickTSPlot(ggamma_sim)
#'
quickTSPlot <- function(TS, xlimdf = NULL, ci = 0.95){

  value <- Lag <- Acf <- NULL
  if(class(TS)[1] %in% c("data.frame", "data.table") ) n <- nrow(TS)
  else n <- length(TS)
  dta <- data.frame(time = c(1:n), value = TS)
  if(is.null(xlimdf)) xlimdf <- range(TS)
  p1 <- ggplot(dta, aes(x = time, y = value)) + geom_line()
  p2 <- ggplot(dta, aes(x=value)) + geom_density() + xlim(xlimdf)

  acf0 <- acf(TS, plot = FALSE)
  acf1 <- data.frame(Lag = acf0$lag, Acf = acf0$acf)
  clim <- qnorm((1 + ci)/2)/sqrt(acf0$n.used)

  p3 <- ggplot(data = acf1, aes(x = Lag, y = Acf)) +
    geom_hline(aes(yintercept = clim), linetype = 2) +
    geom_hline(aes(yintercept = -clim), linetype = 2) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = Lag, yend = 0))

  ggdraw() +
    draw_plot(p1, x = 0, y = .5, width = 1, height = 0.5) +
    draw_plot(p2, x = 0, y = 0, width = .5, height = .5) +
    draw_plot(p3, x = .5, y = 0, width = .5, height = .5)

}




