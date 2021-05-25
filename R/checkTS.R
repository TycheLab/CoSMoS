#' Check generated timeseries
#'
#' Compares generated time series sample statistics with the theoretically expected values.
#'
#' @param TS generated timeseries
#' @inheritParams moments
#'
#' @export
#'
#' @examples
#'
#' library(CoSMoS)
#'
#' ## check your generated timeseries
#' x <- generateTS(margdist = 'burrXII',
#'                 margarg = list(scale = 1,
#'                                shape1 = .75,
#'                                shape2 = .25),
#'                 acsvalue = acs(id = 'weibull',
#'                                t = 0:30,
#'                                scale = 10,
#'                                shape = .75),
#'                 n = 1000, p = 30, p0 = .5, TSn = 5)
#'
#' checkTS(x)
#'
checkTS <- function(TS, distbounds = c(-Inf, Inf)) {

  if (!is.list(TS)) {

    TS <- list(TS)
  }

  att <- attributes(TS[[1]])

  margdist <- att$margdist
  margarg <- att$margarg
  p0 <- att$p0
  acsvalue <- att$acsvalue

  ac <- sapply(TS, function(x) acf(x, plot = F)$acf)

  out <- data.frame(mean = c(popmean(margdist,
                                     margarg,
                                     distbounds = distbounds,
                                     p0 = p0),
                             sapply(TS, mean)),
                    sd = c(popsd(margdist,
                                 margarg,
                                 distbounds = distbounds,
                                 p0 = p0),
                           sapply(TS, sd)),
                    skew = c(popskew(margdist,
                                     margarg,
                                     distbounds = distbounds,
                                     p0 = p0),
                             sapply(TS, function(x) sample.moments(x,
                                                                   raw = F,
                                                                   central = F,
                                                                   coef = T)$coefficients[2])),
                    p0 = c(p0,
                           sapply(TS, function(x) length(which(x == 0))/length(x))),
                    acf_t1 = c(acsvalue[2],
                               ac[2,]),
                    acf_t2 = c(acsvalue[3],
                               ac[3,]),
                    acf_t3 = c(acsvalue[4],
                               ac[4,]))

  row.names(out) <- c('expected', paste0('simulation', seq_along(TS)))

  out <- round(out, 2)

  structure(.Data = as.matrix(out),
            class = c('checkTS', 'matrix'),
            margdist = margdist,
            margarg = margarg,
            p0 = p0)
}

#' Plot method for check results
#'
#' Plot method for check results.
#'
#' @param x check result
#' @param ... other args
#'
#' @export
#' @import ggplot2
#' @method plot checkTS
#'
#' @examples
#'
#' library(CoSMoS)
#'
#' ## check your generated timeseries
#' x <- generateTS(margdist = 'burrXII',
#'                 margarg = list(scale = 1,
#'                                shape1 = .75,
#'                                shape2 = .15),
#'                 acsvalue = acs(id = 'weibull',
#'                                t = 0:30,
#'                                scale = 10,
#'                                shape = .75),
#'                 n = 1000, p = 30, p0 = .25, TSn = 100)
#'
#' chck <- checkTS(x)
#'
#' plot(chck)
#'
plot.checkTS <- function(x, ...) {

  att <- attributes(x)

  margdist <- att$margdist
  margarg <- att$margarg
  p0 <- att$p0

  dta <- melt(as.data.table(x = x,
                            keep.rownames = TRUE),
              id.vars = 'rn')
  dta.e <- dta[which(dta$rn == 'expected'),]
  dta.s <- dta[which(dta$rn != 'expected'),]

  p <- ggplot() +
    geom_boxplot(data = dta.s,
                 aes_string(x = 'variable',
                            y = 'value',
                            group = 'variable')) +
    geom_point(data = dta.e,
               aes_string(x = 'variable',
                          y = 'value',
                          group = 'variable'),
               size = 2,
               colour = 'red1') +
    facet_wrap('variable',
               scales = 'free',
               nrow = 1) +
    labs(x = '',
         y = '',
         title = paste('Marginal =', margdist),
         subtitle = paste(paste(names(margarg),
                                '=',
                                margarg,
                                collapse = '; '),
                          paste('\np0 =',
                                p0),
                          collapse = ' ')) +
    theme_gray() +
    theme(legend.position = 'bottom',
          strip.background = element_rect(fill = 'grey5'),
          strip.text = element_text(colour = 'grey95'),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

  return(p)
}
