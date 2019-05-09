#' The Functions analyzeTS, reportTS, and simulateTS
#'
#' Provide a complete set of tools to make time series analysis a piece of cake - analyzeTS() automatically performs seasonal analysis, fits distributions and correlation structures, reportTS provides visualizations of the fitted distributions and correlation structures, and a table with the values of the fitted parameters and basic descriptive statistics, simulateTS automatically takes the results of the analyseTS and generates syntetic ones.
#'
#' In practice, we usually want to simulate a natural process using some sampled time series. To generate a synthetic time series with similar characteristics to the observed values, we have to determine marginal distribution, autocorrelation structure and probability zero for each individual month. This can is done by fitting distributions and autocorrelation structures with analyzeTS(). Result can be checked with reportTS(). Syynthetic time series with the same statistical properties can be produced with simulateTS().
#'
#' Recomended distributions for variables:
#'  * _precipitation_: ggamma (Generalized Gamma), burr### (Burr type)
#'  * _streamflow_: ggamma (Generalized Gamma), burr### (Burr type)
#'  * _relative humidity_: beta
#'  * _temperature_: norm (Normal distribution)
#'
#' @param TS time series in format - date, value
#' @param season name of the season (e.g. month, week)
#' @param acsID ID of the autocorrelation structure to be fitted
#' @param lag.max max lag for the empirical autocorrelation structure
#' @param aTS analyzed timeseries
#' @param method report method - 'dist' for distribution fits, 'acs' for ACS fits and 'stat' for basic statistical report
#' @param from starting date/time of the simulation
#' @param to end date/time of the simulation
#' @inheritParams N
#' @inheritParams fitDist
#'
#' @rdname analyzeTS
#' @export
#' @import ggplot2 data.table
#'
#' @examples
#'
#' library(CoSMoS)
#'
#' ## Load data included in the package
#' ## (to find out more about the data use ?precip)
#' data('precip')
#' \donttest{
#' ## Fit seasonal ACSs and distributions to the data
#' a <- analyzeTS(precip)
#'
#' reportTS(a, 'dist') ## show seasonal distribution fit
#' reportTS(a, 'acs') ## show seasonal ACS fit
#' reportTS(a, 'stat') ## display basic descriptive statisctics
#'
#' ######################################
#' ## 'duplicate' analyzed time series ##
#' sim <- simulateTS(a)
#'
#' ## plot the result
#' precip[, id := 'observed']
#' sim[, id := 'simulated']
#'
#' dta <- rbind(precip, sim)
#'
#' ggplot(dta) +
#'   geom_line(aes(x = date, y = value)) +
#'   facet_wrap(~id, ncol = 1) +
#'   theme_classic()
#'
#' ################################################
#' ## or simulate timeseries of different length ##
#' sim <- simulateTS(a,
#'                   from = as.POSIXct('1978-12-01 00:00:00'),
#'                   to = as.POSIXct('2008-12-01 00:00:00'))
#'
#' ## and plot the result
#' precip[, id := 'observed']
#' sim[, id := 'simulated']
#'
#' dta <- rbind(precip, sim)
#'
#' ggplot(dta) +
#'   geom_line(aes(x = date, y = value)) +
#'   facet_wrap(~id, ncol = 1) +
#'   theme_classic()
#'}
#' \dontshow{
#' ## test for one month to make it fast
#' precip <- precip[between(date, as.POSIXct('1990-1-01', format('%Y-%m-%d'), tz = 'America/Regina'), as.POSIXct('1990-1-5', format('%Y-%m-%d'), tz = 'America/Regina'))]
#' a <- analyzeTS(precip)
#'}
#'
analyzeTS <- function(TS, season = 'month', dist = 'ggamma', acsID = 'weibull', norm = 'N2', n.points = 30, lag.max = 30, constrain = FALSE){

  ea <- seasonalACF(TS, ## calculate seasonal empirical ACS
                    season = season,
                    lag.max = lag.max)

  a <- lapply(ea, function(x) fitACS(x, ## fit empirical ACS by ID
                                     acsID))

  x <- stratifySeasonData(TS, season) ## split data to seasons

  f <- lapply(x[[2]], function(x) fitDist(x$value, ## fit the seasonal data
                                          dist,
                                          norm = norm,
                                          n.points = n.points,
                                          constrain = constrain))

  structure(.Data = list(data = x, ## send the result out
                         dfits = f,
                         afits = a),
            season = season,
            dist = dist,
            acsID = acsID,
            date = TS[, 'date'])
}

#' @rdname analyzeTS
#' @export
reportTS <- function(aTS, method = 'dist') {

  dist <- attr(aTS, 'dist')
  acsID <- attr(aTS, 'acsID')
  season <- attr(aTS, 'season')

  if((method == 'stat')) { # | (method == 'all')) {

    nz <- aTS$data[[2]]

    dp <- as.data.frame(round(t(sapply(aTS$dfits, function(x) do.call(rbind, x))), 3))
    names(dp) <- getDistArg(dist)

    ap <- as.data.frame(round(t(sapply(aTS$afits, function(x) do.call(rbind, x))), 3))
    names(ap) <- getACSArg(acsID)

    laux <- t(sapply(nz, function(x) {

      lmom(x$value)
    }))

    l <- round(data.frame(l.var = laux[, 1]/laux[, 2],
                          l.skew = laux[, 3],
                          l.kurt = laux[, 4]), 2)

    s <- t(round(sapply(nz, function(x) {

      c(mean = mean(x$value, na.rm = T),
        sd = sd(x$value, na.rm = T),
        min = min(x$value, na.rm = T),
        q = quantile(x$value, na.rm = T, probs = .05),
        q = quantile(x$value, na.rm = T, probs = .25),
        q = quantile(x$value, na.rm = T, probs = .5),
        q = quantile(x$value, na.rm = T, probs = .75),
        q = quantile(x$value, na.rm = T, probs = .95),
        max = max(x$value, na.rm = T),
        skew = sample.moments(x$value,
                              raw = F,
                              central = F,
                              coef = T)$coefficients[2])
    }), 2))

    err <- round(sapply(aTS$dfits, function(x) attr(x, 'err')), 4)

    p0 <- 1 - round(sapply(nz, dim)[1,]/sapply(aTS$data[[1]], dim)[1,], 2)

    a <- round(as.data.frame(t(sapply(aTS$afits, function(x) {attr(x, 'eACS')[2:5]}))), 2)
    names(a) <- paste0('acs.l.', 2:5)

    stat <- cbind(dist = dist, dp,
                  error = err,
                  acsID = acsID, ap,
                  s,
                  l,
                  p0,
                  a)

    rownames(stat) <- paste(season, 1:dim(stat)[1], sep = '_')

    out <- stat
  }

  if((method == 'dist')) { # | (method == 'all')) {

    f <- aTS$dfits

    CDF <- lapply(f, function(x) {

      edf <- attr(x, 'edf')
      dist <- attr(x, 'dist')

      mm <- do.call(paste0('p', dist),
                    c(list(q = range(edf$value)),
                      x))
      p <- exp(seq(ifelse(!is.finite(log(mm[1])),
                          .00001,
                          log(mm[1])),
                   log(mm[2]),
                   length.out = 10000))
      cdf <- data.frame(p = p,
                        value = do.call(paste0('q', dist),
                                        c(list(p = p), x)))
      cdf
    })

    names(CDF) <- paste(season,
                        seq_along(CDF),
                        sep = '_')

    sea <- NULL
    cdf <- rbindlist(CDF,
                     idcol = 'sea')
    cdf[, sea := factor(sea,
                        levels = paste(season,
                                       seq_along(CDF),
                                       sep = '_'))]

    EDF <- lapply(f, function(x) {

      attr(x, 'edf')
    })

    names(EDF) <- paste(season,
                        seq_along(EDF),
                        sep = '_')

    edf <- rbindlist(EDF,
                     idcol = 'sea')
    edf[, sea := factor(sea,
                        levels = paste(season,
                                       seq_along(EDF),
                                       sep = '_'))]

    df <- ggplot() +
      geom_point(data = edf,
                 aes(x = edf$value,
                     y = log(1 - edf$p),
                     shape = factor(2)),
                 colour = 'royalblue4',
                 alpha = .5) +
      scale_shape_manual(values = 19,
                         label = 'Empirical',
                         name = NULL) +
      geom_line(data = cdf,
                aes(x = cdf$value,
                    y = log(1 - cdf$p),
                    linetype = factor(1)),
                colour = 'red4',
                lwd = .5,
                alpha = .75) +
      scale_linetype_manual(values = 1,
                            label = 'Fitted',
                            name = NULL) +
      scale_y_continuous(breaks = seq(-10,
                                      0,
                                      length.out = 5),
                         labels = format(exp(seq(log(.0001),
                                                 log(1),
                                                 length.out = 5)),
                                         scientific = TRUE)) +
      labs(x = 'Nonzero values',
           y = 'Exceedence probability',
           title = 'Probability distribution fit') +
      theme_gray() +
      facet_wrap(~sea,
                 scales = 'free',
                 nrow = 4) +
      theme(legend.position = 'bottom',
            strip.background = element_rect(fill = 'grey5'),
            strip.text = element_text(colour = 'grey95'))

    out <- df
  }

  if((method == 'acs')) { # | (method == 'all')) {

    a <- aTS$afits

    ACS <- lapply(a, function(x) {

      eACS <- attr(x, 'eACS')
      lag <- 0:(length(eACS) - 1)
      id <- attr(x, 'ID')

      ACS <-   ACS <- do.call(acs,
                              c(list(id = id,
                                     t = lag), x))
      ac <- data.frame(lag = lag,
                       ACS = ACS,
                       eACS = eACS)
      ac
    })

    names(ACS) <- paste(season,
                        seq_along(ACS),
                        sep = '_')

    sea <- NULL
    acs <- rbindlist(ACS,
                     idcol = 'sea')
    acs[, sea := factor(sea,
                        levels = paste(season,
                                       seq_along(ACS),
                                       sep = '_'))]

    ac <- ggplot(acs) +
      geom_point(aes(x = as.factor(acs$lag),
                     y = acs$eACS,
                     shape = factor(2)),
                 colour = 'royalblue4',
                 alpha = .5) +
      scale_shape_manual(values = 19,
                         label = 'Empirical',
                         name = NULL) +
      geom_line(aes(x = acs$lag + 1,

                    y = acs$ACS,
                    linetype = factor(1)),
                colour = 'red4',
                lwd = .5,
                alpha = .75) +
      scale_linetype_manual(values = 1,
                            label = 'Fitted',
                            name = NULL) +
      labs(x = bquote(lag ~ tau),
           y = 'Autocorrelation',
           title = 'Autocorrelation structure fit') +
      theme_grey() +
      facet_wrap(~sea,
                 scales = 'free',
                 nrow = 4) +
      theme(legend.position = 'bottom',
            strip.background = element_rect(fill = 'grey5'),
            strip.text = element_text(colour = 'grey95'))

    out <- ac
  }

  return(out)
}

#' @rdname analyzeTS
#' @export
simulateTS <- function(aTS, from = NULL, to = NULL) {

  dist <- attr(aTS, 'dist') ## get necesary info from attributes
  acsID <- attr(aTS, 'acsID')
  season <- attr(aTS, 'season')
  date <- attr(aTS, 'date')

  x <- aTS$data

  r <- reportTS(aTS, method = 'stat')

  f <- aTS$dfits
  a <- aTS$afits

  ACS <- list()

  #################
  if (dist == 'beta') {

    distbounds = c(0, 1)
  } else {

    distbounds = c(-Inf, Inf)
  }
  ################

  for (i in seq_along(x[[1]])) {

    p <- actpnts(margdist = attr(f[[i]], 'dist'), ## caculate acti points for each season
                 margarg = f[[i]],
                 p0 = r[i, 'p0'],
                 distbounds = distbounds)
    fp <- fitactf(p) ## fit acti points
    # plot(fp)

    lag <- 0:(length(attr(a[[i]], 'eACS')) - 1) ## get correct lag
    id <- attr(a[[i]], 'ID')

    as <- do.call(acs, c(list(id = id, t = lag), a[[i]])) ## get ACS
    ACS[[i]] <- actf(as, fp$actfcoef[1], fp$actfcoef[2]) ## transform ACS
  }

  names(ACS) <- names(a)

  p0 <- uval <- gauss <- value <- . <- NULL ## global variable check

  if (is.null(from)) {

    from <- date[1, date]
  }

  if (is.null(to)) {

    to <- date[.N, date]
  }

  by <- (strsplit(format(difftime(date[2, date], ## time seq step
                                  date[1, date])), ' ')[[1]][2])

  gausian <- seasonalAR(x = seq(from = from, ## generate seasonal gaussian process
                                to = to,
                                by = by),
                        ACS = ACS)

  setkey(gausian, season)

  para <- as.data.table(t(sapply(f, function(x) as.matrix(do.call(cbind, x))))) ## get distribution pars
  names(para) <- names(f[[1]])
  para[, season := as.numeric(gsub('data_nz_', '', rownames(para)))]
  para[, p0 := r[, 'p0']]
  setkey(para, season)

  aux <- merge(gausian, para, all.x = T) ## merge gaussian process with parameters
  aux <- aux[order(date)]
  aux[, uval := (pnorm(gauss) - p0)/(1 - p0)] ## calculate intermitent process
  aux[uval < 0, uval := 0]

  d <- getDistArg(dist)

  l <- lapply(seq_along(d), function(i) as.data.frame(aux)[, d[i]]) ## get dist para
  names(l) <- d

  suppressWarnings(x <- do.call(paste0('q', dist), args = c(list(p = aux[, uval]), l))) ## transform gaussian process
  aux[, value := x]

  out <- aux[, .(date, value)] ## select date amd value

  return(out) ## send it out
}
