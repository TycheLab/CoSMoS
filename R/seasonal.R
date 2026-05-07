# --- Internal seasonal helpers ------------------------------------------------

#' Stratify a time series by season
#'
#' Splits a two-column (date, value) time series into a list of seasonal
#' subsets, plus a parallel list of nonzero-only subsets.
#'
#' @param TS data frame or data table with columns \code{date} and \code{value}
#' @param season character; name of a date-component function
#'   (e.g. \code{"month"}, \code{"week"})
#'
#' @return a list of two elements: \code{[[1]]} named list of full seasonal
#'   subsets; \code{[[2]]} named list of nonzero seasonal subsets
#'
#' @seealso \code{\link{analyzeTS}}
#' @keywords internal
stratifySeasonData <- function(TS, season) {

  TS    <- as.data.frame(TS)
  strat <- split(TS, do.call(season, list(x = TS[, "date"])))
  names(strat) <- paste("data", seq_along(strat), sep = "_")

  nz <- lapply(strat, function(x) x[x[, "value"] > 0, ])
  names(nz) <- paste("data_nz", seq_along(nz), sep = "_")

  structure(.Data = list(strat, nz))
}


#' Calculate seasonal autocorrelation function
#'
#' Computes the empirical ACF for each season of a time series.
#'
#' @param TS data frame or data table with columns \code{date} and \code{value}
#' @param season character; name of a date-component function
#'   (e.g. \code{"month"})
#' @param lag.max integer; maximum lag for the ACF (default 50)
#'
#' @return a named list (one element per season) of numeric ACF vectors
#'   starting at lag 0
#'
#' @seealso \code{\link{analyzeTS}}, \code{\link{fitACS}}
#' @keywords internal
seasonalACF <- function(TS, season, lag.max = 50) {

  TS <- as.data.table(TS)
  TS[, season := do.call(season, list(x = date))]

  lag0 <- NULL
  TS[, lag0 := .I]
  for (i in 1:lag.max) TS[, paste0("lag", i) := lag0 - i]

  out <- lapply(unique(TS[, season]), function(j) {

    index <- as.data.frame(
      TS[season == j, .SD, .SDcols = grep("lag", names(TS))]
    )

    as <- sapply(1:lag.max, function(i) {
      xi <- index[index[paste0("lag", i)] > 0, "lag0"]
      yi <- index[index[paste0("lag", i)] > 0, paste0("lag", i)]
      x  <- unlist(TS[xi, "value"])
      y  <- unlist(TS[yi, "value"])
      cor(x, y, use = "p")
    })

    c(1, as)
  })

  names(out) <- paste(season, unique(TS[, season]))
  out
}


#' Seasonal AR model
#'
#' Generates a seasonal Gaussian AR process over a vector of dates, one season
#' at a time, using precomputed ACS-derived AR coefficients.
#'
#' @param x POSIXct vector of dates for which to generate the process
#' @param ACS named list of ACS vectors (one per season, starting at lag 0)
#' @param season character; name of a date-component function (default
#'   \code{"month"})
#'
#' @return a \code{data.table} with columns \code{date}, \code{gauss}
#'   (Gaussian process values), and \code{season} (season index)
#'
#' @seealso \code{\link{simulateTS}}, \code{\link{analyzeTS}}
#' @keywords internal
#' @import data.table
seasonalAR <- function(x, ACS, season = "month") {

  time <- data.table(time = x)
  y <- s <- n <- . <- id <- value <- NULL

  time[, y := year(time)]
  time[, s := do.call(season, .(time))]
  time[, n := .N, by = .(y, s)]

  d     <- as.data.frame(unique(time[, -1]))
  alpha <- lapply(ACS, YW)

  esd <- vapply(seq_along(alpha), function(i) {
    sqrt(1 - sum(alpha[[i]] * ACS[[i]][2:length(ACS[[i]])]))
  }, numeric(1))

  ## initialise with a short AR(1) burn-in
  init_season <- which(gsub(paste(season, ""), "", names(ACS)) == d[1, "s"])
  init_val    <- AR1(max(sapply(alpha, length)),
                     ACS[[init_season]][2])

  ## pre-allocate output list then rbind once
  chunks <- vector("list", nrow(d))

  for (j in seq_len(nrow(d))) {

    s_j  <- d[j, "s"]
    ss   <- which(gsub(paste(season, ""), "", names(alpha)) == s_j)
    p    <- length(alpha[[ss]])
    n_j  <- d[j, "n"]

    if (j == 1L) {
      val <- c(init_val, numeric(n_j))
    } else {
      ## carry over last p values from previous chunk
      prev <- unlist(lapply(chunks[seq_len(j - 1)], function(ch) ch$value))
      val  <- c(tail(prev, p), numeric(n_j))
    }

    gn    <- rnorm(n_j + p, mean = 0, sd = esd[s_j])
    a_rev <- rev(alpha[[ss]])

    for (i in (p + 1):(n_j + p)) {
      val[i] <- sum(val[(i - p):(i - 1)] * a_rev) + gn[i]
    }

    chunks[[j]] <- data.table(value = val[(p + 1):(n_j + p)], id = s_j)
  }

  out <- rbindlist(chunks)

  data.table(date   = x,
             gauss  = out[, value],
             season = out[, id])
}


# --- Public analysis pipeline -------------------------------------------------

#' Analyse, report, and simulate seasonal time series
#'
#' \code{analyzeTS} automatically performs seasonal analysis, fits distributions
#' and correlation structures. \code{reportTS} visualises the fitted
#' distributions and correlation structures, or returns a table of fitted
#' parameters and descriptive statistics. \code{simulateTS} takes the result of
#' \code{analyzeTS} and generates synthetic realisations.
#'
#' @details
#' In practice, we typically want to simulate a natural process from observed
#' data. \code{analyzeTS} fits a marginal distribution and autocorrelation
#' structure for each season; \code{reportTS} lets you inspect the fit;
#' \code{simulateTS} generates synthetic time series with the same seasonal
#' statistical properties.
#'
#' Recommended distributions by variable type:
#' * *precipitation / streamflow*: \code{ggamma}, \code{burrXII}, \code{burrIII}
#' * *relative humidity*: \code{beta}
#' * *temperature*: \code{norm}
#'
#' @param TS data frame or data table with columns \code{date} and \code{value}
#' @param season character; name of a date-component function
#'   (e.g. \code{"month"}, \code{"week"})
#' @param acsID character; ACS identifier passed to \code{\link{fitACS}}
#' @param lag.max integer; maximum lag for the empirical ACF
#' @param aTS an \code{analyzeTS} result object
#' @param method character; report type — \code{"dist"} for distribution fits,
#'   \code{"acs"} for ACS fits, \code{"stat"} for descriptive statistics table
#' @param from POSIXct; start of simulation period (defaults to start of
#'   observed series)
#' @param to POSIXct; end of simulation period (defaults to end of observed
#'   series)
#' @inheritParams N
#' @inheritParams fitDist
#'
#' @return
#' * \code{analyzeTS}: a list with elements \code{data}, \code{dfits},
#'   \code{afits}, and attributes \code{season}, \code{dist}, \code{acsID},
#'   \code{date}
#' * \code{reportTS}: a \code{ggplot} object (\code{"dist"} or \code{"acs"}
#'   method) or a \code{data.frame} (\code{"stat"} method)
#' * \code{simulateTS}: a \code{data.table} with columns \code{date} and
#'   \code{value}
#'
#' @seealso \code{\link{fitDist}}, \code{\link{fitACS}}, \code{\link{generateTS}}
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
#' data("precip")
#' \donttest{
#' ## Fit seasonal ACSs and distributions to the data
#' a <- analyzeTS(precip)
#'
#' reportTS(a, "dist")  ## seasonal distribution fits
#' reportTS(a, "acs")   ## seasonal ACS fits
#' reportTS(a, "stat")  ## descriptive statistics
#'
#' ## Simulate a time series of the same length
#' sim <- simulateTS(a)
#'
#' precip[, id := "observed"]
#' sim[, id := "simulated"]
#' dta <- rbind(precip, sim)
#'
#' ggplot(dta) +
#'   geom_line(aes(x = date, y = value)) +
#'   facet_wrap(~id, ncol = 1) +
#'   theme_classic()
#'
#' ## Simulate a time series of different length
#' sim <- simulateTS(a,
#'                   from = as.POSIXct("1978-12-01 00:00:00"),
#'                   to   = as.POSIXct("2008-12-01 00:00:00"))
#' }
#' \dontshow{
#' precip <- precip[between(date,
#'   as.POSIXct("1990-1-01", format("%Y-%m-%d"), tz = "America/Regina"),
#'   as.POSIXct("1990-1-5",  format("%Y-%m-%d"), tz = "America/Regina"))]
#' a <- analyzeTS(precip, opts = list("algorithm" = "NLOPT_LN_NELDERMEAD",
#'                                    "maxeval" = 10))
#' }
#'
analyzeTS <- function(TS, season = "month", dist = "ggamma", acsID = "weibull",
                      norm = "N1", n.points = 30, lag.max = 30,
                      constrain = FALSE, opts = NULL) {

  if (is.null(opts)) opts <- formals(fitDist)$opts

  ea <- seasonalACF(TS, season = season, lag.max = lag.max)
  a  <- lapply(ea, function(x) fitACS(x, acsID))

  x  <- stratifySeasonData(TS, season)
  f  <- lapply(x[[2]], function(x) {
    fitDist(x$value, dist,
            norm      = norm,
            n.points  = n.points,
            constrain = constrain,
            opts      = opts)
  })

  structure(.Data = list(data = x, dfits = f, afits = a),
            season = season,
            dist   = dist,
            acsID  = acsID,
            date   = TS[, "date"])
}


#' @rdname analyzeTS
#' @export
reportTS <- function(aTS, method = "dist") {

  dist   <- attr(aTS, "dist")
  acsID  <- attr(aTS, "acsID")
  season <- attr(aTS, "season")

  if (method == "stat") {

    nz <- aTS$data[[2]]

    dp  <- as.data.frame(round(t(sapply(aTS$dfits, function(x) do.call(rbind, x))), 3))
    names(dp) <- getDistArg(dist)

    ap  <- as.data.frame(round(t(sapply(aTS$afits, function(x) do.call(rbind, x))), 3))
    names(ap) <- getACSArg(acsID)

    laux <- t(sapply(nz, function(x) lmom(x$value)))

    l <- round(data.frame(l.var  = laux[, 1] / laux[, 2],
                           l.skew = laux[, 3],
                           l.kurt = laux[, 4]), 2)

    s <- t(round(sapply(nz, function(x) {
      c(mean = mean(x$value, na.rm = TRUE),
        sd   = sd(x$value,   na.rm = TRUE),
        min  = min(x$value,  na.rm = TRUE),
        q    = quantile(x$value, na.rm = TRUE, probs = .05),
        q    = quantile(x$value, na.rm = TRUE, probs = .25),
        q    = quantile(x$value, na.rm = TRUE, probs = .5),
        q    = quantile(x$value, na.rm = TRUE, probs = .75),
        q    = quantile(x$value, na.rm = TRUE, probs = .95),
        max  = max(x$value,  na.rm = TRUE),
        skew = sample.moments(x$value, raw = FALSE, central = FALSE,
                              coef = TRUE)$coefficients[2])
    }), 2))

    err <- round(sapply(aTS$dfits, function(x) attr(x, "nfo")$objective), 4)

    p0 <- 1 - round(
      sapply(nz,                  nrow) /
      sapply(aTS$data[[1]], nrow), 2)

    a <- round(as.data.frame(
      t(sapply(aTS$afits, function(x) attr(x, "eACS")[2:5]))), 2)
    names(a) <- paste0("acs.l.", 2:5)

    stat <- cbind(dist = dist, dp,
                  error = err,
                  acsID = acsID, ap,
                  s, l, p0, a)
    rownames(stat) <- paste(season, seq_len(nrow(stat)), sep = "_")

    out <- stat
  }

  if (method == "dist") {

    f <- aTS$dfits

    CDF <- lapply(f, function(x) {
      edf  <- attr(x, "edf")
      dist <- attr(x, "dist")
      mm   <- do.call(paste0("p", dist), c(list(q = range(edf$value)), x))
      p    <- exp(seq(ifelse(!is.finite(log(mm[1])), .00001, log(mm[1])),
                      log(mm[2]), length.out = 10000))
      data.frame(p     = p,
                 value = do.call(paste0("q", dist), c(list(p = p), x)))
    })
    names(CDF) <- paste(season, seq_along(CDF), sep = "_")

    sea <- NULL
    cdf <- rbindlist(CDF, idcol = "sea")
    cdf[, sea := factor(sea, levels = paste(season, seq_along(CDF), sep = "_"))]

    EDF <- lapply(f, function(x) attr(x, "edf"))
    names(EDF) <- paste(season, seq_along(EDF), sep = "_")

    edf <- rbindlist(EDF, idcol = "sea")
    edf[, sea := factor(sea, levels = paste(season, seq_along(EDF), sep = "_"))]

    out <- ggplot() +
      geom_point(data = edf,
                 aes(x = edf$value, y = log(1 - edf$p), shape = factor(2)),
                 colour = "royalblue4", alpha = .5) +
      scale_shape_manual(values = 19, label = "Empirical", name = NULL) +
      geom_line(data = cdf,
                aes(x = cdf$value, y = log(1 - cdf$p), linetype = factor(1)),
                colour = "red4", lwd = .5, alpha = .75) +
      scale_linetype_manual(values = 1, label = "Fitted", name = NULL) +
      scale_y_continuous(
        breaks = seq(-10, 0, length.out = 5),
        labels = format(exp(seq(log(.0001), log(1), length.out = 5)),
                        scientific = TRUE)) +
      labs(x = "Nonzero values", y = "Exceedance probability",
           title = "Probability distribution fit") +
      theme_gray() +
      facet_wrap(~sea, scales = "free", nrow = 4) +
      theme(legend.position  = "bottom",
            strip.background = element_rect(fill = "grey5"),
            strip.text       = element_text(colour = "grey95"))
  }

  if (method == "acs") {

    a <- aTS$afits

    ACS <- lapply(a, function(x) {
      eACS <- attr(x, "eACS")
      lag  <- 0:(length(eACS) - 1)
      id   <- attr(x, "ID")
      ACS  <- do.call(acs, c(list(id = id, t = lag), x))
      data.frame(lag = lag, ACS = ACS, eACS = eACS)
    })
    names(ACS) <- paste(season, seq_along(ACS), sep = "_")

    sea <- NULL
    acs_dt <- rbindlist(ACS, idcol = "sea")
    acs_dt[, sea := factor(sea, levels = paste(season, seq_along(ACS), sep = "_"))]

    out <- ggplot(acs_dt) +
      geom_point(aes(x = as.factor(acs_dt$lag), y = acs_dt$eACS,
                     shape = factor(2)),
                 colour = "royalblue4", alpha = .5) +
      scale_shape_manual(values = 19, label = "Empirical", name = NULL) +
      geom_line(aes(x = acs_dt$lag + 1, y = acs_dt$ACS, linetype = factor(1)),
                colour = "red4", lwd = .5, alpha = .75) +
      scale_linetype_manual(values = 1, label = "Fitted", name = NULL) +
      labs(x = bquote(lag ~ tau), y = "Autocorrelation",
           title = "Autocorrelation structure fit") +
      theme_grey() +
      facet_wrap(~sea, scales = "free", nrow = 4) +
      theme(legend.position  = "bottom",
            strip.background = element_rect(fill = "grey5"),
            strip.text       = element_text(colour = "grey95"))
  }

  out
}


#' @rdname analyzeTS
#' @export
simulateTS <- function(aTS, from = NULL, to = NULL) {

  dist   <- attr(aTS, "dist")
  acsID  <- attr(aTS, "acsID")
  season <- attr(aTS, "season")
  date   <- data.table(attr(aTS, "date"))

  x <- aTS$data
  f <- aTS$dfits
  a <- aTS$afits

  distbounds <- if (dist == "beta") c(0, 1) else c(-Inf, Inf)

  ## compute p0 directly from stratified data — avoids running reportTS
  p0_vec <- 1 - round(
    sapply(x[[2]], nrow) /
    sapply(x[[1]], nrow), 2)

  ACS <- vector("list", length(x[[1]]))
  for (i in seq_along(x[[1]])) {
    p  <- actpnts(margdist   = attr(f[[i]], "dist"),
                  margarg    = f[[i]],
                  p0         = p0_vec[i],
                  distbounds = distbounds)
    fp <- fitactf(p)

    lag  <- 0:(length(attr(a[[i]], "eACS")) - 1)
    id   <- attr(a[[i]], "ID")
    as_i <- do.call(acs, c(list(id = id, t = lag), a[[i]]))
    ACS[[i]] <- actf(as_i, fp$actfcoef[1], fp$actfcoef[2])
  }
  names(ACS) <- names(a)

  p0 <- uval <- gauss <- value <- . <- season_id <- rn <- NULL

  if (is.null(from)) from <- date[1, date]
  if (is.null(to))   to   <- date[.N, date]

  by <- difftime(date[2, date], date[1, date])

  gausian <- seasonalAR(x = seq(from = from, to = to, by = by),
                        ACS = ACS)

  setkey(gausian, season)

  para <- as.data.table(
    t(sapply(f, function(x) as.matrix(do.call(cbind, x))))
  )
  names(para) <- names(f[[1]])
  para[, season := as.numeric(gsub("data_nz_", "", rownames(para)))]
  para[, p0     := p0_vec]
  setkey(para, season)

  aux <- merge(gausian, para, all.x = TRUE)
  aux <- aux[order(date)]
  aux[, uval := (pnorm(q = gauss) - p0) / (1 - p0)]
  aux[uval < 0, uval := 0]

  d <- getDistArg(dist)

  for (i in para[, season]) {
    trans.para <- para[season == i, !c("p0", "season")]
    aux[season == i,
        value := do.call(paste0("q", dist),
                         args = c(list(p = uval), as.list(trans.para)))]
  }

  aux[, .(date, value)]
}
