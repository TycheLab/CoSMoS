#' Seasonal AR model
#'
#' @param x vector of dates for gaussian process generation
#' @param ACS list of ACS for each season
#' @param season season name
#'
#' @keywords internal
#' @import data.table
#' @export
#'
#' @examples
#'
#' data('precip')
#'
#'\dontshow{
#'  precip <- precip[between(date, as.POSIXct('1990-1-01', format('%Y-%m-%d'), tz = 'America/Regina'),
#'  as.POSIXct('1990-1-10', format('%Y-%m-%d'), tz = 'America/Regina'))]
#'}
#' x <- seasonalACF(precip, 'month')
#'
#' seasonalAR(precip$date, x)
#'
seasonalAR <- function(x, ACS, season = 'month') {

  time <- data.table(time = x)

  y <- s <- n <- . <- id <- value <- NULL ## global variable check

  time[, y := year(time)]
  time[, s := do.call(season, .(time))]
  time[, n := .N, by = .(y, s)]

  d <- as.data.frame(unique(time[, -1])) ## main dataframe of seasons and number of values to be generated per season

  alpha <- lapply(ACS, YW) ## pars for the seasonal models

  out <- data.table(value = AR1(max(sapply(alpha, length)),
                                ACS[[which(gsub(paste(season, ''), '', names(ACS)) == d[1, 's'])]][2])) ## overal init values

  out[, id := 0]

  esd <- c()

  for (i in seq_along(alpha)) { ## sd for gaussian noise

    esd[i] <- sqrt(1 - sum(alpha[[i]]*ACS[[i]][2:length(ACS[[i]])]))
  }


  for (j in 1:dim(d)[1]) {

    s <- d[j, 's'] ## season selection

    ss <- which(gsub(paste(season, ''), '', names(alpha)) == s) ## season possition

    p <- length(alpha[[ss]]) ## model order

    val <- out[(dim(out)[1] + 1 - p):dim(out)[1], value] ## initial value

    aux <- length(val)

    n <- unlist(d[j, 'n']) ## number od values to gen in a seasonal run

    gn <- rnorm(n + p, ## gaussian noise
                mean = 0,
                sd = esd[s])

    a.rev <- rev(alpha[[ss]]) ## alpha in the correct order

    for (i in (p + 1):(n + p)) { ## AR
      val[i] <- sum(val[(i - p):(i - 1)]*a.rev) + gn[i]
    }

    out <- rbind(out,
                 data.table(value = val[-1:-aux],
                            id = s)) ## concenate data

  }

  out <- out[id != 0, ]

  dt <- data.table(date = x,
                   gauss = out[, value],
                   season = out[, id])

  return(dt)
}
#' Yule-Walker solver
#'
#' @param ACS vector of ACS values
#'
#' @keywords internal
#' @export
#'
#' @examples
#'
#' YW(rev(exp(seq(-1, 0, .1))))
#'
YW <- function(ACS) {

  p <- length(ACS) - 1

  P <- matrix(NA, p, p) ## cov matrix generation

  for (i in 1:p) {
    P[i, i:p] <- ACS[1:(p - i + 1)]
    P[i, 1:i] <- ACS[i:1]
  }

  rho <- matrix(ACS[2:(p + 1)], p, 1)

  m <- solve(P, rho) ## Yule-Walker

  return(m)
}
