#' Calculate seasonal ACF
#'
#' @param TS time series
#' @param season name of the season
#' @param lag.max max lag for acf
#'
#' @keywords internal
#' @export
#' @import data.table
#'
#' @examples
#'
#' data('precip')
#'
#'\dontshow{
#'  precip <- precip[between(date, as.POSIXct('1990-1-01', format('%Y-%m-%d'), tz = 'America/Regina'),
#'  as.POSIXct('1990-1-10', format('%Y-%m-%d'), tz = 'America/Regina'))]
#'}
#' seasonalACF(precip, 'month')
#'
seasonalACF <- function(TS, season, lag.max = 50) {

  TS <- as.data.table(TS) ## transform TS to data.table for convenience

  TS[, season := do.call(season, list(x = date))] ## call for seasonal index

  lag0 <- NULL ## global variable check
  TS[, lag0 := .I] ## index for lag 0
  x <- sapply(1:lag.max, function(i) TS[, paste0('lag', i) := lag0 - i]) ## index for lag up to lag.max

  out <- lapply(unique(TS[, season]), function(j) {

    index <- as.data.frame(TS[season == j, .SD, .SDcols = c(grep('lag', names(TS)))]) ## seasonal index selection

    as <- sapply(1:lag.max, function(i) {

      xi <- index[index[paste0('lag', i)] > 0, 'lag0'] ## lag 0 index selection
      yi <- index[index[paste0('lag', i)] > 0, paste0('lag', i)] ## desired lag index selection

      x <- unlist(TS[xi, 'value']) ## selection of values by correct index
      y <- unlist(TS[yi, 'value'])

      c <- cor(x, y, use = 'p') ## correlation calculation

    })

    c(1, as) ## addition of lag 0 value ;)
  })

  names(out) <- paste(season, unique(TS[, season])) ## name the results

  return(out) ## and send it out
}
