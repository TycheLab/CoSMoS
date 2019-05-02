#' Stratiffy timeseries by season
#'
#' @inheritParams seasonalACF
#'
#' @keywords internal
#' @export
#'
#' @examples
#'
#' x <- data.frame(date = seq(Sys.Date(), by = 'day', length.out = 1000),
#'                 value = rnorm(1000))
#'
#' stratifySeasonData(x, 'month')
stratifySeasonData <- function(TS, season) {

  TS <- as.data.frame(TS) ## convert TS to dataframe for convenience

  strat <- split(TS, do.call(season, list(x = TS[, 'date']))) ## split TS to list my season
  names(strat) <- paste('data', seq_along(strat), sep = '_') ## name the list elements

  nz <- lapply(strat, function(x) { ## create another list for nonzero values

    x <- x[x[,'value'] > 0,]; x
  })

  names(nz) <- paste('data_nz', seq_along(nz), sep = '_') ## name the nonzero values

  structure(.Data = list(strat, ## send it out
                         nz))
}
