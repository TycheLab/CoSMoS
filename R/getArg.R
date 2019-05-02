#' Get names of distribution function arguments
#'
#' @param dist distribution name
#'
#' @keywords internal
#' @export
#'
#' @examples
#'
#' getDistArg('norm')
#'
getDistArg <- function(dist) {

  if (exists(paste0('d', dist))) { ## see whether CDF is defined

    x <- unlist(strsplit(deparse(args(paste0('q', dist)))[1], ',')) ## get string of all arguments
    out <- gsub(' |= [0-9]+', ## strip the arguments of the balast
                '',
                x[grep('function|lower.tail|log.p = FALSE|/',
                       x,
                       invert = TRUE)])

    return(out)
  } else {

    message('Distribution in not defined')
  }
}

#' Get names of autocorrelation structure (ACS) function arguments
#'
#' @param id ACS id
#'
#' @keywords internal
#' @export
#'
#' @examples
#'
#' getACSArg('weibull')
#'
getACSArg <- function(id) {

  if (exists(paste0('acf', id))) { ## is the ACS function defined?

    x <- unlist(strsplit(deparse(args(paste0('acf', id)))[1], ',')) ## get string of all arguments
    out <- gsub(' |= [0-9]+|)', ## get rid of the balast
                '',
                x[grep('function',
                       x,
                       invert = TRUE)])

    return(out)
  } else {

    message('ACS in not defined')
  }
}
