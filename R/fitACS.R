#' Autocorrelation structure fit
#'
#' @param acf vector of autocorrelation function values from lag 0
#' @param ID ACS id
#' @param start starting parameter value
#' @param lag acf lag
#'
#' @export
#' @import stats
#'
#' @examples
#'
#' x <- AR1(1000, .8)
#'
#' acsfit <- fitACS(acf(x, plot = FALSE)$acf, 'weibull', c(1, 1))
#'
fitACS <- function(acf, ID, start = NULL, lag = NULL) {

  if(is.null(lag)) { ## check if lag is unspecified

    if(length(which(acf <= .01)) == 0) { ## if it is and if it never drops below .01 use the entire ACS

      eACS <- acf
    } else { ## if it drops below .01 use only higher values

      eACS <- acf[1:which(acf <= .01)[1]]
    }
  } else {

    eACS <- acf[1:(lag + 1)]
  }


  if(is.null(start)) { ## if there are no starting parameter values start all from 1

    start <- rep(1, length(getACSArg(ID)))
  }

  par <- optim(par = start, ## startin g parameter value
               fn = optimACS, ## function to be minimized
               id = ID, ## ACS ID
               eACS = eACS, ## empirical ACS
               error = 'MSE')$par ## error crit

  out <- as.list(setNames(par, ## list of fitted pars
                          getACSArg(ID)))

  structure(.Data = out, ## outgoing structure
            ID = ID,
            eACS = eACS,
            class = 'fitACS')
}

#' Auxiliary function passed to fitACS
#'
#' @param par parameter value
#' @param id ACS id
#' @param eACS empirical ACS
#' @param error which error to minimize
#'
#' @rdname fitACS
#' @keywords internal
#' @export
#'
optimACS <- function(par, id, eACS, error = 'MSE') {

  arg <- as.list(setNames(par, ## get list of acs arguments
                          getACSArg(id)))
  ACS <- do.call(acs, ## calculate theoretical acs
                 c(list(id = id,
                        t = 0:(length(eACS) - 1)), arg))

  do.call(error, list(x = ACS, y = eACS)) ## calculate zhe error
}

#' Plot method for fitACS
#'
#' @param x fitACS obbject
#' @param ... other args
#'
#' @keywords internal
#' @export
#' @import ggplot2
#' @method plot fitACS
#'
#' @examples
#'
#' x <- AR1(1000, .8)
#'
#' acsfit <- fitACS(acf(x, plot = FALSE)$acf, 'weibull', c(1, 1))
#'
#' plot(acsfit)
#'
plot.fitACS <- function(x, ...) {

  eACS <- attr(x, 'eACS') ## get empirical ACS
  lag <- 0:(length(eACS) - 1) ## get lag
  id <- attr(x, 'ID') ## get ASC ID

  ACS <- do.call(acs, ## get the theoretical ACS
                 c(list(id = id,
                        t = lag), x))

  df <- data.frame(lag = lag, ## put the empirical and theoretical values into one dataframe
                   ACS = ACS,
                   eACS = eACS)

  p <- ggplot(df) + ## plot it :)
    geom_line(aes(x = lag, ## theoretical
                  y = ACS),
              colour = 'grey25',
              lwd = .75,
              alpha = .75) +
    geom_point(aes(x = lag, ## empirical
                   y = eACS),
               colour = 'red4',
               alpha = .5) +
    labs(x = bquote(lag ~ tau),
         y = 'Autocorrelation',
         title = '') +
    theme_grey()

  return(p)
}
