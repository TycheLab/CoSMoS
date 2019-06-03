#' Norm - function to be minimized during distrubution fit
#'
#' @param par parameter value
#' @param val empirical value
#' @param dist name of the distribution to be fitted
#' @param norm norm used for distribution fitting - id ('N1', 'N2', 'N3', 'N4')
#' @param n.points number of points to be subsetted from ecdf
#'
#' @keywords internal
#' @export
#'
#' @examples
#'
#' N(c(0,1), rnorm(1000), 'norm', 'N2', 30)
#'
N <- function(par, val, dist, norm, n.points) {

  edf <- ECDF(val) ## get ecdf

  aux <- dim(edf)[1] ## get number of values

  if(aux != n.points) { ## select number of values for norms

    if(aux < n.points) {

      n.points <- aux
    }

    edf <- edf[seq(from = 1,
                   to = aux,
                   length.out = n.points),]
  }

  a <- getDistArg(dist) ## get distribution arguments

  arg <- as.list(setNames(par, ## ready dist arg for do.call
                          a))
  if((norm == 'N1') | (norm == 'N2')) {

    FN <- edf$value ## empirical values
    F <- do.call(paste0('q', ## theoretical values
                        dist),
                 args = c(list(edf$p),
                          arg))


    if(norm == 'N1') { ## N1

      out <- rMSE(F, FN)
    }

    if(norm == 'N2') { ## N2

      out <- MSE(F, FN)
    }
  }

  if((norm == 'N3') | (norm == 'N4')) {

    Xi <- edf$p ## empirical probs
    Xu <- do.call(paste0('p', ## theoretical probs
                         dist),
                  args = c(list(edf$value),
                           arg))

    if(norm == 'N3') { ## N3

      out <- rMSE(Xu, Xi)
    }

    if(norm == 'N4') { ## N4

      out <- MSE(Xu, Xi)
    }
  }

  if(!is.finite(out)) {

    out <- 10E9 ## if the norm result is not finite return a big number
  }

  return(out)
}
