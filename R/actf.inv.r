#' Inverse of ACTF auto-correlation transformation function
#'
#' Provides inverse transformation for continuous distributions, based on two parameters
#'
#' @param rhoz marginal correlation value of the parent Gaussian process
#' @param b 1st line parameter
#' @param c 2nd line parameter
#'
#' @export
#'
#' @keywords internal
#'
#' @examples
#'
#' actf.inv(.4, 1, 0)
#'
actf.inv <- function (rhoz, b, c) {

    rhox <- ((rhoz * ((1 + b)^(1 - c) - 1) + 1)^(1/(1 - c)) - 1) / b
    return(rhox)
}


#' Inverse of ACTF auto-correlation transformation function
#'
#' Provides inverse transformation for continuous distributions, based on two parameters
#'
#' @inheritParams actf.inv
#'
#' @export
#'
#' @keywords internal
#'
#' @examples
#'
#' actfdiscrete.inv(.4, .2, 1)
#'
actfdiscrete.inv <- function(rhoz, b, c) {

  rhox <- (1 - (1 - rhoz)^(1/c) )^(1/b)

  return(rhox)
}
