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
#' actfInv(.4, 1, 0)
#'
actfInv <- function (rhoz, b, c) {

    rhox <- ((rhoz * ((1 + b)^(1 - c) - 1) + 1)^(1/(1 - c)) - 1) / b
    return(rhox)
}


#' Inverse of ACTF auto-correlation transformation function
#'
#' Provides inverse transformation for continuous distributions, based on two parameters
#'
#' @inheritParams actfInv
#'
#' @export
#'
#' @keywords internal
#'
#' @examples
#'
#' actfdiscreteInv(.4, .2, 1)
#'
actfdiscreteInv <- function(rhoz, b, c) {

  rhox <- (1 - (1 - rhoz)^(1/c) )^(1/b)

  return(rhox)
}
