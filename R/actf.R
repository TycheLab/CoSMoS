#' ACTF auto-correlation transformation function
#'
#' Provides tranformation for continuous distributions, based on two parameters
#'
#' @param rhox marginal correlation value
#' @param b 1st line parameter
#' @param c 2nd line parameter
#'
#' @export
#'
#' @keywords internal
#'
#' @examples
#'
#' actf(.4, 1, 0)
#'
actf <- function(rhox, b, c) {

  rhoz <- ((1 + b * rhox) ^ (1 - c) - 1)/((1 + b) ^ (1 - c) - 1)

  return(rhoz)
}

#' ACTF auto-correlation transformation function for discrete distributions
#'
#' Provides tranformation for discrete distributions, based on two parameters
#'
#' @inheritParams actf
#'
#' @export
#'
#' @keywords internal
#'
#' @examples
#'
#' actfdiscrete(.4, .2, 1)
#'
actfdiscrete <- function(rhox, b, c) {

  rhoz <- 1 - (1 - rhox ^ b) ^ c

  return(rhoz)
}
