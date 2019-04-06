#' ACTF auto-correlation transformation function
#'
#' @param rhox marginal correlation value
#' @param b 1st line parameter
#' @param c 2nd line parameter
#'
#' @export
#'
#' @keywords internal
#'
actf <- function(rhox, b, c) {

  rhoz <- ((1 + b*rhox)^(1 - c) - 1)/((1 + b)^(1 - c) - 1)

  return(rhoz)
}

#' ACTF auto-correlation transformation function for discrete distributions
#'
#' @inheritParams actf
#'
#' @export
#'
#' @keywords internal
#'
actfdiscrete <- function(rhox, b, c) {

  rhoz <- 1 - (1 - rhox^b)^c

  return(rhoz)
}
