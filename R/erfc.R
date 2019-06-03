#' Complementary (inverse) error function
#'
#' @param x number vector
#'
#' @rdname erfc
#' @export
#' @keywords internal
#'
#' @examples
#'
#' erfc(1)
#' inv.erfc(1)
#'
erfc <- function(x) {

  out <- 2 * pnorm(-sqrt(2) * x)

  return(out)
}

#' @rdname erfc
#' @export
inv.erfc <- function (x) {

  x[x < 0 | x > 2] <- NA
  out <- -qnorm(x / 2) / sqrt(2)

  return(out)
}

