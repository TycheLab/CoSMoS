## getArg.R
## Internal helpers to introspect argument names of distribution and ACS
## functions by name.  Both functions look up targets from the CoSMoS
## namespace directly so the target functions do not need to be exported.


#' Get argument names of a distribution function
#'
#' @description
#' Returns the non-standard argument names of the quantile function
#' \code{q<dist>} (i.e. everything except \code{p}, \code{lower.tail},
#' and \code{log.p}).
#'
#' @param dist  distribution name (character), e.g. \code{"paretoII"}.
#'
#' @return A character vector of argument names, or \code{NULL} (with a
#'   message) if the distribution is not defined.
#'
#' @keywords internal

getDistArg <- function(dist) {

  fn_name <- paste0("q", dist)
  ns      <- asNamespace("CoSMoS")

  if (!existsFunction(fn_name) &&
      !exists(fn_name, envir = ns, inherits = FALSE)) {
    message("Distribution is not defined")
    return(NULL)
  }

  ## prefer the package-namespace copy; fall back to search path
  fn <- if (exists(fn_name, envir = ns, inherits = FALSE)) {
    get(fn_name, envir = ns, inherits = FALSE)
  } else {
    get(fn_name, inherits = TRUE)
  }

  x <- unlist(strsplit(deparse(args(fn))[1L], ","))

  gsub(" |= [0-9]+", "",
       x[grep("function|lower\\.tail|log\\.p = FALSE|/",
               x, invert = TRUE)])
}


#' Get argument names of an ACS function
#'
#' @description
#' Returns the non-\code{t} argument names of the ACS function
#' \code{acf<id>} (i.e. the parameter names only).
#'
#' @param id  ACS id (character), e.g. \code{"weibull"}.
#'
#' @return A character vector of parameter names, or \code{NULL} (with a
#'   message) if the ACS is not defined.
#'
#' @keywords internal

getACSArg <- function(id) {

  fn_name <- paste0("acf", id)
  ns      <- asNamespace("CoSMoS")

  if (!exists(fn_name, envir = ns, inherits = FALSE)) {
    message("ACS is not defined")
    return(NULL)
  }

  fn <- get(fn_name, envir = ns, inherits = FALSE)
  x  <- unlist(strsplit(deparse(args(fn))[1L], ","))

  gsub(" |= [0-9]+|\\)", "",
       x[grep("function", x, invert = TRUE)])
}
