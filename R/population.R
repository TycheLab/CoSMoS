#' Population statistics
#'
#' Computes theoretical descriptive statistics of a distribution.
#'
#' @inheritParams moments
#' @param stat character; statistic to compute — one of \code{"mean"},
#'   \code{"sd"}, \code{"var"}, \code{"cvar"}, \code{"skew"}, \code{"kurt"}
#'
#' @return scalar numeric
#'
#' @seealso \code{\link{moments}}, \code{\link{checkTS}}
#'
#' @name PopulationStat
#'
#' @export
#'
#' @examples
#'
#' library(CoSMoS)
#'
#' populationstat("mean", "norm", list(mean = 2, sd = 1))
#' populationstat("sd",   "norm", list(mean = 2, sd = 1))
#' populationstat("var",  "norm", list(mean = 2, sd = 1))
#' populationstat("cvar", "norm", list(mean = 2, sd = 1))
#' populationstat("skew", "norm", list(mean = 2, sd = 1))
#' populationstat("kurt", "norm", list(mean = 2, sd = 1))
#'
populationstat <- function(stat = "mean", dist, distarg,
                           p0 = 0, distbounds = c(-Inf, Inf)) {
  do.call(paste0("pop", stat),
          args = list(dist = dist, distarg = distarg,
                      distbounds = distbounds, p0 = p0))
}

#' @rdname PopulationStat
#' @keywords internal
popmean <- function(dist, distarg, p0 = 0, distbounds = c(-Inf, Inf)) {
  x <- moments(dist = dist, distarg = distarg, distbounds = distbounds,
               p0 = p0, central = TRUE, raw = FALSE, coef = FALSE)
  unname(x[["mu"]][1])
}

#' @rdname PopulationStat
#' @keywords internal
popsd <- function(dist, distarg, p0 = 0, distbounds = c(-Inf, Inf)) {
  x <- moments(dist = dist, distarg = distarg, distbounds = distbounds,
               p0 = p0, central = TRUE, raw = FALSE, coef = FALSE)
  unname(sqrt(x[["mu"]][2]))
}

#' @rdname PopulationStat
#' @keywords internal
popvar <- function(dist, distarg, p0 = 0, distbounds = c(-Inf, Inf)) {
  x <- moments(dist = dist, distarg = distarg, distbounds = distbounds,
               p0 = p0, central = TRUE, raw = FALSE, coef = FALSE)
  unname(x[["mu"]][2])
}

#' @rdname PopulationStat
#' @keywords internal
popcvar <- function(dist, distarg, p0 = 0, distbounds = c(-Inf, Inf)) {
  x <- moments(dist = dist, distarg = distarg, distbounds = distbounds,
               p0 = p0, coef = TRUE, raw = FALSE, central = FALSE)
  unname(x[["coefficients"]][1])
}

#' @rdname PopulationStat
#' @keywords internal
popskew <- function(dist, distarg, p0 = 0, distbounds = c(-Inf, Inf)) {
  x <- moments(dist = dist, distarg = distarg, distbounds = distbounds,
               p0 = p0, coef = TRUE, raw = FALSE, central = FALSE)
  unname(x[["coefficients"]][2])
}

#' @rdname PopulationStat
#' @keywords internal
popkurt <- function(dist, distarg, p0 = 0, distbounds = c(-Inf, Inf)) {
  x <- moments(dist = dist, distarg = distarg, distbounds = distbounds,
               p0 = p0, coef = TRUE, raw = FALSE, central = FALSE)
  unname(x[["coefficients"]][3])
}
