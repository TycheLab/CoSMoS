#' Population statistics
#'
#' Provides theoretical descriptive statistics
#'
#' @inheritParams moments
#' @param stat define what you what to calculate - possible population desc. statistics ('mean', 'sd', 'var', 'cvar', 'skew', 'kurt')
#'
#' @name PopulationStat
#'
#' @export
#'
#' @examples
#'
#' library(CoSMoS)
#'
#' ## check population statistics
#' populationstat('mean', 'norm', list(mean = 2, sd = 1))
#' populationstat('sd', 'norm', list(mean = 2, sd = 1))
#' populationstat('var', 'norm', list(mean = 2, sd = 1))
#' populationstat('cvar', 'norm', list(mean = 2, sd = 1))
#' populationstat('skew', 'norm', list(mean = 2, sd = 1))
#' populationstat('kurt', 'norm', list(mean = 2, sd = 1))
#'
populationstat <- function(stat = 'mean', dist, distarg, p0 = 0, distbounds = c(-Inf, Inf)) {

  s <- do.call(paste0('pop', stat),
               args = list(dist = dist,
                           distarg = distarg,
                           distbounds = distbounds,
                           p0 = p0))

  return(s)
}

#' @rdname PopulationStat
#' @export
#' @keywords internal

popmean <- function(dist, distarg, p0 = 0, distbounds = c(-Inf, Inf)) {

  x <- moments(dist = dist,
               distarg = distarg,
               distbounds = distbounds,
               p0 = p0,
               central = T,
               raw = F, coef = F)

  return(unname(x[['mu']][1]))
}

#' @rdname PopulationStat
#' @export
#' @keywords internal

popsd <- function(dist, distarg, p0 = 0, distbounds = c(-Inf, Inf)) {

  x <- moments(dist = dist,
               distarg = distarg,
               distbounds = distbounds,
               p0 = p0,
               central = T,
               raw = F, coef = F)

  return(unname(sqrt(x[['mu']][2])))
}

#' @rdname PopulationStat
#' @export
#' @keywords internal

popvar <- function(dist, distarg, p0 = 0, distbounds = c(-Inf, Inf)) {

  x <- moments(dist = dist,
               distarg = distarg,
               distbounds = distbounds,
               p0 = p0,
               central = T,
               raw = F, coef = F)

  return(unname(x[['mu']][2]))
}

#' @rdname PopulationStat
#' @export
#' @keywords internal

popcvar <- function(dist, distarg,p0 = 0, distbounds = c(-Inf, Inf)) {

  x <- moments(dist = dist,
               distarg = distarg,
               distbounds = distbounds,
               p0 = p0,
               coef = T,
               raw = F, central = F)

  return(unname(x[['coefficients']][1]))
}

#' @rdname PopulationStat
#' @export
#' @keywords internal

popskew <- function(dist, distarg, p0 = 0, distbounds = c(-Inf, Inf)) {

  x <- moments(dist = dist,
               distarg = distarg,
               distbounds = distbounds,
               p0 = p0,
               coef = T,
               raw = F, central = F)

  return(unname(x[['coefficients']][2]))
}

#' @rdname PopulationStat
#' @export
#' @keywords internal

popkurt <- function(dist, distarg, p0 = 0, distbounds = c(-Inf, Inf)) {

  x <- moments(dist = dist,
               distarg = distarg,
               distbounds = distbounds,
               p0 = p0,
               coef = T,
               raw = F, central = F)

  return(unname(x[['coefficients']][3]))
}
