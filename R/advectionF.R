#' Advection fields
#'
#' Provides parametric functions that describe different types of advection fields.
#'
#' @param id advection type id (\code{uniform}, \code{rotation}, \code{spiral}, \code{spiralCE}, \code{radial}, and \code{hyperbolic})
#' @param ... other arguments (vector of coordinates and parameters of advection field functions)
#'
#' @name advectionF
#'
#' @import ggplot2 ggquiver
#' @export
#'
#' @references Papalexiou, S.M., Serinaldi, F., Porcu, E. (2021). Advancing
#' Space-Time Simulation of Random Fields: From Storms to Cyclones and Beyond.
#' Water Resources Research, 57, e2020WR029466, \doi{10.1029/2020WR029466}
#'
#' @examples
#'
#' library(ggquiver)
#' library(ggplot2)
#'
#' ## specify coordinates
#' m = 25
#' aux <- seq(0, m - 1, length = m)
#' coord <- expand.grid(aux, aux)
#'
#' ## get the advection field
#' af <- advectionF('spiral',
#'                  spacepoints = coord,
#'                  x0 = floor(m / 2),
#'                  y0 = floor(m / 2),
#'                  a = 3,
#'                  b = 2,
#'                  rotation = 1)
#'
#' ## visualize advection field
#' dta <- data.frame(lon = coord[ ,1], lat = coord[ ,2], u = af[ ,1], v = af[ ,2])
#' ggplot(dta, aes(x = lon, y = lat, u = u, v = v)) +
#' geom_quiver() +
#' theme_light()
#'
advectionF <- function(id, ...) {

  .args <- list(...)
  do.call(paste0('advectionF', id), args = .args)
}


#' Advection fields
#'
#' Provides parametric functions that describe different types of advection fields.
#'
#' @param id advection type id
#' @param arglist list of additional arguments (vector of coordinates and parameters of advection field functions)
#'
#' @name advectionF2
#'
#' @references Papalexiou, S.M., Serinaldi, F., Porcu, E. (2021). Advancing
#' Space-Time Simulation of Random Fields: From Storms to Cyclones and Beyond.
#' Water Resources Research, 57, e2020WR029466, \doi{10.1029/2020WR029466}
#'
#' @keywords internal
#'
advectionF2 <- function(id, arglist) {

  do.call(paste0('advectionF', id), args = arglist)
}


#' Uniform advection field
#'
#' Provides an advection field with constant orthogonal (u and v) components at each grid point. This mimics rigid translation in a given direction according to the components u and v of the velocity vector.
#'
#' @param spacepoints vector of coordinates (2 x d), where d is the number of locations/grid points
#' @param u velocity component along the x axis
#' @param v velocity component along the y axis
#'
#' @name advectionFuniform
#'
#' @import ggplot2
#' @export
#'
#' @references Papalexiou, S.M., Serinaldi, F., Porcu, E. (2021). Advancing
#' Space-Time Simulation of Random Fields: From Storms to Cyclones and Beyond.
#' Water Resources Research, 57, e2020WR029466, \doi{10.1029/2020WR029466}
#'
#' @examples
#'
#' library(ggquiver)
#' library(ggplot2)
#' ## specify coordinates
#' m = 25
#' aux <- seq(0, m - 1, length = m)
#' coord <- expand.grid(aux, aux)
#'
#' af <- advectionFuniform(spacepoints = coord,
#'                        u = 2,
#'                        v = 6)
#'
#' ## visualize advection field
#' dta <- data.frame(lon = coord[ ,1], lat = coord[ ,2], u = af[ ,1], v = af[ ,2])
#' ggplot(dta, aes(x = lon, y = lat, u = u, v = v)) +
#' geom_quiver() +
#' theme_light()
#'
advectionFuniform <- function(spacepoints, u, v) {

  d <- nrow(spacepoints)
  return( data.frame(u = rep(u, d), v = rep(v, d)) )
}


#' Rotational advection field
#'
#' Provides an advection field corresponding to rotation around a specified center.
#'
#' @param spacepoints vector of coordinates (2 x d), where d is the number of locations/grid points
#' @param x0 x coordinate of the center of rotation
#' @param y0 y coordinate of the center of rotation
#' @param a parameter controlling the x component of rotational velocity
#' @param b parameter controlling the y component of rotational velocity
#'
#' @details # Note
#' \itemize{
#' \item if a > 0, b > 0: clockwise rotation around (x0, y0)
#' \item if a < 0, b < 0: counter-clockwise rotation around (x0, y0)
#' }
#'
#' @name advectionFrotation
#'
#' @import ggplot2
#' @export
#'
#' @references Papalexiou, S.M., Serinaldi, F., Porcu, E. (2021). Advancing
#' Space-Time Simulation of Random Fields: From Storms to Cyclones and Beyond.
#' Water Resources Research, 57, e2020WR029466, \doi{10.1029/2020WR029466}
#'
#' @examples
#'
#' library(ggquiver)
#' library(ggplot2)
#' ## specify coordinates
#' m = 25
#' aux <- seq(0, m - 1, length = m)
#' coord <- expand.grid(aux, aux)
#'
#' af <- advectionFrotation(spacepoints = coord,
#'                         x0 = floor(m / 2),
#'                         y0 = floor(m / 2),
#'                         a = 3,
#'                         b = 2)
#'
#' ## visualize advection field
#' dta <- data.frame(lon = coord[ ,1], lat = coord[ ,2], u = af[ ,1], v = af[ ,2])
#' ggplot(dta, aes(x = lon, y = lat, u = u, v = v)) +
#' geom_quiver() +
#' theme_light()
#
advectionFrotation <- function(spacepoints, x0, y0, a, b) {

  x <- spacepoints[,1] - x0
  y <- spacepoints[,2] - y0
  u <-    a * y
  v <- (-b * x)
  return( data.frame(u = u, v = v) )
}


#' Spiraling advection field
#'
#' Provides an advection field corresponding to a spiral motion to/from a specified reference point (sink).
#'
#' @param spacepoints vector of coordinates (2 x d), where d is the number of locations/grid points
#' @param x0 x coordinate of reference point (sink)
#' @param y0 y coordinate of reference point (sink)
#' @param a parameter controlling the x component of rotational velocity
#' @param b parameter controlling the y component of rotational velocity
#' @param rotation parameter controlling the rotational direction. The following combinations hold: \cr
#' \itemize{
#' \item if a > 0, b > 0, and direction = 1: spiraling CLOCKWISE TO (x0, y0)
#' \item if a < 0, b < 0, and direction = 1: spiraling COUNTER-CLOCKWISE FROM (x0, y0)
#' \item if a > 0, b > 0, and direction = 2: spiraling COUNTER-CLOCKWISE TO (x0, y0)
#' \item if a < 0, b < 0, and direction = 2: spiraling CLOCKWISE FROM (x0, y0)
#' }
#'
#' @name advectionFspiral
#'
#' @import ggplot2
#' @export
#'
#' @references Papalexiou, S.M., Serinaldi, F., Porcu, E. (2021). Advancing
#' Space-Time Simulation of Random Fields: From Storms to Cyclones and Beyond.
#' Water Resources Research, 57, e2020WR029466, \doi{10.1029/2020WR029466}
#'
#' @examples
#'
#' library(ggquiver)
#' library(ggplot2)
#' ## specify coordinates
#' m = 25
#' aux <- seq(0, m - 1, length = m)
#' coord <- expand.grid(aux, aux)
#'
#' af <- advectionFspiral(spacepoints = coord,
#'                         x0 = floor(m / 2),
#'                         y0 = floor(m / 2),
#'                         a = 3,
#'                         b = 2,
#'                         rotation = 1)
#'
#' ## visualize advection field
#' dta <- data.frame(lon = coord[ ,1], lat = coord[ ,2], u = af[ ,1], v = af[ ,2])
#' ggplot(dta, aes(x = lon, y = lat, u = u, v = v)) +
#' geom_quiver() +
#' theme_light()
#
advectionFspiral <- function(spacepoints, x0, y0, a, b, rotation = 1) {

  x <- spacepoints[,1] - x0
  y <- spacepoints[,2] - y0
  if (rotation == 1) {
    u <-  ( b * y - a * x)
    v <-  (-a * x - b * y)
  }
  if (rotation == 2) {
    u <-  (-a * x - b * y)
    v <-  ( a * x - b * y)
  }
  return( data.frame(u = u, v = v) )
}



#' Spiraling advection field satisfying continuity equation
#'
#' Provides an advection field corresponding to a spiral motion to/from a specified reference point (sink) satisfying continuity equation (from \href{https://people.sc.fsu.edu/~jburkardt/f77_src/spiral_data/spiral_data.html}{John Burkardt's website}).
#'
#' @param spacepoints vector of coordinates (2 x d), where d is the number of locations/grid points
#' @param a parameter controlling the intensity of rotational velocity (a > 0 clokwise; a < 0 conter-clockwise)
#' @param C parameter ranging in (0, 2*pi)
#'
#' @name advectionFspiralCE
#'
#' @import ggplot2
#' @export
#'
#' @references Papalexiou, S.M., Serinaldi, F., Porcu, E. (2021). Advancing
#' Space-Time Simulation of Random Fields: From Storms to Cyclones and Beyond.
#' Water Resources Research, 57, e2020WR029466, \doi{10.1029/2020WR029466}
#'
#' @examples
#' library(ggquiver)
#' library(ggplot2)
#' ## specify coordinates
#' m = 25
#' aux <- seq(0, m - 1, length = m)
#' coord <- expand.grid(aux, aux)
#'
#' af <- advectionFspiralCE(spacepoints = coord,
#'                         a = 5,
#'                         C = 1)
#'
#' ## visualize advection field
#' dta <- data.frame(lon = coord[ ,1], lat = coord[ ,2], u = af[ ,1], v = af[ ,2])
#' ggplot(dta, aes(x = lon, y = lat, u = u, v = v)) +
#' geom_quiver() +
#' theme_light()
#
advectionFspiralCE <- function(spacepoints, a, C) {

  xaux <- spacepoints[,1]
  yaux <- spacepoints[,2]
  x <- (xaux - min(xaux)) / diff(range(xaux))
  y <- (yaux - min(yaux)) / diff(range(yaux))
  dPPdx <- (sin(C * pi * x) * (C * pi) * (1 - x)^2 - (1 - cos(C * pi * x)) *
              (2 * (1 - x))) * ((1 - cos(C * pi * y)) * (1 - y)^2)
  dPPdy <- ((1 - cos(C * pi * x)) * (1 - x)^2) * (sin(C * pi * y) * (C * pi) *
                                                    (1 - y)^2 - (1 - cos(C * pi * y)) * (2 * (1 - y)))
  u <- -a * dPPdy
  v <- a * dPPdx

  return( data.frame(u = u, v = v) )
}


#' Radial advection field
#'
#' Provides an advection field corresponding to radial motion from or towards a specified reference point.
#'
#' @param spacepoints vector of coordinates (2 x d), where d is the number of locations/grid points
#' @param x0 x coordinate of the center of radial motion
#' @param y0 y coordinate of the center of radial motion
#' @param a parameter controlling the x component of radial velocity
#' @param b parameter controlling the y component of radial velocity
#'
#' @details # Note
#' \itemize{
#' \item if a > 0, b > 0: divergence from (x0, y0) (source point effect)
#' \item if a < 0, b < 0: convergence to (x0, y0) (sink effect)
#' }
#'
#' @name advectionFradial
#'
#' @import ggplot2
#' @export
#'
#' @references Papalexiou, S.M., Serinaldi, F., Porcu, E. (2021). Advancing
#' Space-Time Simulation of Random Fields: From Storms to Cyclones and Beyond.
#' Water Resources Research, 57, e2020WR029466, \doi{10.1029/2020WR029466}
#'
#' @examples
#' library(ggquiver)
#' library(ggplot2)
#'
#' ## specify coordinates
#' m = 25
#' aux <- seq(0, m - 1, length = m)
#' coord <- expand.grid(aux, aux)
#'
#' af <- advectionFradial(spacepoints = coord,
#'                         x0 = floor(m / 2),
#'                         y0 = floor(m / 2),
#'                         a = 3,
#'                         b = 2)
#'
#' ## visualize advection field
#' dta <- data.frame(lon = coord[ ,1], lat = coord[ ,2], u = af[ ,1], v = af[ ,2])
#' ggplot(dta, aes(x = lon, y = lat, u = u, v = v)) +
#' geom_quiver() +
#' theme_light()
#
advectionFradial <- function(spacepoints, x0, y0, a, b) {

  x <- spacepoints[,1] - x0
  y <- spacepoints[,2] - y0
  u <- a * x
  v <- b * y
  return( data.frame(u = u, v = v) )
}



#' Hyperbolic advection field
#'
#' Provides an advection field with hyperbolic trajectories.
#'
#' @param spacepoints vector of coordinates (2 x d), where d is the number of locations/grid points
#' @param x0 x coordinate of the center of hyperbola
#' @param y0 y coordinate of the center of hyperbola
#' @param a parameter controlling the x component of rotational velocity
#' @param b parameter controlling the y component of rotational velocity
#'
#' @details # Note
#' \itemize{
#' \item if a > 0, b > 0: toward bottom-left and top-right corner
#' \item if a < 0, b < 0: toward top-left and bottom-right corner
#' }
#'
#' @name advectionFhyperbolic
#'
#' @import ggplot2
#' @export
#'
#' @references Papalexiou, S.M., Serinaldi, F., Porcu, E. (2021). Advancing
#' Space-Time Simulation of Random Fields: From Storms to Cyclones and Beyond.
#' Water Resources Research, 57, e2020WR029466, \doi{10.1029/2020WR029466}
#'
#' @examples
#' library(ggquiver)
#' library(ggplot2)
#' ## specify coordinates
#' m = 25
#' aux <- seq(0, m - 1, length = m)
#' coord <- expand.grid(aux, aux)
#'
#' af <- advectionFhyperbolic(spacepoints = coord,
#'                            x0 = floor(m / 2),
#'                            y0 = floor(m / 2),
#'                            a = 3,
#'                            b = 2)
#'
#' ## visualize advection field
#' dta <- data.frame(lon = coord[ ,1], lat = coord[ ,2], u = af[ ,1], v = af[ ,2])
#' ggplot(dta, aes(x = lon, y = lat, u = u, v = v)) +
#' geom_quiver() +
#' theme_light()
#
advectionFhyperbolic <- function(spacepoints, x0, y0, a, b) {

  x <- spacepoints[,1] - x0
  y <- spacepoints[,2] - y0
  u <- a * y
  v <- b * x
  return( data.frame(u = u, v = v) )
}

