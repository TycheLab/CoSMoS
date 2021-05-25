#' Anisotropy transformation
#'
#' Provides parametric functions that describe different types of planar deformation fields, including affine (rotation and stretching), and swirl-like deformation. For more details see Papalexiou et al.(2021) and references therein.
#'
#' @param id anisotropy type id (\code{affine}, \code{swirl}, and \code{wave})
#' @param ... additional arguments (vector of coordinates and parameters of the anisotropy transformations)
#'
#' @name anisotropyT
#'
#' @import ggplot2
#' @export
#'
#' @references Papalexiou, S. M., Serinaldi, F., Porcu, E. (2021). Advancing
#' Space-Time Simulation of Random Fields: From Storms to Cyclones and Beyond,
#' Water Resources Research, \doi{10.1029/2020WR029466}
#'
#' @examples
#'
#' library(CoSMoS)
#'
#' ## specify coordinates
#' m = 25
#' aux <- seq(0, m - 1, length = m)
#' coord <- expand.grid(aux, aux)
#'
#' ## get the anisotropy field
#' at1 <- anisotropyT('affine',
#'                  spacepoints = coord,
#'                  phi1 = 0.5,
#'                  phi2 = 2,
#'                  phi12 = 0,
#'                  theta = -pi/3)
#' at2 <- anisotropyT('swirl',
#'                  spacepoints = coord,
#'                  x0 = floor(m / 2),
#'                  y0 = floor(m / 2),
#'                  b = 10,
#'                  alpha = 1.5 * pi)
#' at3 <- anisotropyT('wave',
#'                  spacepoints = coord,
#'                  phi1 = 0.5,
#'                  phi2 = 2,
#'                  beta = 3,
#'                  theta = 0)
#'
#' ## visualize anisotropy field
#' aux = data.frame(lon = at2[ ,1], lat = at2[ ,2], id1 = rep(1:m, each = m), id2 = rep(1:m, m))
#' ggplot(aux, aes(x = lon, y = lat)) +
#' geom_path(aes(group = id1)) +
#' geom_path(aes(group = id2)) +
#' geom_point(col = 2) +
#' theme_light()
#'
anisotropyT <- function(id, ...) {

    .args <- list(...)
    do.call(paste0('anisotropyT', id), args = .args)
}


#' Anisotropy transformation
#'
#' Provides parametric functions that describe different types of planar deformation fields, including affine (rotation and stretching), and swirl-like deformation. For more details see Papalexiou et al.(2021) and references therein.
#'
#' @param id anisotropy type id
#' @param arglist list of additional arguments (vector of coordinates and parameters of the anisotropy transformations)
#'
#' @name anisotropyT2
#'
#' @references Papalexiou, S.M., Serinaldi, F., Porcu, E. (2021). Advancing
#' Space-Time Simulation of Random Fields: From Storms to Cyclones and Beyond.
#' Water Resources Research, 57, e2020WR029466, \doi{10.1029/2020WR029466}
#'
#' @keywords internal
#'
anisotropyT2 <- function(id, arglist) {

  do.call(paste0('anisotropyT', id), args = arglist)
}


#' Affine anisotropy transformation
#'
#' Affine anisotropy transformation.
#'
#' @param spacepoints vector of coordinates (2 x d), where d is the number of locations/grid points
#' @param phi1 stretching parameter along the x axis
#' @param phi2 stretching parameter along the y axis
#' @param phi12 shear effect
#' @param theta rotation angle
#'
#' @name anisotropyTaffine
#'
#' @import ggplot2
#' @export
#'
#' @references Allard, D., Senoussi, R., Porcu, E. (2016). Anisotropy Models for
#' Spatial Data. Mathematical Geosciences, 48(3), 305-328, \doi{10.1007/s11004-015-9594-x}
#' @references Papalexiou, S.M., Serinaldi, F., Porcu, E. (2021). Advancing
#' Space-Time Simulation of Random Fields: From Storms to Cyclones and Beyond.
#' Water Resources Research, 57, e2020WR029466, \doi{10.1029/2020WR029466}
#'
#' @examples
#'
#' ## specify coordinates
#' m = 25
#' aux <- seq(0, m - 1, length = m)
#' coord <- expand.grid(aux, aux)
#'
#' at <- anisotropyTaffine(spacepoints = coord,
#'                        phi1 = 0.5,
#'                        phi2 = 2,
#'                        phi12 = 0,
#'                        theta = -pi/3)
#'
#' ## visualize transformed coordinate system
#' aux = data.frame(lon = at[ ,1], lat = at[ ,2], id1 = rep(1:m, each = m), id2 = rep(1:m, m))
#' ggplot(aux, aes(x = lon, y = lat)) +
#' geom_path(aes(group = id1)) +
#' geom_path(aes(group = id2)) +
#' geom_point(col = 2) +
#' theme_light()
#'
anisotropyTaffine <- function(spacepoints, phi1, phi2, phi12, theta) {

  The <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), ncol = 2)
  Phi <- matrix(c(phi1, phi12, phi12, phi2), ncol = 2)
  return( t(Phi %*%  The %*% t(spacepoints)) )
}

#' Swirl anisotropy transformation
#'
#' Swirl anisotropy transformation.
#'
#' @param spacepoints vector of coordinates (2 x d), where d is the number of locations/grid points
#' @param x0 x coordinate of the center of the swirl deformation
#' @param y0 y coordinate of the center of the swirl deformation
#' @param b scaling parameter controlling the swirl deformation
#' @param alpha rotation angle
#'
#' @name anisotropyTswirl
#'
#' @import ggplot2
#' @export
#'
#' @references Ligas, M., Banas, M., Szafarczyk, A. (2019). A method for local
#' approximation of a planar deformation field. Reports on Geodesy and
#' Geoinformatics, 108(1), 1-8, \doi{10.2478/rgg-2019-0007}
#' @references Papalexiou, S.M., Serinaldi, F., Porcu, E. (2021). Advancing
#' Space-Time Simulation of Random Fields: From Storms to Cyclones and Beyond.
#' Water Resources Research, 57, e2020WR029466, \doi{10.1029/2020WR029466}
#'
#' @examples
#'
#' ## specify coordinates
#' m = 25
#' aux <- seq(0, m - 1, length = m)
#' coord <- expand.grid(aux, aux)
#'
#' at <- anisotropyTswirl(spacepoints = coord,
#'                       x0 = floor(m / 2),
#'                       y0 = floor(m / 2),
#'                       b = 10,
#'                       alpha = 1.5 * pi)
#'
#' ## visualize transformed coordinate system
#' aux = data.frame(lon = at[ ,1], lat = at[ ,2], id1 = rep(1:m, each = m), id2 = rep(1:m, m))
#' ggplot(aux, aes(x = lon, y = lat)) +
#' geom_path(aes(group = id1)) +
#' geom_path(aes(group = id2)) +
#' geom_point(col = 2) +
#' theme_light()
#
anisotropyTswirl <- function(spacepoints, x0, y0, b, alpha) {

  x <- spacepoints[,1]
  y <- spacepoints[,2]
  r <- sqrt((x - x0)^2 + (y - y0)^2)
  The <- alpha * exp(-(r / b)^2)
  xani <- (x - x0) * cos(The) - (y - y0) * sin(The) + x0
  yani <- (x - x0) * sin(The) + (y - y0) * cos(The) + y0
  return(cbind(xani, yani))
}

#' Wave anisotropy transformation
#'
#' Wave anisotropy transformation.
#'
#' @param spacepoints vector of coordinates (2 x d), where d is the number of locations/grid points
#' @param phi1 stretching parameter along the x axis
#' @param phi2 stretching parameter along the y axis
#' @param beta amplitude of sinusoidal wave
#' @param theta rotation angle
#'
#' @name anisotropyTwave
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
#' ## specify coordinates
#' m = 25
#' aux <- seq(0, m - 1, length = m)
#' coord <- expand.grid(aux, aux)
#'
#' at <- anisotropyTwave(spacepoints = coord,
#'                      phi1 = 0.5,
#'                      phi2 = 2,
#'                      beta = 3,
#'                      theta = 0)
#'
#' ## visualize transformed coordinate system
#' aux = data.frame(lon = at[ ,1], lat = at[ ,2], id1 = rep(1:m, each = m), id2 = rep(1:m, m))
#' ggplot(aux, aes(x = lon, y = lat)) +
#' geom_path(aes(group = id1)) +
#' geom_path(aes(group = id2)) +
#' geom_point(col = 2) +
#' theme_light()
#'
anisotropyTwave <- function(spacepoints, phi1, phi2, beta, theta) {

  x <- spacepoints[,1]
  y <- spacepoints[,2]
  xani <- phi1 * x * cos(theta) - y * sin(theta)
  yani <- x * sin(theta) + phi1 * y * cos(theta) + beta * sin(x)
  return(cbind(xani, yani))
}
