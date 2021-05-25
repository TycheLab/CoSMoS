#' VAR model parameters to simulate correlated parent Gaussian random vectors and fields
#'
#' Compute VAR model parameters to simulate parent Gaussian random vectors with specified spatiotemporal correlation structure using the method described by Biller and Nelson (2003).
#'
#' @param spacepoints it can be a numeric integer, which is interpreted as the side length m of the square field (m x m), or a matrix (d x 2) of coordinates (e.g. longitude and latitude) of d spatial locations (e.g. d gauge stations)
#' @param p order of VAR(p) model
#' @param margdist target marginal distribution of the field
#' @param margarg  list of marginal distribution arguments. Please consult the documentation of the selected marginal distribution indicated in the argument \code{margdist} for the list of required parameters
#' @param p0 probability zero
#' @param distbounds distribution bounds (default set to \code{c(-Inf, Inf)})
#' @param stcsid spatiotemporal correlation structure ID
#' @param stcsarg  list of spatiotemporal correlation structure arguments. Please consult the documentation of the selected spatiotemporal correlation structure indicated in the argument \code{stcsid} for the list of required parameters
#' @param scalefactor factor specifying the distance between the centers of two pixels (default set to 1)
#' @param anisotropyid spatial anisotropy ID (\code{affine} by default, \code{swirl} or \code{wave})
#' @param anisotropyarg list of arguments characterizing the spatial anisotropy according to the syntax of the function \code{\link{anisotropyT}}. Isotropic fields by default
#' @param advectionid advection field ID (\code{uniform} by default, \code{rotation}, \code{spiral}, \code{spiralCE}, \code{radial}, or \code{hyperbolic})
#' @param advectionarg list of arguments characterizing the advection field according to the syntax of the function \code{\link{advectionF}}. No advection by default
#'
#' @name fitVAR
#'
#' @importFrom Matrix nearPD
#' @importFrom pracma distmat
#'
#' @export
#'
#' @details
#' The fitting algorithm has \eqn{O(m*m)^3} complexity for a \eqn{(m*m)} field
#' or equivalently \eqn{O(d^3)} complexity for a \eqn{d}-dimensional vector.
#' Very large values of \eqn{(m*m)} (or \eqn{d}) and high order AR correlation
#' structures can be unpractical on standard machines.
#'
#' Here, we give indicative CPU times for some settings, referring to a
#' Windows 10 Pro x64 laptop with Intel(R) Core(TM) i7-6700HQ CPU @ 2.60GHz,
#' 4-core, 8 logical processors, and 32GB RAM. \cr:
#' CPU time:\cr
#' d = 100 or m = 10, p = 1: ~  0.4s \cr
#' d = 900 or m = 30, p = 1: ~  6.0s \cr
#' d = 900 or m = 30, p = 5: ~ 47.0s \cr
#' d = 2500 or m = 50, p = 1: ~100.0s
#'
#' @note
#' While all the advection types can be applied to isotropic random fields,
#' anisotropic random fields require more care. We suggest combining \code{affine}
#' anysotropy with \code{uniform} advection, and \code{swirl} anisotropy
#' with \code{rotation} or \code{spiral} advection with the same rotation center.
#'
#' @references Biller, B., Nelson, B.L. (2003). Modeling and generating multivariate
#' time-series input processes using a vector autoregressive technique.
#' ACM Trans. Model. Comput. Simul. 13(3), 211-237, \doi{10.1145/937332.937333}
#' @references Papalexiou, S.M. (2018). Unified theory for stochastic modelling of
#' hydroclimatic processes: Preserving marginal distributions, correlation structures,
#' and intermittency. Advances in Water Resources, 115, 234-252,
#' \doi{10.1016/j.advwatres.2018.02.013}
#' @references Papalexiou, S.M.,  Serinaldi, F. (2020). Random Fields Simplified:
#' Preserving Marginal Distributions, Correlations, and Intermittency,
#' With Applications From Rainfall to Humidity. Water Resources Research, 56(2),
#' e2019WR026331, \doi{10.1029/2019WR026331}
#' @references Papalexiou, S.M., Serinaldi, F., Porcu, E. (2021). Advancing
#' Space-Time Simulation of Random Fields: From Storms to Cyclones and Beyond.
#' Water Resources Research, 57, e2020WR029466, \doi{10.1029/2020WR029466}
#'
#' @examples
#' ## for multivariate simulation
#' coord <- cbind(runif(4)*30, runif(4)*30)
#'
#' fit <- fitVAR(
#'   spacepoints = coord,
#'   p = 1,
#'   margdist ='burrXII',
#'   margarg = list(scale = 3,
#'                  shape1 = .9,
#'                  shape2 = .2),
#'   p0 = 0.8,
#'   stcsid = "clayton",
#'   stcsarg = list(scfid = "weibull",
#'                  tcfid = "weibull",
#'                  copulaarg = 2,
#'                  scfarg = list(scale = 20,
#'                                shape = 0.7),
#'                  tcfarg = list(scale = 1.1,
#'                                shape = 0.8))
#' )
#'
#' dim(fit$alpha)
#' dim(fit$res.cov)
#'
#' fit$m
#' fit$margarg
#' fit$margdist
#'
#' ## for random fields simulation
#' fit <- fitVAR(
#'   spacepoints = 10,
#'   p = 1,
#'   margdist ='burrXII',
#'   margarg = list(scale = 3, shape1 = .9, shape2 = .2),
#'   p0 = 0.8,
#'   stcsid = "clayton",
#'   stcsarg = list(scfid = "weibull", tcfid = "weibull",
#'                  copulaarg = 2,
#'                  scfarg = list(scale = 20, shape = 0.7),
#'                  tcfarg = list(scale = 1.1, shape = 0.8))
#' )
#'
#' dim(fit$alpha)
#' dim(fit$res.cov)
#'
#' fit$m
#' fit$margarg
#' fit$margdist
#'
fitVAR <- function(spacepoints, p, margdist, margarg, p0, distbounds = c(-Inf, Inf), stcsid, stcsarg, scalefactor = 1, anisotropyid = "affine", anisotropyarg = list(phi1 = 1, phi2 = 1, phi12 = 0, theta = 0), advectionid = "uniform", advectionarg = list(u = 0, v = 0)) {

  if(class(spacepoints)[1] == "numeric"){
    mm <- spacepoints^2
    d <- seq(0, spacepoints - 1, length = spacepoints) * scalefactor
    dd0 <- expand.grid(d, d)
  } else {
    if(class(spacepoints)[1] %in% c("matrix", "array")){
      mm <- dim(spacepoints)[1]
      dd0 <- spacepoints
    } else stop("spacepoints should be an integer or a matrix (dx2)")
  }

  anisotropyarg$spacepoints <- dd0
  dd <- anisotropyT2(anisotropyid, anisotropyarg)
  dis <- dist(dd)

  advectionarg$spacepoints <- dd0
  uv <- advectionF2(advectionid, advectionarg)

  if (anisotropyid == "swirl") {
    advectionarg$spacepoints <- cbind(anisotropyarg$x0, anisotropyarg$y0)
    uv1 <- advectionF2(advectionid, advectionarg)
    x0 =  anisotropyarg$x0
    y0 =  anisotropyarg$y0
  }

  pnts <- actpnts(margdist = margdist, margarg = margarg, p0 = p0, distbounds = distbounds)
  actfpara <- fitactf(pnts)

  sigma.i <- array(NA,c(mm, mm, p+1))
  for (i in 1:(p + 1)){
    stcsarg1 <- stcsarg
    stcsarg1$t <- i - 1
    anisotropyarg$spacepoints <- dd0 + uv * stcsarg1$t
    if (anisotropyid == "swirl") {
      anisotropyarg$x0 =  x0 + uv1[ ,1] * stcsarg1$t
      anisotropyarg$y0 =  y0 + uv1[ ,2] * stcsarg1$t
    }
    dd1 <- anisotropyT2(anisotropyid, anisotropyarg)
    stcsarg1$s <- distmat(dd,dd1)
    sigma.i[ , ,i] <- actf( stcs2(stcsid, stcsarg1), b = actfpara$actfcoef[1], c = actfpara$actfcoef[2])
  }

  sigma.Z <- matrix(NA,p * mm,p * mm)
  sigma <- matrix(NA, mm, p * mm)
  for (i in 1:p) {
    for (j in 1:p) {
      if (j >= i) sigma.Z[seq(1 + mm * (i - 1), mm * i), seq(1 + mm * (j - 1), mm * j)] =   sigma.i[,,j-i+1]
      else        sigma.Z[seq(1 + mm * (i - 1), mm * i), seq(1 + mm * (j - 1), mm * j)] = t(sigma.i[,,i-j+1])
    }
    sigma[ , seq(1 + mm * (i - 1), mm * i)] = sigma.i[ , , i+1]
  }

  alpha <- t(solve(sigma.Z, t(sigma)))

  sums <- matrix(0, mm, mm)
  for (i in 1:p) {
    sums <- sums  +  ( alpha[ , seq(1 + mm * (i - 1), mm * i)] %*% t(sigma.i[,,i+1]) )
  }
  sigma.u <- nearPD( sigma.i[ , , 1]  - sums )$mat

  d <- seq(0, max(dis), length = mm)
  u <- 0:200
  ud <- expand.grid(u, d)
  stcsarg$t <- ud[ ,1]
  stcsarg$s <- ud[ ,2]
  fout <- (matrix(stcs2(stcsid, stcsarg) , nrow = length(u)))

  return(list(alpha = alpha, res.cov = sigma.u, spacepoints = spacepoints, mm = mm, p = p, margdist = margdist, margarg = margarg, p0 = p0,
              stcs = fout, s.lags = d, t.lags=u, grid.lags = dis, distbounds = distbounds))
}


