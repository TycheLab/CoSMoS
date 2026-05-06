# Internal helper shared by generateRF() and generateMTS()
# Runs the VAR simulation and applies the marginal back-transformation.
# Not exported.
simulate_var <- function(n, STmodel) {
  mm      <- STmodel$mm
  sim.out <- mAr.sim(w = rep(0, mm),
                     A = STmodel$alpha,
                     C = STmodel$res.cov,
                     N = n)
  x <- as.matrix(sim.out)

  if (STmodel$dsid == "gauss") {
    x <- pnorm(x, mean(x), sd(x))
  }
  if (STmodel$dsid == "student") {
    nu <- STmodel$dsarg
    S  <- rchisq(n, nu)
    x  <- pt(x * sqrt(nu / S), nu)
  }
  if (STmodel$dsid == "bardossy") {
    m <- STmodel$dsarg
    set.seed(666)
    x <- apply((x - mean(x) + m)^2, 2,
                function(x) jitter(rank(x) / (n + 1), amount = 1 / (n + 1)))
  }
  if (STmodel$dsid == "bardossyF") {
    m <- STmodel$dsarg
    set.seed(666)
    x <- 1 - apply((x - mean(x) + m)^2, 2,
                   function(x) jitter(rank(x) / (n + 1), amount = 1 / (n + 1)))
  }

  p0       <- STmodel$p0
  margdist <- STmodel$margdist
  margarg  <- STmodel$margarg

  x[x <= p0] <- 0
  x[x > p0]  <- do.call(paste0("q", margdist),
                         args = c(list(p = (x[x > p0] - p0) / (1 - p0)),
                                  margarg))

  structure(.Data = x, class = "matrix", STmodel = STmodel)
}


#' Simulation of random fields with given marginals and spatiotemporal properties
#'
#' Generates a random field with given marginals and spatiotemporal properties.
#' Provide (1) the output of \code{\link{fitVAR}} and (2) the number of time
#' steps to simulate.
#'
#' @param n number of fields (time steps) to simulate
#' @param STmodel list of arguments from \code{\link{fitVAR}}
#'
#' @return A matrix of class \code{"matrix"} with attribute \code{STmodel}.
#'   Rows correspond to spatial locations and columns to time steps.
#'
#' @seealso \code{\link{fitVAR}}, \code{\link{checkRF}}, \code{\link{generateMTS}}
#'
#' @name generateRF
#'
#' @import mAr
#'
#' @export
#'
#' @details
#' Referring to the documentation of \code{\link{fitVAR}} for details on
#' computational complexity, here we report indicative simulation CPU times,
#' assuming model parameters are already evaluated. CPU times refer to a
#' Windows 10 Pro x64 laptop with Intel(R) Core(TM) i7-6700HQ CPU @ 2.60GHz,
#' 4-core, 8 logical processors, and 32 GB RAM. \cr
#' CPU time: \cr
#' m = 30, p = 1, n = 1000: ~17s \cr
#' m = 30, p = 1, n = 10000: ~75s \cr
#' m = 30, p = 5, n = 100: ~280s \cr
#' m = 30, p = 5, n = 1000: ~302s \cr
#' m = 50, p = 1, n = 1000: ~160s \cr
#' m = 50, p = 1, n = 10000: ~570s \cr
#' where m denotes the side length of a square field (m x m).
#'
#' @examples
#' ## The example below simulates few random fields of size 10x10 with AR(1)
#' ## temporal correlation for illustration. For reliable performance assessment
#' ## generate a larger number of fields (e.g. 100 or more) of size ~30x30.
#' ## See 'Details' for running times with different settings.
#'
#' fit <- fitVAR(
#'   spacepoints = 10,
#'   p = 1,
#'   margdist = "burrXII",
#'   margarg = list(scale = 3, shape1 = .9, shape2 = .2),
#'   p0 = 0.8,
#'   stcsid = "clayton",
#'   stcsarg = list(scfid = "weibull", tcfid = "weibull",
#'                  copulaarg = 2,
#'                  scfarg = list(scale = 20, shape = 0.7),
#'                  tcfarg = list(scale = 1.1, shape = 0.8))
#' )
#'
#' sim <- generateRF(n = 12, STmodel = fit)
#' checkRF(sim, lags = 10, nfields = 12)
#'
generateRF <- function(n, STmodel) {
  simulate_var(n, STmodel)
}


#' Simulation of multiple time series with given marginals and spatiotemporal
#' properties
#'
#' Generates multiple time series with given marginals and spatiotemporal
#' properties. Provide (1) the output of \code{\link{fitVAR}} and (2) the number
#' of time steps to simulate.
#'
#' @param n number of time steps to simulate
#' @param STmodel list of arguments from \code{\link{fitVAR}}
#'
#' @return A matrix of class \code{"matrix"} with attribute \code{STmodel}.
#'   Rows correspond to time steps and columns to spatial locations.
#'
#' @seealso \code{\link{fitVAR}}, \code{\link{generateRF}}, \code{\link{generateMTSFast}}
#'
#' @name generateMTS
#'
#' @import mAr
#'
#' @export
#'
#' @details
#' Referring to the documentation of \code{\link{fitVAR}} for details on
#' computational complexity, here we report indicative simulation CPU times,
#' assuming model parameters are already evaluated. CPU times refer to a
#' Windows 10 Pro x64 laptop with Intel(R) Core(TM) i7-6700HQ CPU @ 2.60GHz,
#' 4-core, 8 logical processors, and 32 GB RAM. \cr
#' CPU time: \cr
#' d = 900, p = 1, n = 1000: ~17s \cr
#' d = 900, p = 1, n = 10000: ~75s \cr
#' d = 900, p = 5, n = 100: ~280s \cr
#' d = 900, p = 5, n = 1000: ~302s \cr
#' d = 2500, p = 1, n = 1000: ~160s \cr
#' d = 2500, p = 1, n = 10000: ~570s \cr
#' where \eqn{d} denotes the number of spatial locations.
#'
#' @examples
#' ## Simulation of a 4-dimensional vector with VAR(1) correlation structure
#' coord <- cbind(runif(4) * 30, runif(4) * 30)
#'
#' fit <- fitVAR(
#'   spacepoints = coord,
#'   p = 1,
#'   margdist = "burrXII",
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
#' sim <- generateMTS(n = 100, STmodel = fit)
#'
generateMTS <- function(n, STmodel) {
  simulate_var(n, STmodel)
}


#' Faster simulation of random fields with approximately separable
#' spatiotemporal correlation structure
#'
#' For more details see section 6 in Serinaldi and Kilsby (2018) and section 2.4
#' in Papalexiou and Serinaldi (2020).
#'
#' @param n number of fields (time steps) to simulate
#' @param spacepoints side length m of the square field (m x m)
#' @param margdist target marginal distribution of the field
#' @param margarg list of marginal distribution arguments; consult the
#'   documentation of the selected distribution for the required parameters
#' @param p0 probability zero
#' @param distbounds distribution bounds (default \code{c(-Inf, Inf)})
#' @param stcsarg list of spatiotemporal correlation structure arguments;
#'   consult the documentation of the selected structure for required parameters
#' @param scalefactor factor specifying the distance between pixel centres
#'   (default 1)
#' @param anisotropyid spatial anisotropy ID (\code{"affine"} by default;
#'   \code{"swirl"} or \code{"wave"} also available)
#' @param anisotropyarg list of arguments for \code{\link{anisotropyT}};
#'   isotropic fields by default
#' @param dsid dependence structure ID (\code{"gauss"} by default;
#'   \code{"student"}, \code{"bardossy"}, or \code{"bardossyF"})
#' @param dsarg argument for the dependence structure: \code{NULL} for
#'   \code{"gauss"}, degrees of freedom for \code{"student"}, or parameter
#'   \code{m} in \eqn{(-\infty, \infty)} for \code{"bardossy"}
#'
#' @return A matrix of class \code{c("matrix", "cosmosts")} with attribute
#'   \code{STmodel} containing the fitted model components.
#'
#' @seealso \code{\link{generateRF}}, \code{\link{generateMTSFast}},
#'   \code{\link{checkRF}}, \code{\link{fitVAR}}
#'
#' @name generateRFFast
#'
#' @import mvtnorm
#'
#' @importFrom Matrix nearPD
#' @importFrom matrixcalc is.positive.definite
#'
#' @export
#'
#' @references Serinaldi, F., Kilsby, C.G. (2018). Unsurprising Surprises:
#' The Frequency of Record-breaking and Overthreshold Hydrological Extremes Under
#' Spatial and Temporal Dependence. Water Resources Research, 54(9), 6460-6487,
#' \doi{10.1029/2018WR023055}
#' @references Papalexiou, S.M., Serinaldi, F. (2020). Random Fields Simplified:
#' Preserving Marginal Distributions, Correlations, and Intermittency,
#' With Applications From Rainfall to Humidity. Water Resources Research, 56(2),
#' e2019WR026331, \doi{10.1029/2019WR026331}
#'
#' @details
#' \code{\link{generateRFFast}} provides faster RF simulation than
#' \code{\link{generateRF}} by exploiting circulant-embedding fast Fourier
#' transformation. This approach is feasible only for approximately separable
#' target spatiotemporal correlation functions.
#' \code{\link{generateRFFast}} combines fitting and simulation in a single call.
#' Indicative CPU times (Windows 10 Pro x64, Intel Core i7-6700HQ, 32 GB RAM): \cr
#' m = 50, n = 1000: ~58s \cr
#' m = 50, n = 10000: ~160s \cr
#' m = 100, n = 1000: ~2955s (~50 min) \cr
#'
#' @examples
#'
#' sim <- generateRFFast(
#'     n = 50,
#'     spacepoints = 3,
#'     p0 = 0.7,
#'     margdist = "paretoII",
#'     margarg = list(scale = 1,
#'                    shape = .3),
#'     stcsarg = list(scfid = "weibull",
#'                    tcfid = "weibull",
#'                    scfarg = list(scale = 20,
#'                                  shape = 0.7),
#'                    tcfarg = list(scale = 1.1,
#'                                  shape = 0.8))
#' )
#'
#' checkRF(sim, lags = 10, nfields = 49)
#'
generateRFFast <- function(n, spacepoints, margdist, margarg, p0,
                           distbounds = c(-Inf, Inf), stcsarg,
                           scalefactor = 1,
                           anisotropyid = "affine",
                           anisotropyarg = list(phi1 = 1, phi2 = 1,
                                                phi12 = 0, theta = 0),
                           dsid = "gauss", dsarg = NULL) {

  if (!is.numeric(spacepoints))
    stop("spacepoints should be an integer")

  mm <- spacepoints^2
  d  <- seq(0, spacepoints - 1, length = spacepoints) * scalefactor
  dd <- expand.grid(d, d)

  anisotropyarg$spacepoints <- dd
  dd  <- anisotropyT2(anisotropyid, anisotropyarg)
  dis <- as.matrix(dist(dd, diag = TRUE, upper = TRUE))

  if (dsid %in% c("gauss", "student")) {
    pnts <- actpnts(margdist = margdist, margarg = margarg,
                    p0 = p0, distbounds = distbounds)
  }
  if (dsid %in% c("bardossy", "bardossyF")) {
    pnts <- actpntsB6(margdist = margdist, margarg = margarg,
                      m = dsarg, p0 = p0)
  }
  actfpara <- fitactf(pnts)

  scf <- actf(acs(id = stcsarg$scfid, t = dis,
                  scale = stcsarg$scfarg$scale,
                  shape = stcsarg$scfarg$shape),
              b = actfpara$actfcoef[1],
              c = actfpara$actfcoef[2])
  if (!is.positive.definite(round(as.matrix(scf), 6)))
    scf <- as.matrix(nearPD(scf, corr = TRUE, conv.tol = 1e-04)$mat)

  tcf <- actf(acs(id = stcsarg$tcfid, t = 0:n,
                  scale = stcsarg$tcfarg$scale,
                  shape = stcsarg$tcfarg$shape),
              b = actfpara$actfcoef[1],
              c = actfpara$actfcoef[2])

  Sj  <- DHMgenSj(tcf)
  rn  <- rmvnorm(2 * length(Sj) - 2, sigma = scf)
  out <- matrix(NA, nrow = mm, ncol = n)
  for (i in 1:mm) {
    out[i, ] <- DHMgenSim(Sj = Sj, rn = rn[, i])
  }

  x <- t(as.matrix(out))

  if (dsid == "gauss") {
    x <- pnorm(x, mean(x), sd(x))
  }
  if (dsid == "student") {
    nu <- dsarg
    S  <- rchisq(n, nu)
    x  <- pt(x * sqrt(nu / S), nu)
  }
  if (dsid == "bardossy") {
    m <- dsarg
    set.seed(666)
    x <- apply((x - mean(x) + m)^2, 2,
               function(x) jitter(rank(x) / (n + 1), amount = 1 / (n + 1)))
  }
  if (dsid == "bardossyF") {
    m <- dsarg
    set.seed(666)
    x <- 1 - apply((x - mean(x) + m)^2, 2,
                   function(x) jitter(rank(x) / (n + 1), amount = 1 / (n + 1)))
  }

  x[x <= p0] <- 0
  x[x > p0]  <- do.call(paste0("q", margdist),
                         args = c(list(p = (x[x > p0] - p0) / (1 - p0)),
                                  margarg))

  d    <- seq(0, max(dis), length = 200)
  u    <- 0:n
  scf0 <- actf(acs(id = stcsarg$scfid, t = d,
                   scale = stcsarg$scfarg$scale,
                   shape = stcsarg$scfarg$shape),
               b = actfpara$actfcoef[1],
               c = actfpara$actfcoef[2])
  tcf0 <- actf(acs(id = stcsarg$tcfid, t = u,
                   scale = stcsarg$tcfarg$scale,
                   shape = stcsarg$tcfarg$shape),
               b = actfpara$actfcoef[1],
               c = actfpara$actfcoef[2])

  fout <- actfInv(tcf0 %*% t(scf0),
                  b = actfpara$actfcoef[1],
                  c = actfpara$actfcoef[2])

  structure(.Data = x,
            class = c("matrix", "cosmosts"),
            STmodel = list(margdist = margdist,
                           margarg  = margarg,
                           p0       = p0,
                           stcs     = fout,
                           s.lags   = d,
                           t.lags   = u,
                           grid.lags = dis))
}


#' Faster simulation of multiple time series with approximately separable
#' spatiotemporal correlation structure
#'
#' For more details see section 6 in Serinaldi and Kilsby (2018) and section 2.4
#' in Papalexiou and Serinaldi (2020).
#'
#' @param n number of time steps to simulate
#' @param spacepoints matrix (d x 2) of coordinates (e.g. longitude and
#'   latitude) for d spatial locations (e.g. gauge stations)
#' @param margdist target marginal distribution
#' @param margarg list of marginal distribution arguments; consult the
#'   documentation of the selected distribution for the required parameters
#' @param p0 probability zero
#' @param distbounds distribution bounds (default \code{c(-Inf, Inf)})
#' @param stcsarg list of spatiotemporal correlation structure arguments;
#'   consult the documentation of the selected structure for required parameters
#' @param scalefactor factor specifying the distance between pixel centres
#'   (default 1)
#' @param anisotropyid spatial anisotropy ID (\code{"affine"} by default;
#'   \code{"swirl"} or \code{"wave"} also available)
#' @param anisotropyarg list of arguments for \code{\link{anisotropyT}};
#'   isotropic fields by default
#' @param dsid dependence structure ID (\code{"gauss"} by default;
#'   \code{"student"}, \code{"bardossy"}, or \code{"bardossyF"})
#' @param dsarg argument for the dependence structure: \code{NULL} for
#'   \code{"gauss"}, degrees of freedom for \code{"student"}, or parameter
#'   \code{m} in \eqn{(-\infty, \infty)} for \code{"bardossy"}
#'
#' @return A matrix of class \code{c("matrix", "cosmosts")} with attribute
#'   \code{STmodel} containing the fitted model components.
#'
#' @seealso \code{\link{generateMTS}}, \code{\link{generateRFFast}},
#'   \code{\link{fitVAR}}
#'
#' @name generateMTSFast
#'
#' @import mvtnorm
#'
#' @importFrom Matrix nearPD
#' @importFrom matrixcalc is.positive.definite
#'
#' @export
#'
#' @references Serinaldi, F., Kilsby, C.G. (2018). Unsurprising Surprises:
#' The Frequency of Record-breaking and Overthreshold Hydrological Extremes Under
#' Spatial and Temporal Dependence. Water Resources Research, 54(9), 6460-6487,
#' \doi{10.1029/2018WR023055}
#' @references Papalexiou, S.M., Serinaldi, F. (2020). Random Fields Simplified:
#' Preserving Marginal Distributions, Correlations, and Intermittency,
#' With Applications From Rainfall to Humidity. Water Resources Research, 56(2),
#' e2019WR026331, \doi{10.1029/2019WR026331}
#'
#' @details
#' \code{\link{generateMTSFast}} provides faster multivariate simulation than
#' \code{\link{generateMTS}} by exploiting circulant-embedding fast Fourier
#' transformation. This approach is feasible only for approximately separable
#' target spatiotemporal correlation functions.
#' \code{\link{generateMTSFast}} combines fitting and simulation in a single call.
#' Indicative CPU times (Windows 10 Pro x64, Intel Core i7-6700HQ, 32 GB RAM): \cr
#' d = 2500, n = 1000: ~58s \cr
#' d = 2500, n = 10000: ~160s \cr
#' d = 10000, n = 1000: ~2955s (~50 min) \cr
#' where \eqn{d} denotes the number of spatial locations.
#'
#' @examples
#' coord <- cbind(runif(4) * 30, runif(4) * 30)
#'
#' sim <- generateMTSFast(
#'     n = 50,
#'     spacepoints = coord,
#'     p0 = 0.7,
#'     margdist = "paretoII",
#'     margarg = list(scale = 1,
#'                    shape = .3),
#'     stcsarg = list(scfid = "weibull",
#'                    tcfid = "weibull",
#'                    scfarg = list(scale = 20,
#'                                  shape = 0.7),
#'                    tcfarg = list(scale = 1.1,
#'                                  shape = 0.8))
#' )
#'
generateMTSFast <- function(n, spacepoints, margdist, margarg, p0,
                            distbounds = c(-Inf, Inf), stcsarg,
                            scalefactor = 1,
                            anisotropyid = "affine",
                            anisotropyarg = list(phi1 = 1, phi2 = 1,
                                                 phi12 = 0, theta = 0),
                            dsid = "gauss", dsarg = NULL) {

  if (!(is.matrix(spacepoints) || is.array(spacepoints)))
    stop("spacepoints should be a matrix (d x 2)")

  mm <- nrow(spacepoints)
  dd <- spacepoints

  anisotropyarg$spacepoints <- dd
  dd  <- anisotropyT2(anisotropyid, anisotropyarg)
  dis <- as.matrix(dist(dd, diag = TRUE, upper = TRUE))

  if (dsid %in% c("gauss", "student")) {
    pnts <- actpnts(margdist = margdist, margarg = margarg,
                    p0 = p0, distbounds = distbounds)
  }
  if (dsid %in% c("bardossy", "bardossyF")) {
    pnts <- actpntsB6(margdist = margdist, margarg = margarg,
                      m = dsarg, p0 = p0)
  }
  actfpara <- fitactf(pnts)

  scf <- actf(acs(id = stcsarg$scfid, t = dis,
                  scale = stcsarg$scfarg$scale,
                  shape = stcsarg$scfarg$shape),
              b = actfpara$actfcoef[1],
              c = actfpara$actfcoef[2])
  if (!is.positive.definite(round(as.matrix(scf), 6)))
    scf <- as.matrix(nearPD(scf, corr = TRUE, conv.tol = 1e-04)$mat)

  tcf <- actf(acs(id = stcsarg$tcfid, t = 0:n,
                  scale = stcsarg$tcfarg$scale,
                  shape = stcsarg$tcfarg$shape),
              b = actfpara$actfcoef[1],
              c = actfpara$actfcoef[2])

  Sj  <- DHMgenSj(tcf)
  rn  <- rmvnorm(2 * length(Sj) - 2, sigma = scf)
  out <- matrix(NA, nrow = mm, ncol = n)
  for (i in 1:mm) {
    out[i, ] <- DHMgenSim(Sj = Sj, rn = rn[, i])
  }

  x <- t(as.matrix(out))

  if (dsid == "gauss") {
    x <- pnorm(x, mean(x), sd(x))
  }
  if (dsid == "student") {
    nu <- dsarg
    S  <- rchisq(n, nu)
    x  <- pt(x * sqrt(nu / S), nu)
  }
  if (dsid == "bardossy") {
    m <- dsarg
    set.seed(666)
    x <- apply((x - mean(x) + m)^2, 2,
               function(x) jitter(rank(x) / (n + 1), amount = 1 / (n + 1)))
  }
  if (dsid == "bardossyF") {
    m <- dsarg
    set.seed(666)
    x <- 1 - apply((x - mean(x) + m)^2, 2,
                   function(x) jitter(rank(x) / (n + 1), amount = 1 / (n + 1)))
  }

  x[x <= p0] <- 0
  x[x > p0]  <- do.call(paste0("q", margdist),
                         args = c(list(p = (x[x > p0] - p0) / (1 - p0)),
                                  margarg))

  d    <- seq(0, max(dis), length = 200)
  u    <- 0:n
  scf0 <- actf(acs(id = stcsarg$scfid, t = d,
                   scale = stcsarg$scfarg$scale,
                   shape = stcsarg$scfarg$shape),
               b = actfpara$actfcoef[1],
               c = actfpara$actfcoef[2])
  tcf0 <- actf(acs(id = stcsarg$tcfid, t = u,
                   scale = stcsarg$tcfarg$scale,
                   shape = stcsarg$tcfarg$shape),
               b = actfpara$actfcoef[1],
               c = actfpara$actfcoef[2])

  fout <- actfInv(tcf0 %*% t(scf0),
                  b = actfpara$actfcoef[1],
                  c = actfpara$actfcoef[2])

  structure(.Data = x,
            class = c("matrix", "cosmosts"),
            STmodel = list(margdist  = margdist,
                           margarg   = margarg,
                           p0        = p0,
                           stcs      = fout,
                           s.lags    = d,
                           t.lags    = u,
                           grid.lags = dis))
}
