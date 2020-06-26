#' VAR model parameters to simulate parent Gaussian random fields
#'
#' Compute VAR model parameters to simulate parent Gaussian random fields with specified spatiotemporal correlation structure using the method described by \href{https://doi.org/10.1145/937332.937333}{Biller and Nelson (2003)}
#'
#' @param m side length of the square field (m x m)
#' @param p order of VAR(p) model
#' @param margdist target marginal distribution of the field
#' @param margarg list of marginal distribution arguments
#' @param p0 probability zero
#' @param distbounds distribution bounds (default set to c(-Inf, Inf))
#' @param stcsid spatiotemporal correlation structure ID
#' @param stcsarg list of spatiotemporal correlation structure arguments
#' @param scalefactor factor specifying the distance between the centers of two pixels (default set to 1)
#' @param anisotropyarg list of arguments establishing spatial anisotropy. phi1 and phi2 control the stretch in two orthogonal directions (e.g., longitude and latitude) while the angle theta controls a counterclockwise rotation (default set to list(phi1 = 1, phi2 = 1 , theta = 0) for isotropic fields)
#'
#' @name fitVARSTRF
#'
#' @importFrom Matrix nearPD
#'
#' @export
#'
#' @examples
#'
#' fit <- fitVARSTRF(
#'   m = 10,
#'   p = 1,
#'   margdist ='burrXII',
#'   margarg = list(scale = 3, shape1 = .9, shape2 = .2),
#'   p0 = 0.8,
#'   stcsid = "clayton",
#'   stcsarg = list(scfid = "weibull", tcfid = "weibull",
#'                  copulaarg = 2,
#'                  scfarg = list(scale = 20, shape = 0.7),
#'                 tcfarg = list(scale = 1.1, shape = 0.8))
#' )
#'
#' dim(fit$alpha)
#' dim(fit$res.cov)
#'
#' fit$m
#' fit$margarg
#' fit$margdist
#'
fitVARSTRF <- function(m, p, margdist, margarg, p0, distbounds = c(-Inf, Inf), stcsid, stcsarg, scalefactor = 1, anisotropyarg = list(phi1 = 1, phi2 = 1 , theta = 0)) {

    mm <- m^2
    d <- seq(0, m - 1,length = m) * scalefactor
    dd <- expand.grid(d, d)

    phi1 <- anisotropyarg$phi1
    phi2 <- anisotropyarg$phi2
    theta <- anisotropyarg$theta
    The <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), ncol = 2)
    Phi <- matrix(c(1/phi1, 0, 0, 1/phi2), ncol = 2)
    dd <- t(Phi %*%  The %*% t(dd))
    dis <- as.matrix(dist(dd, diag = T, upper = T))

    pnts <- actpnts(margdist = margdist, margarg = margarg, p0 = p0, distbounds = distbounds)
    actfpara <- fitactf(pnts)

    sigma.i <- array(NA, c(mm, mm, p+1))
    for (i in 1:(p+1)){
        stcsarg1 <- stcsarg
        stcsarg1$t <- i - 1
        stcsarg1$s <- dis
        sigma.i[,,i] <- actf( stcs2(stcsid, stcsarg1), b = actfpara$actfcoef[1], c = actfpara$actfcoef[2])
    }

    sigma.Z <- matrix(NA, p * mm, p * mm)
    sigma <- matrix(NA, mm, p * mm)
    for (i in 1:p) {
      for (j in 1:p) {
           if (j >= i) sigma.Z[seq(1 + mm * (i - 1), mm * i), seq(1 + mm * (j - 1), mm * j)] <-   sigma.i[,,j-i+1]
           else        sigma.Z[seq(1 + mm * (i - 1), mm * i), seq(1 + mm * (j - 1), mm * j)] <- t(sigma.i[,,i-j+1])
      }
      sigma[ , seq(1 + mm * (i - 1), mm * i)] <- sigma.i[ , , i+1]
    }

    alpha <- t(solve(sigma.Z, t(sigma)))

    sums <- matrix(0, mm, mm)
    for (i in 1:p) {
      sums <- sums  +  ( alpha[ , seq(1 + mm * (i - 1), mm * i)] %*% t(sigma.i[,,i+1]) )
    }
    sigma.u <- nearPD( sigma.i[ , , 1]  - sums )$mat

    d <- seq(0, max(dis), length = m)
    u <- 0:200
    ud <- expand.grid(u, d)
    stcsarg$t <- ud[ ,1]
    stcsarg$s <- ud[ ,2]
    fout <- (matrix(stcs2(stcsid, stcsarg) , nrow = length(u)))

    return(list(alpha = alpha, res.cov = sigma.u, m = m, p = p, margdist = margdist, margarg = margarg, p0 = p0,
                stcs = fout, s.lags = d, t.lags=u, grid.lags = dis, distbounds = distbounds))
}



