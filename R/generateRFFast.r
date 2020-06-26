#' Faster simulation of random fields with approximately separable spatiotemporal correlation structure
#'
#' For more details see section 6 in \href{https://doi.org/10.1029/2018WR023055}{Serinaldi and Kilsby (2018)}, and section 2.4 in \href{https://doi.org/10.1029/2019WR026331}{Papalexiou and Serinaldi (2020)}
#'
#' @param n number of fields (time steps) to simulate
#' @param spacepoints side length m of the square field (m x m)
#' @param margdist target marginal distribution of the field
#' @param margarg list of marginal distribution arguments
#' @param p0 probability zero
#' @param distbounds distribution bounds (default set to c(-Inf, Inf))
#' @param stcsid spatiotemporal correlation structure ID
#' @param stcsarg list of spatiotemporal correlation structure arguments
#' @param scalefactor factor specifying the distance between the centers of two pixels (default set to 1)
#' @param anisotropyarg list of arguments establishing spatial anisotropy. phi1 and phi2 control the stretch in two orthogonal directions (e.g., longitude and latitude) while the angle theta controls a counterclockwise rotation (default set to list(phi1 = 1, phi2 = 1 , theta = 0) for isotropic fields)
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
#' @details
#' \code{\link{generateRFFast}} provides a faster approach to RF simulation
#' compared to \code{\link{generateRF}} by exploiting circulant embedding
#' fast Fourier transformation.
#' However, this approach is feasible only for approximately
#' separable target spatiotemporal correlation functions (see section 6 in
#' \href{https://doi.org/10.1029/2018WR023055}{Serinaldi and Kilsby (2018)}).
#' \code{\link{generateRFFast}} comprises fitting and simulation in a single function.
#' Here, we give indicative CPU times for some settings, referring to a
#' Windows 10 Pro x64 laptop with Intel(R) Core(TM) i7-6700HQ CPU @ 2.60GHz,
#' 4-core, 8 logical processors, and 32GB RAM. \cr
#' CPU time:\cr
#' m = 50, n = 1000: ~58s \cr
#' m = 50, n = 10000: ~160s \cr
#' m = 100, n = 1000: ~2955s (~50min) \cr
#'
#' @examples
#'
#' sim <- generateRFFast(
#'     n = 50,
#'     spacepoints = 3,
#'     p0 = 0.7,
#'     margdist ='paretoII',
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
#' checkRF(sim,
#'           lags = 10,
#'           nfields = 49)
#'
generateRFFast <- function(n, spacepoints, margdist, margarg, p0, distbounds = c(-Inf, Inf), stcsid, stcsarg, scalefactor = 1, anisotropyarg = list(phi1 = 1, phi2 = 1 , theta = 0)) {

    DHMgenSj <- function(acvf) {
        MM <- length(acvf) - 1
        Sj <- Re(fft(c(acvf, rev(acvf[2:MM]))))[1:(MM + 1)]
        Sj <- abs(Sj)
        if (!(all(Sj >= 0)))
            stop("some of the S_j's are negative")
        Sj
    }

    DHMgenSim <- function(Sj, rn = NULL) {
        if (is.null(rn))
            rn <- rnorm(2 * length(Sj) - 2)
        M <- 2 * length(Sj) - 2
        N <- M/2
        Yj <- 0:N
        Yj[1] <- sqrt(Sj[1]) * rn[1]
        Yj[N + 1] <- sqrt(Sj[N + 1]) * rn[M]
        js <- 2:N
        Yj[js] <- sqrt(0.5 * Sj[js]) * complex(real = rn[2 *
                                                             (1:(N - 1))], imaginary = rn[2 * (1:(N - 1)) + 1])
        Re(fft(c(Yj, Conj(rev(Yj[js])))))[1:N]/sqrt(M)
    }

    if(class(spacepoints)[1] == "numeric"){
        mm <- spacepoints^2
        d <- seq(0, spacepoints - 1,length = spacepoints) * scalefactor
        dd <- expand.grid(d, d)
    } else {
        stop("spacepoints should be an integer")
    }

    phi1 <- anisotropyarg$phi1
    phi2 <- anisotropyarg$phi2
    theta <- anisotropyarg$theta
    The <- matrix(c(cos(theta) ,sin(theta), -sin(theta), cos(theta)), ncol = 2)
    Phi <- matrix(c(1/phi1, 0, 0, 1/phi2), ncol = 2)
    dd <- t(Phi %*%  The %*% t(dd))
    dis <- as.matrix(dist(dd, diag = T, upper = T))

    pnts <- actpnts(margdist = margdist, margarg = margarg, p0 = p0,distbounds = distbounds)
    actfpara <- fitactf(pnts)
    scf <- actf( acs(id=stcsarg$scfid, t = dis, scale = stcsarg$scfarg$scale,
                     shape = stcsarg$scfarg$shape),
                 b = actfpara$actfcoef[1], c = actfpara$actfcoef[2])
    if (!is.positive.definite(round(as.matrix(scf),6)))  scf = as.matrix(nearPD(scf, corr = T,conv.tol = 1e-04)$mat)

    tcf <- actf( acs(id=stcsarg$tcfid, t = 0:n, scale = stcsarg$tcfarg$scale,
                     shape = stcsarg$tcfarg$shape),
                 b = actfpara$actfcoef[1], c = actfpara$actfcoef[2])

    Sj <- DHMgenSj(tcf)
    rn <- rmvnorm(2 * length(Sj) - 2, sigma = scf)
    out <- matrix(NA, nrow = mm, ncol = n)
    for (i in 1:mm){
        out[i,] <- DHMgenSim(Sj = Sj, rn = rn[ ,i])
    }

    x <- as.matrix(out)
    x <- pnorm(x, mean(x), sd(x))
    x[x < p0] <- 0
    x[x > p0] <- do.call(paste0("q", margdist), args = c(list(p = (x[x > p0] - p0)/(1-p0)), margarg))

    d <- seq(0, max(dis), length = 200)
    u <- 0:n
    scf0 <- actf( acs(id=stcsarg$scfid, t = d, scale = stcsarg$scfarg$scale,
                      shape = stcsarg$scfarg$shape),
                  b = actfpara$actfcoef[1], c = actfpara$actfcoef[2])
    tcf0 <- actf( acs(id=stcsarg$tcfid, t = u, scale = stcsarg$tcfarg$scale,
                      shape = stcsarg$tcfarg$shape),
                  b = actfpara$actfcoef[1], c = actfpara$actfcoef[2])

    fout <- actfInv(tcf0%*%t(scf0),b = actfpara$actfcoef[1], c = actfpara$actfcoef[2])

    structure(.Data = t(x), class = c("matrix", "cosmosts"), STmodel=list(margdist = margdist, margarg = margarg,
                                                                          p0 = p0, stcs = fout, s.lags = d, t.lags=u, grid.lags = dis))
}



