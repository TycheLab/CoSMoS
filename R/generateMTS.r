#' Simulation of multiple time series with given marginals and spatiotemporal properties
#'
#' Generates multiple time series with given marginals and spatiotemporal properties,
#' just provide (1) the output of \code{\link{fitVAR}} function, and (2) the number of time
#' steps to simulate.
#'
#' @param n number of fields (time steps) to simulate
#' @param STmodel list of arguments resulting from \code{\link{fitVAR}} function
#'
#' @name generateMTS
#'
#' @import mAr
#'
#' @export
#'
#' @details
#' Referring to the documentation of \code{\link{fitVAR}} for details on
#' computational complexity of the fitting algorithm, here we report indicative
#' simulation CPU times for some settings, assuming that the model parameters are
#' already evaluated.
#' CPU times refer to a Windows 10 Pro x64 laptop with
#' Intel(R) Core(TM) i7-6700HQ CPU @ 2.60GHz, 4-core, 8 logical processors, and 32GB RAM. \cr
#' CPU time:\cr
#' d = 900, p = 1, n = 1000: ~17s \cr
#' d = 900, p = 1, n = 10000: ~75s \cr
#' d = 900, p = 5, n = 100: ~280s \cr
#' d = 900, p = 5, n = 1000: ~302s \cr
#' d = 2500, p = 1, n = 1000 : ~160s \cr
#' d = 2500, p = 1, n = 10000 : ~570s \cr
#' where \eqn{d} denotes the number of spatial locations
#'
#' @examples
#' ## Simulation of a 4-dimensional vector with VAR(1) correlation structure
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
#' sim <- generateMTS(n = 100,
#'                      STmodel = fit)
#'
generateMTS <- function(n, STmodel) {

    mm <- STmodel$mm
    sim.out <- mAr.sim(w = rep(0, mm), A = STmodel$alpha, C = STmodel$res.cov, N = n)
    x <- as.matrix(sim.out)
    x <- pnorm(x,mean(x),sd(x))
    p0 <- STmodel$p0
    margdist <- STmodel$margdist
    margarg <- STmodel$margarg
    x[x < p0] <- 0
    x[x > p0] <- do.call(paste0("q", margdist), args = c(list(p = (x[x > p0] - p0)/(1-p0)), margarg))
    structure(.Data = x, class = c("matrix"),  STmodel = STmodel)
}
