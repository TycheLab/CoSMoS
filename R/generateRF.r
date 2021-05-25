#' Simulation of random field with given marginals and spatiotemporal properties
#'
#' Generates random field with given marginals and spatiotemporal properties,
#' just provide (1) the output of \code{\link{fitVAR}} function, and (2) the number of time
#' steps to simulate.
#'
#' @param n number of fields (time steps) to simulate
#' @param STmodel list of arguments resulting from \code{\link{fitVAR}} function
#'
#' @name generateRF
#'
#' @import mAr
#'
#' @export
#'
#' @details
#' Referring to the documentation of \code{\link{fitVAR}} for details on
#' computational complexity of the fitting algorithm, here we report indicative
#' simulation CPU times for some settings, assuming that the model parameters are
#' already evaluated. CPU times refer to a Windows 10 Pro x64 laptop with
#' Intel(R) Core(TM) i7-6700HQ CPU @ 2.60GHz, 4-core, 8 logical processors, and 32GB RAM. \cr
#' CPU time:\cr
#' m = 30, p = 1, n = 1000: ~17s \cr
#' m = 30, p = 1, n = 10000: ~75s \cr
#' m = 30, p = 5, n = 100: ~280s \cr
#' m = 30, p = 5, n = 1000: ~302s \cr
#' m = 50, p = 1, n = 1000 : ~160s \cr
#' m = 50, p = 1, n = 10000 : ~570s
#' where m denotes the side length of a square field (mxm)
#'
#' @examples
#' ## The example below refers to the simulation of few random fields of
#' ## size 10x10 with AR(1) temporal correlation for the sake of illustration.
#' ## For a more effective visualization and reliable performance assessment,
#' ## we suggest to generate a larger number of fields (e.g. 100 or more)
#' ## of size about 30X30.
#' ## See section 'Details' for additional information on running times
#' ## with different settings.
#'
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
#'                 tcfarg = list(scale = 1.1, shape = 0.8))
#' )
#'
#' sim <- generateRF(n = 12,
#'                     STmodel = fit)
#' checkRF(sim,
#'           lags = 10,
#'           nfields = 12)
#'
generateRF <- function(n, STmodel) {

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
