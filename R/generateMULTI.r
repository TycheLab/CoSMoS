#' Simulation of spatiotemporal random vectors with given properties
#'
#' Generates spatiotemporal random vectors with given properties, just provide (1) the output of 'fitVARMULTI' function, and (2) the number of time steps to simulate
#'
#' @param n number of fields (time steps) to simulate
#' @param STmodel list of arguments resulting from 'fitVARMULTI' function
#'
#' @name generateMULTI
#'
#' @import mAr
#'
#' @export
#'
#' @examples
#'
#' coord <- cbind(runif(4)*30, runif(4)*30)
#'
#' fit <- fitVARMULTI(
#'   coord = coord,
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
#' sim <- generateMULTI(n = 100,
#'                      STmodel = fit)
#'
generateMULTI <- function(n, STmodel) {

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

