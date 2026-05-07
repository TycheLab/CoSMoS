#' Check generated time series
#'
#' Compares sample statistics of generated time series against theoretically
#' expected values.
#'
#' @param TS a \code{cosmosts} object from \code{\link{generateTS}}, or a list
#'   of numeric vectors, or a single numeric vector
#' @inheritParams moments
#'
#' @return An object of class \code{c("checkTS", "matrix")} with rows
#'   \code{"expected"} and one row per simulated series, and columns for
#'   \code{mean}, \code{sd}, \code{skew}, \code{p0}, \code{acf_t1},
#'   \code{acf_t2}, \code{acf_t3}. Attributes \code{margdist}, \code{margarg},
#'   and \code{p0} are attached for use by \code{\link{plot.checkTS}}.
#'
#' @seealso \code{\link{generateTS}}, \code{\link{plot.checkTS}},
#'   \code{\link{moments}}
#'
#' @export
#'
#' @examples
#'
#' library(CoSMoS)
#'
#' x <- generateTS(margdist = "burrXII",
#'                 margarg = list(scale = 1,
#'                                shape1 = .75,
#'                                shape2 = .25),
#'                 acsvalue = acs(id = "weibull",
#'                                t = 0:30,
#'                                scale = 10,
#'                                shape = .75),
#'                 n = 1000, p = 30, p0 = .5, TSn = 5)
#'
#' checkTS(x)
#'
checkTS <- function(TS, distbounds = c(-Inf, Inf)) {

  if (!is.list(TS)) TS <- list(TS)

  att      <- attributes(TS[[1]])
  margdist <- att$margdist
  margarg  <- att$margarg
  p0       <- att$p0
  acsvalue <- att$acsvalue

  ac <- sapply(TS, function(x) acf(x, plot = FALSE)$acf)

  out <- data.frame(
    mean = c(popmean(margdist, margarg, distbounds = distbounds, p0 = p0),
             sapply(TS, mean)),
    sd   = c(popsd(margdist, margarg, distbounds = distbounds, p0 = p0),
             sapply(TS, sd)),
    skew = c(popskew(margdist, margarg, distbounds = distbounds, p0 = p0),
             sapply(TS, function(x) sample.moments(x, raw = FALSE,
                                                   central = FALSE,
                                                   coef = TRUE)$coefficients[2])),
    p0     = c(p0, sapply(TS, function(x) length(which(x == 0)) / length(x))),
    acf_t1 = c(acsvalue[2], ac[2, ]),
    acf_t2 = c(acsvalue[3], ac[3, ]),
    acf_t3 = c(acsvalue[4], ac[4, ])
  )

  row.names(out) <- c("expected", paste0("simulation", seq_along(TS)))

  structure(.Data    = as.matrix(round(out, 2)),
            class    = c("checkTS", "matrix"),
            margdist = margdist,
            margarg  = margarg,
            p0       = p0)
}
