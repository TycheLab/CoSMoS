#' SpatioTemporal Correlation Structure
#'
#' Provides a parametric function that describes the values of the linear spatiotemporal autocorrelation up to desired lags. For more details on the parametric spatiotemporal correlation structures see section 2.3 and 2.4 in \href{https://doi.org/10.1029/2019WR026331}{Papalexiou and Serinaldi (2020)}
#'
#' @param id spatiotemporal correlation structure ID
#' @param ... additional arguments (t as time lag, s as spatial lag (distance),  and stcs parameters)
#'
#' @name stcs
#'
#' @export
#'
#' @examples
#'
#' library(CoSMoS)
#' library(plotly)
#'
#' ## specify grid of spatial and temporal lags
#' d <- 31
#' st <- expand.grid(0:(d-1),
#'                   0:(d-1))
#'
#' ## get the STCS
#' wc <- stcs("clayton",
#'            t = st[, 1],
#'            s = st[, 2],
#'            scfid = 'weibull',
#'            tcfid = 'weibull',
#'            copulaarg = 2,
#'            scfarg = list(scale = 20,
#'                          shape = 0.7),
#'            tcfarg = list(scale = 1.1,
#'                          shape = 0.8))
#'
#' g14 <- stcs("gneiting14",
#'             t = st[, 1],
#'             s = st[, 2],
#'             a = 1/50,
#'             c = 1/10,
#'             alpha = 1,
#'             beta = 1,
#'             gamma = 0.5,
#'             tau = 1)
#'
#' g16 <- stcs("gneiting16",
#'             t = st[, 1],
#'             s = st[, 2],
#'             a = 1/50,
#'             c = 1/10,
#'             alpha = 1,
#'             beta = 1,
#'             nu = 0.5,
#'             tau = 1)
#'
#' ## note: for nu = 0.5 stcfgneiting16 is equivalent to
#' ## stcfgneiting14 with gamma = 0.5
#'
#' ## visualize the STCS
#'
#' wc.m <- matrix(wc,
#'                nrow = d)
#'
#' plot_ly(z = ~wc.m) %>%
#'     add_surface() %>%
#'     layout(
#'         scene = list(
#'             xaxis = list(title = "Time lag"),
#'             yaxis = list(title = "Distance"),
#'             zaxis = list(title = "STCF")
#'         )
#'     ) %>%
#'     hide_colorbar()
#'
#' g14.m <- matrix(g14,
#'                 nrow = d)
#'
#' plot_ly(z = ~g14.m) %>%
#'     add_surface() %>%
#'     layout(
#'         scene = list(
#'             xaxis = list(title = "Time lag"),
#'             yaxis = list(title = "Distance"),
#'             zaxis = list(title = "STCF")
#'         )
#'     ) %>%
#'     hide_colorbar()
#'
stcs <- function (id, ...) {

    .args <- list(...)
    do.call(paste0("stcf", id), args = .args)
}

#' SpatioTemporal Correlation Structure
#'
#' Provides a parametric function that describes the values of the linear spatiotemporal autocorrelation up to desired lags. For more details on the parametric spatiotemporal correlation structures see section 2.3 and 2.4 in \href{https://doi.org/10.1029/2019WR026331}{Papalexiou and Serinaldi (2020)}
#'
#' @param id spatiotemporal correlation structure ID
#' @param arglist list of additional arguments (t as time lag, s as spatial lag (distance),  and stcs parameters)
#'
#' @keywords internal
#'
stcs2 <- function (id, arglist) {

    do.call(paste0("stcf", id), args = arglist)
}



#' Clayton SpatioTemporal Correlation Structure
#'
#' Provides spatiotemporal correlation structure function based on Clayton copula. For more details on the parametric spatiotemporal correlation structures see section 2.3 and 2.4 in \href{https://doi.org/10.1029/2019WR026331}{Papalexiou and Serinaldi (2020)}
#'
#' @param t time lag
#' @param s spatial lag (distance)
#' @param scfid ID of the spatial (marginal) correlation structure (e.g. weibull)
#' @param tcfid ID of the temporal (marginal) correlation structure (e.g. weibull)
#' @param copulaarg parameter of the Clayton copula linking the marginal correlation structures
#' @param scfarg parameters of spatial (marginal) correlation structure
#' @param tcfarg parameters of temporal (marginal) correlation structure
#'
#' @name stcfclayton
#'
#' @export
#'
#' @examples
#'
#' library(CoSMoS)
#' library(plotly)
#'
#' ## specify grid of spatial and temporal lags
#' d <- 31
#' st <- expand.grid(0:(d - 1),
#'                   0:(d - 1))
#'
#' ## get the STCS
#' wc <- stcfclayton(t = st[, 1],
#'                   s = st[, 2],
#'                   scfid = 'weibull',
#'                   tcfid = 'weibull',
#'                   copulaarg = 2,
#'                   scfarg = list(scale = 20,
#'                                 shape = 0.7),
#'                   tcfarg = list(scale = 1.1,
#'                                 shape = 0.8))
#'
#' ## visualize the STCS
#' wc.m <- matrix(wc,
#'                nrow = d)
#'
#' plot_ly(z = ~wc.m) %>%
#'     add_surface() %>%
#'     layout(
#'         scene = list(
#'             xaxis = list(title = "Time lag"),
#'             yaxis = list(title = "Distance"),
#'             zaxis = list(title = "STCF")
#'         )
#'     ) %>%
#'     hide_colorbar()
#'
stcfclayton <- function (t, s, scfid, tcfid, copulaarg, scfarg, tcfarg) {

    if (copulaarg < 0) {
        return(NaN)
    }
    else {
        scfarg$t <- s
        tcfarg$t <- t
        scf <- do.call(paste0("acf", scfid), args = scfarg)
        tcf <- do.call(paste0("acf", tcfid), args = tcfarg)
        if(copulaarg > 0) stcf <- (tcf^(-copulaarg) + scf^(-copulaarg) - 1 )^(-1/copulaarg)
        if(copulaarg == 0) stcf <- tcf * scf
        return(stcf)
    }
}

#' Gneiting-14 SpatioTemporal Correlation Structure
#'
#' Provides spatiotemporal correlation structure function proposed by \href{https://doi.org/10.1198/016214502760047113}{Gneiting (2002)} (Eq.14 at p. 593)
#' @param t time lag
#' @param s spatial lag (distance)
#' @param a nonnegative scaling parameter of time
#' @param c nonnegative scaling parameter of space
#' @param alpha smoothness parameter of time. Valid range: \eqn{(0,1]}
#' @param gamma smoothness parameter of space. Valid range: \eqn{(0,1]}
#' @param beta space-time interaction parameter. Valid range: \eqn{[0,1]}
#' @param tau space-time interaction parameter. Valid range: \eqn{\ge 1} (for 2-dimensional fields)
#'
#' @name stcfgneiting14
#'
#' @export
#'
#' @examples
#'
#' library(CoSMoS)
#' library(plotly)
#'
#' ## specify grid of spatial and temporal lags
#' d <- 31
#' st <- expand.grid(0:(d - 1),
#'                   0:(d - 1))
#'
#' ## get the STCS
#' g14 <- stcfgneiting14(t = st[, 1],
#'                       s = st[, 2],
#'                       a = 1/50,
#'                       c = 1/10,
#'                       alpha = 1,
#'                       beta = 1,
#'                       gamma = 0.5,
#'                       tau = 1)
#'
#' ## visualize the STCS
#'
#' g14.m <- matrix(g14,
#'                 nrow = d)
#'
#' plot_ly(z = ~g14.m) %>%
#'     add_surface() %>%
#'     layout(
#'         scene = list(
#'             xaxis = list(title = "Time lag"),
#'             yaxis = list(title = "Distance"),
#'             zaxis = list(title = "STCF")
#'         )
#'     ) %>%
#'     hide_colorbar()
#'
stcfgneiting14 <- function (t, s, a, c, alpha, beta, gamma, tau) {

    if((a < 0) | (c < 0) | (alpha <= 0) | (alpha > 1) |
       (beta < 0) | (beta > 1) | (gamma <= 0) | (gamma > 1) |
       (tau < 1)) {

        return(NaN)
    } else {

        stcf <- 1 / ((a * abs(t)^(2 * alpha) + 1)^tau) *
            exp( - (c * abs(s)^(2 * gamma)) / ((a * abs(t)^(2 * alpha) + 1)^(beta * gamma)) )

        return(stcf)
    }
}

#' Gneiting-16 SpatioTemporal Correlation Structure
#'
#' Provides spatiotemporal correlation structure function proposed by \href{https://doi.org/10.1198/016214502760047113}{Gneiting (2002)} (Eq.16 at p. 594)
#'
#' @param t time lag
#' @param s spatial lag (distance)
#' @param a nonnegative scaling parameter of time
#' @param c nonnegative scaling parameter of space
#' @param alpha smoothness parameter of time. Valid range: \eqn{(0,1]}
#' @param nu smoothness parameter of space. Valid range: \eqn{>0}
#' @param beta space-time interaction parameter. Valid range: \eqn{[0,1]}
#' @param tau space-time interaction parameter. Valid range: \eqn{\ge 1} (for 2-dimensional fields)
#'
#' @name stcfgneiting16
#'
#' @export
#'
#' @examples
#'
#' library(CoSMoS)
#' library(plotly)
#'
#' ## specify grid of spatial and temporal lags
#' d <- 31
#' st <- expand.grid(0:(d - 1),
#'                   0:(d - 1))
#'
#' ## get the STCS
#' g16 <- stcfgneiting16(t = st[, 1],
#'                       s = st[, 2],
#'                       a = 1/50,
#'                       c = 1/10,
#'                       alpha = 1,
#'                       beta = 1,
#'                       nu = 0.5, tau = 1)
#'
#' ## visualize the STCS
#'
#' g16.m <- matrix(g16,
#'                 nrow = d)
#'
#' plot_ly(z = ~g16.m) %>%
#'     add_surface() %>%
#'     layout(
#'         scene = list(
#'             xaxis = list(title = "Time lag"),
#'             yaxis = list(title = "Distance"),
#'             zaxis = list(title = "STCF")
#'         )
#'     ) %>%
#'     hide_colorbar()
#'
stcfgneiting16 <- function (t, s, a, c, alpha, beta, nu, tau) {

    if((a < 0) | (c < 0) | (alpha <= 0) | (alpha > 1) |
       (beta < 0) | (beta > 1) | (nu <= 0) | (tau < 1)) {
        return(NaN)
    }
    else {
        stcf <- 1 / (2^(nu - 1)  * gamma(nu) * (a * abs(t)^(2 * alpha) + 1)^tau) *
            ((c * abs(s)) / ((a * abs(t)^(2 * alpha) + 1)^(beta / 2) ))^nu *
            besselK((c * abs(s)) / ((a * abs(t)^(2 * alpha) + 1)^(beta / 2) ), nu)
        return(stcf)
    }
}

