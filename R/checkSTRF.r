#' Numerical and visual check of generated random fields
#'
#' Compares generated random fields sample statistics with the theoretically expected values (similar to checkTS). It also returns graphical output for visual check
#'
#' @param STRF output of 'generateSTRF' or 'generateSTRFsepfast'
#' @param lags number of lags of empirical STCF to be considered in the graphical output (default set to 30)
#'
#' @name checkSTRF
#'
#' @import viridis  grDevices
#'
#' @importFrom MBA mba.surf
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
#' sim <- generateSTRF(n = 10,
#'                     STmodel = fit)
#' checkSTRF(STRF = sim,
#'           lags = 10,
#'           nfields = 10)
#'
#'
checkSTRF <- function(STRF, lags = 30, nfields = 50) {

    att <- attributes(STRF)$STmodel
    margdist <- att$margdist
    margarg <- att$margarg
    p0 <- att$p0
    stcs <- att$stcs
    s.lags <- att$s.lags
    t.lags <- att$t.lags
    grid.lags <- att$grid.lags
    distbounds <-  att$distbounds
    m <- sqrt(ncol(STRF))
    if (m^2 < 20){
        npoints <- m^2
    } else {
        npoints <- 20
    }
    stcs.sim <- acf(as.matrix(STRF)[ , 1:npoints], lags, plot=F)[[1]][ , , 1]

    # dev.new(w = 15, h = 11)
    op <- par(no.readonly = TRUE) # the whole list of settable par's.

    par(mar = c(4, 4, .5, 2), mgp = c(2.5, 1, 0), las = 1)

    layout(mat = matrix(data = 1:3,
                        nrow = 3,
                        byrow = TRUE))

    plot(grid.lags[1:npoints], stcs.sim[1, ], type = "n", xlab=expression(paste('Distance ', delta)),
         ylab="Correlation", ylim=c(0,1))
    for (i in 1:3) {
        points(grid.lags[1:npoints], stcs.sim[i, ], col = i, pch = 19)
        lines(s.lags, stcs[i, ], col = i, lwd = 2)
    }
    legend("topright", legend=c(expression(paste(tau," = 0")),
                                expression(paste(tau," = 1")),
                                expression(paste(tau," = 2")),
                                "Simulated"),
           col = c(1:3, 1), lwd = 2, pch = c(rep(NA, 3), 19), lty = c(rep(1, 3), NA), bg = "darkgrey", box.col = "darkgrey", text.col = "white")
    box()

    contour(x = s.lags, y = 1:lags - 1, z = t(stcs[1:lags, ]), col = viridis_pal(option = 'D')(10),
            zlim = c(0,1), xlim = c(0, tail(s.lags[stcs[1, ] > 0.05], 1)), ylim = c(0, tail((1:length(stcs[ ,1]))[stcs[ ,1] > 0.05],1)),
            labcex = 0.8, ylab = expression(paste('Time ', tau)), xlab = expression(paste('Distance ', delta)), lwd = 2)
    legend("topright", 'Target STCS', adj = c(0.2, 0.25), bg = "darkgrey", box.col = "darkgrey", text.col = "white")
    box()

    if(p0 != 0) xpos <- as.vector(STRF[STRF > 0])
    else xpos <- as.vector(STRF)
    plot(sort(xpos), ppoints(xpos,0), type = 'l', lwd = 2,  xlab = 'X', ylab = 'CDF', xlim = c(min(xpos), quantile(xpos, 0.99)))
    aux <- do.call(paste0("q", margdist), args = c(list(p = seq(0,0.999, by = 0.001)), margarg))
    lines(aux, do.call(paste0("p", margdist), args = c(list(q = aux), margarg)) ,col = 2, lwd = 2, lty = 2)
    legend("bottomright", legend=c("Simulated ECDF", "Target CDF"),
           col = c(1:2), lwd = 2, pch = c(rep(NA,2)), lty = c(1,2), bg = "darkgrey", box.col = "darkgrey", text.col = "white")
    box()


    # dev.new(w = 11 * 1.2, h = 13 * 1.2)
    par(mar = c(0, 0, 0, 0),
        oma = c(0.5, 0.5, 0.5, 0.5))

    layout(mat = matrix(data = 1:nfields,
                        nrow = n2mfrow(nfields)[1],
                        byrow = TRUE))

    if(p0 > 0 & sample.moments(STRF[1:nfields, ], raw = F, central = F, coef = T)$coefficients[2] > 1) {

        aux  <-  sqrt(STRF[1:nfields, ])
    }  else {

        aux <- STRF[1:nfields, ]
    }
    for(i in 1:nfields){
        image(matrix(aux[i, ], m, m), axes=F, col = my.colors(20), zlim = range(sqrt(STRF[1:nfields,])))
        legend("topleft", legend = i, adj = c(1.3, 0.1), cex = 1.2,
               pch = NA, lty = 0, bty = 'n', text.col = "darkgrey", text.font = 2)
        box(col = "darkgrey")
    }

    par(op)

    out <- data.frame(mean = c(popmean(margdist, margarg, distbounds = distbounds,
                      p0 = p0), apply(STRF[,1:npoints], 2, mean)), sd = c(popsd(margdist, margarg,
                      distbounds = distbounds, p0 = p0), apply(STRF[,1:npoints], 2, sd)), skew = c(popskew(margdist,
                      margarg, distbounds = distbounds, p0 = p0), apply(STRF[,1:npoints], 2,
                      function(x) sample.moments(x, raw = F, central = F, coef = T)$coefficients[2])),
                      p0 = c(p0, apply(STRF[,1:npoints], 2, function(x) length(which(x == 0))/length(x))))
                      row.names(out) <- c("expected", paste0("sample location ", 1:npoints))
    out <- round(out, 2)

    structure(.Data = as.matrix(out), class = c("matrix"),
              margdist = margdist, margarg = margarg, p0 = p0)
}


