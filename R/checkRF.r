#' Numerical and visual check of generated random fields
#'
#' Compares generated random fields sample statistics with the theoretically
#' expected values (similar to checkTS). It also returns graphical output for
#' visual check
#'
#' @param RF output of \code{\link{generateRF}}
#' @param lags number of lags of empirical STCF to be considered in the
#' graphical output (default set to 30)
#' @param nfields number of fields to be used in the numerical and graphical
#' output (default set to 49). As the plots are arranged in a matrix with nrows as close as possible to ncol, we suggest using values such as 3x3, 3x4, 7x8, etc.
#' @param method report method - 'stat' for basic statistical report,
#' 'statplot'for graphical check of lagged SCS, target STCS, and marginal
#' distribution, 'field' for plotting a matrix of the first 'nfields',
#' and 'movie' to save the first 'nfields' as a GIF file named "movieRF.gif"
#' in the current working directory

#' @name checkRF
#'
#' @import  grDevices  directlabels  cowplot  animation
#'
#' @importFrom MBA mba.surf
#'
#' @export
#'
#' @examples
#' ## The example below refers to the fitting and simulation of 10 random fields
#' ## of size 10x10 with AR(1) temporal correlation. As the fitting algorithm has
#' ## O((mxm)^3) complexity for a mxm field, this setting allows for quick fitting
#' ## and simulation (short CPU time). However, for a more effective visualization
#' ## and reliable performance assessment, we suggest to generate a larger number
#' ## of fields (e.g. 100 or more) of size about 30X30. This setting needs more
#' ## CPU time but enables more effective comparison of theoretical and
#' ## empirical statistics. Sizes larger than about 50x50 can be unpractical
#' ##  on standard machines.
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
#' checkRF(RF = sim,
#'           lags = 10,
#'           nfields = 12)
#'
#'
checkRF <- function(RF, lags = 30, nfields = 49, method = 'stat') {

  att <- attributes(RF)$STmodel
  margdist <- att$margdist
  margarg <- att$margarg
  p0 <- att$p0
  stcs <- att$stcs
  s.lags <- att$s.lags
  t.lags <- att$t.lags
  grid.lags <- att$grid.lags
  distbounds <-  att$distbounds
  m <- sqrt(ncol(RF))
  if (m^2 < 20){
    npoints <- m^2
  } else {
    npoints <- 20
  }
  stcs.sim <- acf(as.matrix(RF)[ , 1:npoints], lag.max = lags, plot=F)[[1]][ , , 1]


  if(method == 'stat') {
    out <- data.frame(mean = c(popmean(margdist, margarg, distbounds = distbounds,
                                       p0 = p0), apply(RF[,1:npoints], 2, mean)), sd = c(popsd(margdist, margarg,
                                                                                               distbounds = distbounds, p0 = p0), apply(RF[,1:npoints], 2, sd)), skew = c(popskew(margdist,
                                                                                                                                                                                  margarg, distbounds = distbounds, p0 = p0), apply(RF[,1:npoints], 2,
                                                                                                                                                                                                                                    function(x) sample.moments(x, raw = F, central = F, coef = T)$coefficients[2])),
                      p0 = c(p0, apply(RF[,1:npoints], 2, function(x) length(which(x == 0))/length(x))))
    row.names(out) <- c("expected", paste0("sample location ", 1:npoints))
    out <- round(out, 2)

    structure(.Data = as.matrix(out), class = c("matrix"),
              margdist = margdist, margarg = margarg, p0 = p0)
    return(out)

  }


  if(method == 'statplot') {

    ## STCS  slices

    DataS <- data.table(grid.lags[1:npoints], t(as.matrix(stcs.sim[1:3, ])))
    colnames(DataS) <- c('Lag', 'scf1', 'scf2', 'scf3')
    DataS <- melt(DataS, id.vars = 'Lag')

    i.lag <- which(s.lags <= max(grid.lags[1:npoints]) )
    DataT <- data.table(s.lags[i.lag], t(as.matrix(stcs[1:3, i.lag])))
    colnames(DataT) <- c('Lag', 'scf1', 'scf2', 'scf3')
    DataT <- melt(DataT, id.vars = 'Lag')

    ran <- range(c(DataT$Lag, DataS$Lag))

    Lag <- value <- variable <- NULL

    p1 <- ggplot() +
      geom_point(data = DataS,
                 aes(x = Lag,
                     y = value,
                     colour = (variable)),
                 alpha = .5) +
      geom_line(data = DataT,
                aes(x = Lag,
                    y = value,
                    colour = (variable)),
                lwd = .5,
                alpha = 1) +
      xlim(ran) +
      scale_color_viridis_d(labels = c(expression(paste(tau, ' = 0')),
                                      expression(paste(tau, ' = 1')),
                                      expression(paste(tau, ' = 2'))
                                      ),
                           option = 'D') +
      labs(x = expression(paste('Distance ', delta)),
           y = 'SCS',
           title = 'Lagged spatial correlation structure',
           color = "Time lag") +
      theme_light() +
      theme(legend.position=c(0.99,0.99),
            legend.justification=c(1,1),
            strip.background = element_rect(fill = 'grey5'),
            strip.text = element_text(colour = 'grey95'))


    ## STCS plot
    Data0 <- cbind(expand.grid(s.lags, (1:lags) - 1), as.numeric(t(stcs[1:lags, ])) )
    colnames(Data0) <- c('x', 'y', 'z')

    p2.0 <- ggplot() +
      stat_contour(data = Data0, aes_string(x = "x", y = "y",
                                    z = "z", colour = "..level.."),
                   breaks = seq(0, 1, 0.1),
                   size = 1) +
      scale_color_viridis_c(option = 'D') +
      labs(x = expression(paste('Distance ', delta)),
           y = expression(paste('Time lag ', tau)),
           title = 'Target STCS') +
      theme_light() +
      theme(legend.position= 'none',
            legend.justification=c(1,1),
            strip.background = element_rect(fill = 'grey5'),
            strip.text = element_text(colour = 'grey95'))

    p2 <- directlabels::direct.label(p2.0, list("bottom.pieces", cex = 0.8)) #


    ## CDF plot
    if(p0 != 0) xpos <- as.vector(RF[RF > 0])
    else xpos <- as.vector(RF)
    if(length(xpos) > 1000) xpos <- xpos[1:1000]

    Data.simu <- data.table(x = sort(xpos), cdf = ppoints(xpos,0), type = "simu")

    Data.theo <- data.table(
      x = do.call(paste0("q", margdist),
                  args = c(list(p = seq(0.001, 0.999, by = 0.001)), margarg)),
      cdf = seq(0.001, 0.999, by = 0.001),
      type = 'theo')

    Data <- rbind(Data.simu, Data.theo)

    p3 <- ggplot(data = Data) +
      geom_line(aes(x = Data$x,
                    y = Data$cdf,
                    colour = factor(Data$type),
                    linetype  = factor(Data$type)),
                lwd = 1) +
      xlim(range(Data$x)) +
      ## scale_color_viridis_d(labels = c('Simulated', 'Fitted')) +
      scale_color_manual(labels = c('Simulated', 'Fitted'),
                            values = c('royalblue4', 'red4')) +
      scale_linetype(labels = c('Simulated', 'Fitted')) +
      labs(x = 'Nonzero values',
           y = 'Nonexceedence probability',
           title = 'Probability distribution fit') +
      guides(color = guide_legend(title=NULL),
             linetype = guide_legend(title=NULL)) +
      theme_light() +
      theme(legend.position=c(0.99,0.01),
            legend.justification=c(1,0),
            strip.background = element_rect(fill = 'grey5'),
            strip.text = element_text(colour = 'grey95'))

    out <- ggdraw() +
      draw_plot(p1, x = 0, y = .5, width = 1, height = .5) +
      draw_plot(p2, x = 0, y = 0, width = .5, height = .5) +
      draw_plot(p3, x = .5, y = 0, width = .5, height = .5)

    return(out)
  }

  my.colors <- colorRampPalette( c(
    "#ffffff","#610B11", "#80110A", "#AF2F03", "#D05201", "#E87D0A", "#EE9E1B", "#F3C03B",
    "#F7DE68", "#FAF397", "#E8D794", "#DABD77", "#CBA265", "#BC8854", "#A96E48",
    "#965F4A", "#7C4F47", "#322221", "#3D2E2C", "#4F413B", "#60534B", "#72675B",
    "#877E6F", "#A39E8A", "#C0BDA8", "#E2E2CD", "#E0CD99", "#C0AE6F", "#A09253",
    "#81783C", "#645F27", "#545316", "#43440A", "#303401", "#173F3B", "#1D6447",
    "#1E803A", "#169517", "#42B218", "#8DCC3C", "#C8E065"))

  if(method == 'field') {

    ## Matrix of fields
    op <- par(no.readonly = TRUE) # the whole list of settable par's.
    par(mar = c(0, 0, 0, 0),
        oma = c(0.5, 0.5, 0.5, 0.5))
    graphics::layout(mat = matrix(data = 1:nfields,
                                  nrow = n2mfrow(nfields)[1],
                                  byrow = TRUE))

    if(p0 > 0 & sample.moments(RF[1:nfields, ], raw = F, central = F, coef = T)$coefficients[2] > 1) {

      aux  <-  sqrt(RF[1:nfields, ])
    }  else {

      aux <- RF[1:nfields, ]
    }
    d <- u <- seq(0, m - 1, length = m)
    ud <- expand.grid(u, d)
    for(i in 1:nfields){
      z <- mba.surf(data.frame(x = ud[ ,1], y = ud[ ,2], z = aux[i, ]), 100, 100)$xyz.est$z
      if(p0 > 0) z[z < 0] <- 0
      image(matrix(z, 100, 100), axes=F, col = my.colors(40), zlim = range(aux))
      legend("topleft", legend = i, adj = c(1.3, 0.1), cex = 1.2,
             pch = NA, lty = 0, bty = 'n', text.col = "darkgrey", text.font = 2)
      box(col = "darkgrey")
    }

    par(op)

  }


  if(method == 'movie') {

    d <- u <- seq(0, m - 1, length = m)
    ud <- expand.grid(u, d)
    op <- par(mar = c(0, 0, 0, 0), oma = c(0.5, 0.5, 0.5, 0.5))
    if(p0 > 0 & sample.moments(RF[1:nfields, ], raw = F, central = F, coef = T)$coefficients[2] > 1) aux  <-  sqrt(RF[1:nfields, ])
    else aux <- RF[1:nfields, ]
    animation::ani.options(loop=0)
    animation::saveGIF({
      for(i in 1:nfields){
        z <- mba.surf(data.frame(x = ud[ ,1], y = ud[ ,2], z = aux[i, ]), 100, 100)$xyz.est$z
        if(p0 > 0) z[z < 0] <- 0
        image(matrix(z, 100, 100), axes=F, col = my.colors(40), zlim = range(sqrt(RF[1:nfields,])))
        legend("topleft", legend = i, adj = c(1.3, 0.1), cex = 1.2,
               pch = NA, lty = 0, bty = 'n', text.col = "darkgrey", text.font = 2)
        box(col = "darkgrey")
      }
    }, movie.name = "movieRF.gif",
    interval = 0.5,  ani.width = 480, ani.height = 480)
    par(op)

  }


}


