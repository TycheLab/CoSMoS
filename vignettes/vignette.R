## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(eval = FALSE, ###
                      echo = TRUE,
                      fig.width = 7, 
                      warning = FALSE,
                      message = FALSE)

## ---- message = FALSE---------------------------------------------------------
#  library(CoSMoS)

## -----------------------------------------------------------------------------
#  marginaldist <- "ggamma"
#  param <- list(scale = 1,
#                shape1 = .8,
#                shape2 = .8)

## -----------------------------------------------------------------------------
#  no <- 1000
#  ggamma_sim <- rggamma(n = no,
#                        scale = 1,
#                        shape1 = 1,
#                        shape2 = .5)
#  plot.ts(x = ggamma_sim,
#          main = "")
#  
#  par(mfrow = c(1, 2))
#  
#  plot(x = density(x = ggamma_sim),
#       main = "")
#  acf(x = ggamma_sim,
#      main = "")

## -----------------------------------------------------------------------------
#  acf <- c(1, 0.8)
#  ggamma_sim <- generateTS(n = no,
#                           margdist = marginaldist,
#                           margarg = param,
#                           acsvalue = acf)
#  plot(x = ggamma_sim,
#       main = "")
#  
#  par(mfrow = c(1, 2))
#  
#  plot(x = density(ggamma_sim[[1]]),
#       main = "")
#  acf(x = ggamma_sim[[1]],
#      main = "")

## -----------------------------------------------------------------------------
#  acf <- c(1, 0.5, 0.5, 0.4, 0.4) #up to lag-4
#  ggamma_sim <- generateTS(n = no,
#                           margdist = marginaldist,
#                           margarg = param,
#                           acsvalue = acf)
#  plot(x = ggamma_sim,
#       main = "")
#  
#  par(mfrow = c(1, 2))
#  
#  plot(density(x = ggamma_sim[[1]]),
#       main = "")
#  acf(x = ggamma_sim[[1]],
#      main = "")

## -----------------------------------------------------------------------------
#  ## specify lag
#  lags <- 0:10
#  
#  ## get the ACS
#  f <- acs(id = "fgn",
#           t = lags,
#           H = .75)
#  b <- acs(id = "burrXII",
#           t = lags,
#           scale = 1,
#           shape1 = .6,
#           shape2 = .4)
#  w <- acs(id = "weibull",
#           t = lags,
#           scale = 2,
#           shape = 0.8)
#  p <- acs(id = "paretoII",
#           t = lags,
#           scale = 3,
#           shape = 0.3)
#  
#  ## visualize the ACS
#  dta <- data.table(lags, f, b, w, p)
#  
#  m.dta <- melt(data = dta,
#                id.vars = "lags")
#  
#  ggplot(data = m.dta,
#         mapping = aes(x = lags,
#                       y = value,
#                       group = variable,
#                       colour = variable)) +
#    geom_point(size = 2.5) +
#    geom_line(lwd = 1) +
#    scale_color_manual(values = c("steelblue4", "red4", "green4", "darkorange"),
#                       labels = c("FGN", "Burr XII", "Weibull", "Pareto II"),
#                       name = "") +
#    labs(x = bquote(lag ~ tau),
#         y = "Acf") +
#    scale_x_continuous(breaks = lags) +
#    theme_grey()

## -----------------------------------------------------------------------------
#  acf = acs(id = "paretoII",
#            t = 0:30,
#            scale = 1,
#            shape = .75)
#  
#  ggamma_sim <- generateTS(n = no,
#                           margdist = marginaldist,
#                           margarg = param,
#                           acsvalue = acf)
#  
#  plot(x = ggamma_sim,
#       main = "")

## -----------------------------------------------------------------------------
#  my_acf <- exp(seq(0, -2, -0.1))
#  ggamma_sim <- generateTS(n = no,
#                           margdist = marginaldist,
#                           margarg = param,
#                           acsvalue = my_acf)
#  
#  plot(x = ggamma_sim,
#       main = "")
#  
#  acf(x = ggamma_sim[[1]])

## ---- fig.height = 7----------------------------------------------------------
#  prob_zero <- .9
#  ggamma_sim <- generateTS(n = no,
#                           margdist = marginaldist,
#                           margarg = param,
#                           acsvalue = acf,
#                           p0 = prob_zero,
#                           TSn = 5)
#  
#  par(mfrow = c(1, 1))
#  
#  plot(x = ggamma_sim,
#       main = "")

## -----------------------------------------------------------------------------
#  checkTS(ggamma_sim)

## ---- fig.height = 7----------------------------------------------------------
#  library(plotly)
#  
#  ## specify grid of spatial and temporal lags
#  d <- 31
#  st <- expand.grid(0:(d - 1),
#                    0:(d - 1))
#  
#  ## get the STCS
#  wc <- stcfclayton(t = st[, 1],
#                    s = st[, 2],
#                    scfid = "weibull",
#                    tcfid = "weibull",
#                    copulaarg = 2,
#                    scfarg = list(scale = 20,
#                                  shape = 0.7),
#                    tcfarg = list(scale = 1.1,
#                                  shape = 0.8))
#  
#  ## visualize the STCS
#  wc.m <- matrix(data = wc,
#                 nrow = d)
#  
#  plot_ly(z = ~wc.m) %>%
#    add_surface() %>%
#    layout(
#        scene = list(
#            xaxis = list(title = "Time lag"),
#            yaxis = list(title = "Distance"),
#            zaxis = list(title = "STCF")
#        )
#    ) %>%
#    hide_colorbar()

## ---- fig.height = 7----------------------------------------------------------
#  ## set a sequence of hypothetical coordinates
#  coord <- cbind(runif(8)*30,
#                 runif(8)*30)
#  
#  ## compute VAR model parameters
#  fit <- fitVARMULTI(
#    coord = coord,
#    p = 5,
#    margdist ="burrXII",
#    margarg = list(scale = 3,
#                   shape1 = .9,
#                   shape2 = .2),
#    p0 = 0.8,
#    stcsid = "clayton",
#    stcsarg = list(scfid = "weibull",
#                   tcfid = "weibull",
#                   copulaarg = 2,
#                   scfarg = list(scale = 25,
#                                 shape = 0.7),
#                   tcfarg = list(scale = 3.1,
#                                 shape = 0.8))
#  )
#  
#  ## generate correlated random vectors
#  sim <- generateMULTI(n = 500,
#                       STmodel = fit)
#  
#  ## visualize simulated timeseries
#  dta <- melt(data = data.table(time = 1:nrow(sim),
#                                sim[,1:8]),
#              id.vars = "time")
#  
#  ggplot(data = dta,
#         mapping = aes(x = time,
#                       y = value)) +
#         geom_line() +
#         facet_grid(facets = variable ~ .,
#                    scales = "free_y")

## -----------------------------------------------------------------------------
#  ## compute VAR model parameters
#  fit <- fitVARSTRF(
#    m = 20,
#    p = 10,
#    margdist ="burrXII",
#    margarg = list(scale = 3,
#                   shape1 = .9,
#                   shape2 = .2),
#    p0 = 0.8,
#    stcsid = "clayton",
#    stcsarg = list(scfid = "weibull",
#                   tcfid = "weibull",
#                   copulaarg = 2,
#                   scfarg = list(scale = 20,
#                                 shape = 0.7),
#                   tcfarg = list(scale = 1.1,
#                                 shape = 0.8))
#  )
#  ## generate random fields with nonseparable correlation structure
#  sim1 <- generateSTRF(n = 1000,
#                       STmodel = fit)
#  
#  ## fast simulation of random fields with separable correlation structure
#  sim2 <- generateSTRFsepfast(
#      n = 1000,
#      m = 20,
#      p0 = 0.7,
#      margdist ="paretoII",
#      margarg = list(scale = 1,
#                     shape = .3),
#      stcsarg = list(scfid = "weibull",
#                     tcfid = "weibull",
#                     scfarg = list(scale = 20,
#                                   shape = 0.7),
#                     tcfarg = list(scale = 1.1,
#                                   shape = 0.8))
#  )
#  
#  checkSTRF(STRF = sim2)

## -----------------------------------------------------------------------------
#  ## check random fields
#  checkSTRF(sim1)
#  checkSTRF(sim2)

## -----------------------------------------------------------------------------
#  data("precip")
#  plot(x = precip,
#       type = "l")
#  
#  par(mfrow = c(1, 2))
#  
#  plot(x = density(x = precip$value),
#       main = "",
#       xlim = c(0, 10)) #Does not plot extreme values
#  acf(x = precip$value,
#      main = "", lag.max = 20)

## ---- fig.height = 9----------------------------------------------------------
#  precip_ggamma <- analyzeTS(TS = precip,
#                                  season = "month",
#                                  dist = "ggamma",
#                                  acsID = "weibull",
#                                  lag.max = 12)
#  
#  reportTS(aTS = precip_ggamma,
#           method = "dist")
#  reportTS(aTS = precip_ggamma,
#           method = "acs")
#  reportTS(aTS = precip_ggamma,
#           method = "stat")

## ---- warning = FALSE, fig.height = 9-----------------------------------------
#  precip_pareto <- analyzeTS(TS = precip,
#                             season = "month",
#                             dist = "paretoII",
#                             acsID = "fgn",
#                             lag.max = 12)
#  
#  reportTS(aTS = precip_pareto,
#           method = "dist")
#  reportTS(aTS = precip_pareto,
#           method = "acs")

## -----------------------------------------------------------------------------
#  sim_precip <- simulateTS(aTS = precip_ggamma,
#                           from = as.POSIXct(x = "1978-12-01 00:00:00"),
#                           to = as.POSIXct(x = "2008-12-01 00:00:00"))
#  dta <- precip
#  dta[, id := "observed"]
#  sim_precip[, id := "simulated"]
#  
#  dta <- rbind(dta, sim_precip)
#  
#  ggplot(data = dta) +
#    geom_line(mapping = aes(x = date,
#                            y = value)) +
#    facet_wrap(facets = ~id,
#               ncol = 1) +
#    theme_classic()

## -----------------------------------------------------------------------------
#  # id <- "06606600"
#  # dta_raw <- as.data.table(read.fwf(file = sprintf(
#  #   "https://hydrology.nws.noaa.gov/pub/gcip/mopex/US_Data/Us_438_Daily/%s.dly", id),
#  #                                   widths = c(8,10,10,10,10,10),
#  #                                   col.names = c("date", "P", "E", "value", "Tmax", "Tmin")
#  #   )
#  # )
#  #
#  # dta_raw[, date := as.POSIXct(gsub(pattern = " ",
#  #                                   replacement = "0",
#  #                                   x = date),
#  #                              format = "%Y%m%d")]
#  #
#  # daily_streamflow <- dta_raw[value >= 0, .(date, value)]
#  # daily_streamflow <- daily_streamflow[1:1000,]
#  #
#  # str <- analyzeTS(TS = daily_streamflow,
#  #                  dist = "paretoII",
#  #                  norm = "N2",
#  #                  acsID = "paretoII",
#  #                  lag.max = 20,
#  #                  constrain = TRUE)
#  #
#  # reportTS(aTS = str)
#  #
#  # sim_str <- simulateTS(aTS = str)
#  #
#  # dta <- daily_streamflow
#  # dta[, id := "observed"]
#  # sim_str[, id := "simulated"]
#  #
#  # dta <- rbind(dta, sim_str)
#  #
#  # ggplot(data = dta) +
#  #   geom_line(mapping = aes(x = date,
#  #                           y = value)) +
#  #   facet_wrap(facets = ~id,
#  #              ncol = 1) +
#  #   theme_classic()

