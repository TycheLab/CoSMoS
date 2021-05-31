## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(eval = TRUE,
                      echo = TRUE,
                      fig.width = 7, 
                      warning = FALSE,
                      message = FALSE)
library(CoSMoS)
library(plot3D)

## ---- warning=FALSE, message=FALSE--------------------------------------------
## (i) specifying the sample size
no <- 1000
## (ii) defining the type of marginal distribution and its parameters
marginaldist <- "ggamma"
param <- list(scale = 1, shape1 = .8, shape2 = .8)
## (iii) defining the desired autocorrelation
acf.my <- c(1, 0.8)
## (iv) simulating
ggamma_sim <- generateTS(n = no, margdist = marginaldist, margarg = param, acsvalue = acf.my)
## and (v) visually checking the generated time series
quickTSPlot(ggamma_sim[[1]]) 


## ---- warning=FALSE, message=FALSE--------------------------------------------
acf <- c(1, 0.6, 0.5, 0.4, 0.3) #up to lag-4
ggamma_sim <- generateTS(n = no, margdist = marginaldist, margarg = param, acsvalue = acf)
quickTSPlot(ggamma_sim[[1]])

## ---- warning=FALSE, message=FALSE--------------------------------------------
## specify lag
lags <- 0:10

## get the ACS
f <- acs(id = "fgn", t = lags, H = .75)
b <- acs(id = "burrXII", t = lags, scale = 1, shape1 = .6, shape2 = .4)
w <- acs(id = "weibull", t = lags, scale = 2, shape = 0.8)
p <- acs(id = "paretoII", t = lags, scale = 3, shape = 0.3)

## visualize the ACS
dta <- data.table(lags, f, b, w, p)
m.dta <- melt(data = dta, id.vars = "lags")

ggplot(data = m.dta, mapping = aes(x = lags, y = value, group = variable, colour = variable)) + 
  geom_point(size = 2.5) + geom_line(lwd = 1) +
  scale_color_manual(values = c("steelblue4", "red4", "green4", "darkorange"), 
  labels = c("FGN", "Burr XII", "Weibull", "Pareto II"), name = "") +
  labs(x = bquote(lag ~ tau), y = "ACS") + scale_x_continuous(breaks = lags) + theme_light()

## ---- warning=FALSE, message=FALSE--------------------------------------------
acf <- acs(id = "paretoII", t = 0:30, scale = 1, shape = .75)
ggamma_sim <- generateTS(n = no, margdist = marginaldist, margarg = param, acsvalue = acf)
dta <- data.frame(time = 1:no, value = ggamma_sim[[1]])

quickTSPlot(dta$value)

## ---- warning=FALSE, message=FALSE--------------------------------------------
my_acf <- exp(seq(0, -2, -0.1))
ggamma_sim <- generateTS(n = no, margdist = marginaldist, margarg = param, acsvalue = my_acf)
quickTSPlot(ggamma_sim[[1]])


## ---- fig.height = 7, warning=FALSE, message=FALSE----------------------------
prob_zero <- .9
## the argument `TSn = 5` enables the simulation of 5 timeseries.
ggamma_sim <- generateTS(n = no, margdist = marginaldist, margarg = param, acsvalue = acf, 
                         p0 = prob_zero, TSn = 5)
plot(x = ggamma_sim, main = "") + theme_light()

## ---- warning=FALSE, message=FALSE--------------------------------------------
checkTS(ggamma_sim)

## ---- fig.height = 7, warning=FALSE, message=FALSE----------------------------

## specify grid of spatial and temporal lags
d <- 51
st <- expand.grid(0:(d - 1), 0:(d - 1))

## get the STCS
wc <- stcfclayton(t = st[, 1], s = st[, 2], scfid = "weibull", tcfid = "weibull", copulaarg = 2,
                  scfarg = list(scale = 20, shape = 0.7), tcfarg = list(scale = 5.1,shape = 0.8))

## visualize the STCS
wc.m <- matrix(data = wc, nrow = d)
j <- tail(which(wc.m[1, ] > 0.15), 1)
i <- tail(which(wc.m[, 1] > 0.15), 1)
wc.m <- wc.m[1:i, 1:j]

persp3D(z = wc.m, x = 1: nrow(wc.m), y = 1:ncol(wc.m),
expand = 1, main = "", scale = TRUE, facets = TRUE,
xlab="Time lag", ylab = "Distance", zlab = "STCF", colkey = list(side = 4, length = 0.5),
phi = 20, theta = 120, resfac = 5,  col= gg2.col(100))

## ---- fig.height = 5, warning=FALSE, message=FALSE----------------------------
## set a sequence of hypothetical coordinates 
d <- 5
coord <- cbind(runif(d)*30, runif(d)*30)

## compute VAR model parameters  
fit <- fitVAR(spacepoints = coord, 
              p = 4, 
              margdist = "burrXII",
              margarg = list(scale = 3, shape1 = .9, shape2 = .2), 
              p0 = 0.8, 
              stcsid = "clayton",
              stcsarg = list(scfid = "weibull", tcfid = "weibull", copulaarg = 2,
                   scfarg = list(scale = 25, shape = 0.7), 
                   tcfarg = list(scale = 3.1, shape = 0.8) ) )

## generate correlated timeseries  
sim <- generateMTS(n = 500, STmodel = fit)

## visualize simulated timeseries
dta <- melt(data = data.table(time = 1:nrow(sim), sim[,1:d]), id.vars = "time")

ggplot(data = dta, mapping = aes(x = time, y = value)) + geom_line() +
       facet_grid(facets = variable ~ ., scales = "free_y") + theme_light()

## ---- fig.height = 5, warning=FALSE, message=FALSE----------------------------
## set a sequence of hypothetical coordinates
d <- 5
coord <- cbind(runif(d)*30, runif(d)*30)

## fit and generate correlated timeseries
sim <- generateMTSFast(n = 500, 
                       spacepoints = coord, 
                       p0 = 0.7, 
                       margdist ="burrXII",
                       margarg = list(scale = 3, shape1 = .9, shape2 = .2),  
                       stcsarg = list(scfid = "weibull", tcfid = "weibull",
                       scfarg = list(scale = 25, shape = 0.7),
                       tcfarg = list(scale = 3.1, shape = 0.8)) )

## visualize simulated timeseries
dta <- melt(data = data.table(time = 1:nrow(sim), sim[,1:d]), id.vars = "time")

ggplot(data = dta, mapping = aes(x = time, y = value)) + geom_line() +
       facet_grid(facets = variable ~ ., scales = "free_y") + theme_light()

## ---- warning=FALSE, message=FALSE--------------------------------------------
## compute VAR model parameters
## CPU time: ~15s 
fit <- fitVAR(spacepoints = 20, p = 3, margdist ="burrXII",
              margarg = list(scale = 3, shape1 = .9, shape2 = .2), p0 = 0.8, stcsid = "clayton",
              stcsarg = list(scfid = "weibull", tcfid = "weibull", copulaarg = 2,
              scfarg = list(scale = 20, shape = 0.7), tcfarg = list(scale = 1.1, shape = 0.8) ) )

## generate isotropic random fields with nonseparable correlation structure 
sim1 <- generateRF(n = 1000, STmodel = fit)

## fast simulation of isotropic random fields with separable correlation structure
sim2 <- generateRFFast(n = 1000, spacepoints = 20, p0 = 0.7, margdist ="paretoII",
                       margarg = list(scale = 1, shape = .3),
                       stcsarg = list(scfid = "weibull", tcfid = "weibull",
                       scfarg = list(scale = 20, shape = 0.7), 
                       tcfarg = list(scale = 1.1, shape = 0.8)))


## ---- fig.height = 6.5, warning=FALSE, message=FALSE--------------------------
## check random fields
## CPU time: ~20s 
checkRF(RF = sim1, nfields = 9*9, method = "stat")
checkRF(RF = sim1, nfields = 9*9, method = "statplot")
checkRF(RF = sim1, nfields = 9*9, method = "field")

## ---- fig.width = 4, fig.height = 4, warning=FALSE, message=FALSE-------------
## specify a grid of coordinates
m = 30
aux <- seq(0, m - 1, length = m)
coord <- expand.grid(aux, aux)

## transform the coordinate system
at <- anisotropyTswirl(spacepoints = coord, x0 = floor(m / 2), y0 = floor(m / 2), 
                       b = 10, alpha = 1.8 * pi)

## visualize transformed coordinate system
  aux = data.frame(lon = at[ ,1], lat = at[ ,2], id1 = rep(1:m, each = m), id2 = rep(1:m, m))
  ggplot(aux, aes(x = lon, y = lat)) + 
    geom_path(aes(group = id1)) + geom_path(aes(group = id2)) + geom_point(col = 2) + theme_light()


## ---- fig.height = 7, warning=FALSE, message=FALSE----------------------------
## compute VAR model parameters
fit <- fitVAR(spacepoints = m, p = 1, margdist = 'burrXII', 
              margarg = list(scale = 3, shape1 = .9, shape2 = .2), p0 = 0.1, stcsid = "clayton", 
              stcsarg = list(scfid = "weibull", tcfid = "weibull", copulaarg = 2,
                        scfarg = list(scale = 2, shape = 0.7), tcfarg = list(scale = 20.1, shape = 0.8)),
              anisotropyid = "swirl", 
              anisotropyarg = list(x0 = floor(m / 2), y0 = floor(m / 2), b = 10, alpha = 1.8 * pi) )

## generate
set.seed(9)
sim3 <- generateRF(n=25, STmodel = fit)

## check
checkRF(RF = sim3, nfields = 5*5, method = "field")


## ---- fig.height = 7.0, warning=FALSE, message=FALSE--------------------------
## compute VAR model parameters
fit <- fitVAR(spacepoints = 30, p = 1, margdist = 'burrXII', 
              margarg = list(scale = 3, shape1 = .9, shape2 = .2), p0 = 0.8, stcsid = "clayton", 
              stcsarg = list(scfid = "weibull", tcfid = "weibull", copulaarg = 2, 
                      scfarg = list(scale = 20, shape = 0.7), tcfarg = list(scale = 30.1, shape = 0.8)),
              anisotropyid = "affine", 
              anisotropyarg = list(phi1 = 2., phi2 = 0.5, phi12 = 0, theta = -pi/3),
              advectionid = "uniform", advectionarg = list(u = 2.5, v = 1.5) )

## generate
set.seed(10) # 1 # 5
sim3 <- generateRF(n=81, STmodel = fit)

## check
checkRF(RF = sim3, nfields = 9*9, method = "field")


## ---- warning=FALSE, message=FALSE--------------------------------------------
data("precip")
quickTSPlot(precip$value)

## ---- fig.height = 9, warning=FALSE, message=FALSE----------------------------
## CPU time: ~75s 
precip_ggamma <- analyzeTS(TS = precip, season = "month", dist = "ggamma",
                           acsID = "weibull", lag.max = 12)

reportTS(aTS = precip_ggamma, method = "dist") + theme_light()
reportTS(aTS = precip_ggamma, method = "acs") + theme_light()
reportTS(aTS = precip_ggamma, method = "stat")

## ---- fig.height = 9, warning=FALSE, message=FALSE----------------------------
precip_pareto <- analyzeTS(TS = precip, season = "month", dist = "paretoII", acsID = "fgn", lag.max = 12)

reportTS(aTS = precip_pareto, method = "dist")+ theme_light()
reportTS(aTS = precip_pareto, method = "acs") + theme_light()

## ---- warning=FALSE, message=FALSE--------------------------------------------
sim_precip <- simulateTS(aTS = precip_ggamma, from = as.POSIXct(x = "1978-12-01 00:00:00"),
                         to = as.POSIXct(x = "2008-12-01 00:00:00"))
dta <- precip
dta[, id := "observed"]
sim_precip[, id := "simulated"]
dta <- rbind(dta, sim_precip)

ggplot(data = dta) + geom_line(mapping = aes(x = date, y = value)) + 
  facet_wrap(facets = ~id, ncol = 1) + theme_light()

## ---- fig.height = 9, warning=FALSE, message=FALSE----------------------------
## CPU time: ~240s
data("disch")

str <- analyzeTS(TS = disch, dist = "lnorm", norm = "N2", acsID = "paretoII", 
                 lag.max = 20, constrain = TRUE, season = "month")

reportTS(aTS = str) + theme_light()
reportTS(aTS = str, method = "stat")

## ---- fig.height = 3.5, warning=FALSE, message=FALSE--------------------------
sim_str <- simulateTS(aTS = str)

dta <- disch
dta[, id := "observed"]
sim_str[, id := "simulated"]
dta <- rbind(dta, sim_str)

ggplot(data = dta) + geom_line(mapping = aes(x = date, y = value)) + 
  facet_wrap(facets = ~id, ncol = 1) + theme_light()

