## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(eval    = TRUE,
                      echo    = TRUE,
                      fig.width  = 7,
                      fig.retina = 1,
                      warning = FALSE,
                      message = FALSE)
library(CoSMoS)
library(plot3D)
library(patchwork)
ggplot2::theme_set(ggplot2::theme_light())

## ----warning=FALSE, message=FALSE---------------------------------------------
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

## ----warning=FALSE, message=FALSE---------------------------------------------
acf <- c(1, 0.6, 0.5, 0.4, 0.3) # up to lag-4
ggamma_sim <- generateTS(n = no, margdist = marginaldist, margarg = param, acsvalue = acf)
quickTSPlot(ggamma_sim[[1]])

## ----warning=FALSE, message=FALSE---------------------------------------------
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
  geom_point(size = 2.5) +
  geom_line(lwd = 1) +
  scale_color_manual(
    values = c("steelblue4", "red4", "green4", "darkorange"),
    labels = c("FGN", "Burr XII", "Weibull", "Pareto II"), name = ""
  ) +
  labs(x = bquote(lag ~ tau), y = "ACS") +
  scale_x_continuous(breaks = lags)

## ----warning=FALSE, message=FALSE---------------------------------------------
acf <- acs(id = "paretoII", t = 0:30, scale = 1, shape = .75)
ggamma_sim <- generateTS(n = no, margdist = marginaldist, margarg = param, acsvalue = acf)
dta <- data.frame(time = 1:no, value = ggamma_sim[[1]])

quickTSPlot(dta$value)

## ----warning=FALSE, message=FALSE---------------------------------------------
my_acf <- exp(seq(0, -2, -0.1))
ggamma_sim <- generateTS(n = no, margdist = marginaldist, margarg = param, acsvalue = my_acf)
quickTSPlot(ggamma_sim[[1]])

## ----fig.height = 7, fig.cap = "Five simulated intermittent Generalized Gamma time series (p0 = 0.9).", warning=FALSE, message=FALSE----
prob_zero <- .9
## the argument `TSn = 5` enables the simulation of 5 timeseries.
ggamma_sim <- generateTS(
  n = no, margdist = marginaldist, margarg = param, acsvalue = acf,
  p0 = prob_zero, TSn = 5
)
plot(x = ggamma_sim, main = "") 

## ----warning=FALSE, message=FALSE---------------------------------------------
checkTS(ggamma_sim)

## ----fig.height = 7, fig.cap = "Clayton-Weibull spatiotemporal correlation structure.", warning=FALSE, message=FALSE----
d <- 51
st <- expand.grid(0:(d - 1), 0:(d - 1))

## get the STCS
wc <- stcfclayton(
  t = st[, 1], s = st[, 2], scfid = "weibull", tcfid = "weibull", copulaarg = 2,
  scfarg = list(scale = 20, shape = 0.7), tcfarg = list(scale = 5.1, shape = 0.8)
)

## visualize the STCS
wc.m <- matrix(data = wc, nrow = d)
j <- tail(which(wc.m[1, ] > 0.15), 1)
i <- tail(which(wc.m[, 1] > 0.15), 1)
wc.m <- wc.m[1:i, 1:j]

persp3D(
  z = wc.m, x = 1:nrow(wc.m), y = 1:ncol(wc.m),
  expand = 1, main = "", scale = TRUE, facets = TRUE,
  xlab = "Time lag", ylab = "Distance", zlab = "STCF", colkey = list(side = 4, length = 0.5),
  phi = 20, theta = 120, resfac = 10, col = gg2.col(100)
)

## ----fig.height = 5, fig.cap = "Five spatiotemporally correlated time series simulated with `generateMTS()`.", warning=FALSE, message=FALSE----
## set a sequence of hypothetical coordinates
d <- 5
coord <- cbind(runif(d) * 30, runif(d) * 30)

## compute VAR model parameters
fit <- fitVAR(
  spacepoints = coord,
  p = 4,
  margdist = "burrXII",
  margarg = list(scale = 3, shape1 = .9, shape2 = .2),
  p0 = 0.8,
  stcsid = "clayton",
  stcsarg = list(
    scfid = "weibull", tcfid = "weibull", copulaarg = 2,
    scfarg = list(scale = 25, shape = 0.7),
    tcfarg = list(scale = 3.1, shape = 0.8)
  )
)

## generate correlated timeseries
sim <- generateMTS(n = 500, STmodel = fit)

## visualize simulated timeseries
dta <- melt(data = data.table(time = 1:nrow(sim), sim[, 1:d]), id.vars = "time")

ggplot(data = dta, mapping = aes(x = time, y = value)) +
  geom_line() +
  facet_grid(rows = vars(variable), scales = "free_y")

## ----fig.height = 5, fig.cap = "Five spatiotemporally correlated time series simulated with `generateMTSFast()`.", warning=FALSE, message=FALSE----
## set a sequence of hypothetical coordinates
d <- 5
coord <- cbind(runif(d) * 30, runif(d) * 30)

## fit and generate correlated timeseries
sim <- generateMTSFast(
  n = 500,
  spacepoints = coord,
  p0 = 0.7,
  margdist = "burrXII",
  margarg = list(scale = 3, shape1 = .9, shape2 = .2),
  stcsarg = list(
    scfid = "weibull", tcfid = "weibull",
    scfarg = list(scale = 25, shape = 0.7),
    tcfarg = list(scale = 3.1, shape = 0.8)
  )
)

## visualize simulated timeseries
dta <- melt(data = data.table(time = 1:nrow(sim), sim[, 1:d]), id.vars = "time")

ggplot(data = dta, mapping = aes(x = time, y = value)) +
  geom_line() +
  facet_grid(rows = vars(variable), scales = "free_y")

## ----warning=FALSE, message=FALSE, error=TRUE---------------------------------
try({
## fit VAR model parameters for a 20x20 grid
fit <- fitVAR(spacepoints = 20, p = 3, margdist = "burrXII",
              margarg = list(scale = 3, shape1 = .9, shape2 = .2), p0 = 0.8, stcsid = "clayton",
              stcsarg = list(scfid = "weibull", tcfid = "weibull", copulaarg = 2,
              scfarg = list(scale = 20, shape = 0.7), tcfarg = list(scale = 1.1, shape = 0.8)))

## generate isotropic random fields with nonseparable correlation structure
sim1 <- generateRF(n = 1000, STmodel = fit)

## fast simulation of isotropic random fields with separable correlation structure
sim2 <- generateRFFast(n = 1000, spacepoints = 20, p0 = 0.7, margdist = "paretoII",
                       margarg = list(scale = 1, shape = .3),
                       stcsarg = list(scfid = "weibull", tcfid = "weibull",
                       scfarg = list(scale = 20, shape = 0.7),
                       tcfarg = list(scale = 1.1, shape = 0.8)))
})

## ----fig.height = 6.5, fig.cap = "Validation of simulated random fields: empirical moments (top), correlation functions and CDF (middle), spatial snapshots (bottom).", warning=FALSE, message=FALSE----
## check random fields
checkRF(RF = sim1, nfields = 9*9, method = "stat")
checkRF(RF = sim1, nfields = 9*9, method = "statplot")
checkRF(RF = sim1, nfields = 9*9, method = "field")

## ----fig.width = 4, fig.height = 4, fig.cap = "Swirl-transformed coordinate system (b = 10, α = 1.8π).", warning=FALSE, message=FALSE----
## specify a grid of coordinates
m <- 30
aux <- seq(0, m - 1, length = m)
coord <- expand.grid(aux, aux)

## transform the coordinate system
at <- anisotropyTswirl(
  spacepoints = coord, x0 = floor(m / 2), y0 = floor(m / 2),
  b = 10, alpha = 1.8 * pi
)

## visualize transformed coordinate system
aux <- data.frame(lon = at[, 1], lat = at[, 2], id1 = rep(1:m, each = m), id2 = rep(1:m, m))
ggplot(aux, aes(x = lon, y = lat)) +
  geom_path(aes(group = id1)) +
  geom_path(aes(group = id2)) +
  geom_point(col = 2)

## ----fig.height = 7, fig.cap = "Simulated swirl-like (cyclone-shaped) random fields over a 30×30 grid.", warning=FALSE, message=FALSE----
## compute VAR model parameters
fit <- fitVAR(
  spacepoints = m, p = 1, margdist = "burrXII",
  margarg = list(scale = 3, shape1 = .9, shape2 = .2), p0 = 0.1, stcsid = "clayton",
  stcsarg = list(
    scfid = "weibull", tcfid = "weibull", copulaarg = 2,
    scfarg = list(scale = 2, shape = 0.7), tcfarg = list(scale = 20.1, shape = 0.8)
  ),
  anisotropyid = "swirl",
  anisotropyarg = list(x0 = floor(m / 2), y0 = floor(m / 2), b = 10, alpha = 1.8 * pi)
)

## generate
set.seed(9)
sim3 <- generateRF(n = 25, STmodel = fit)

## check
checkRF(RF = sim3, nfields = 5 * 5, method = "field")

## ----fig.height = 7, fig.cap = "Simulated Lagrangian random fields with affine anisotropy and uniform north-easterly advection.", warning=FALSE, message=FALSE----
## compute VAR model parameters
fit <- fitVAR(
  spacepoints = 30, p = 1, margdist = "burrXII",
  margarg = list(scale = 3, shape1 = .9, shape2 = .2), p0 = 0.8, stcsid = "clayton",
  stcsarg = list(
    scfid = "weibull", tcfid = "weibull", copulaarg = 2,
    scfarg = list(scale = 20, shape = 0.7), tcfarg = list(scale = 30.1, shape = 0.8)
  ),
  anisotropyid = "affine",
  anisotropyarg = list(phi1 = 2., phi2 = 0.5, phi12 = 0, theta = -pi / 3),
  advectionid = "uniform", advectionarg = list(u = 2.5, v = 1.5)
)

## generate
set.seed(10) # 1 # 5
sim3 <- generateRF(n = 81, STmodel = fit)

## check
checkRF(RF = sim3, nfields = 9 * 9, method = "field")

## ----fig.height = 3.5, fig.cap = "Bivariate samples (n = 5000) from Gaussian, Student (df = 4), and flipped Bárdossy copulas with identical linear correlation.", warning=FALSE, message=FALSE----
## set the hypothetical coordinates of two locations
coord <- matrix(rep(c(0, 1.65), times = 2), ncol = 2)

## compute VAR model parameters
## Gaussian dependence structure
fit1 <- fitVAR(
  spacepoints = coord, p = 1, margdist = "norm",
  margarg = list(mean = 0, sd = 1), p0 = 0.0, stcsid = "clayton",
  stcsarg = list(
    scfid = "weibull", tcfid = "weibull", copulaarg = 2,
    scfarg = list(scale = 20, shape = 0.7),
    tcfarg = list(scale = 30.1, shape = 0.8)
  )
)
## Student dependence structure
fit2 <- fitVAR(
  spacepoints = coord, p = 1, margdist = "norm",
  margarg = list(mean = 0, sd = 1), p0 = 0.0, stcsid = "clayton",
  stcsarg = list(
    scfid = "weibull", tcfid = "weibull", copulaarg = 2,
    scfarg = list(scale = 20, shape = 0.7),
    tcfarg = list(scale = 30.1, shape = 0.8)
  ),
  dsid = "student", dsarg = 4
)
## Bardossy dependence structure
fit3 <- fitVAR(
  spacepoints = coord, p = 1, margdist = "norm",
  margarg = list(mean = 0, sd = 1), p0 = 0.0, stcsid = "clayton",
  stcsarg = list(
    scfid = "weibull", tcfid = "weibull", copulaarg = 2,
    scfarg = list(scale = 20, shape = 0.7),
    tcfarg = list(scale = 30.1, shape = 0.8)
  ),
  dsid = "bardossyF", dsarg = 0.5
)
## generate
set.seed(10)
sim1 <- data.frame(generateRF(n = 5000, STmodel = fit1))
sim2 <- data.frame(generateRF(n = 5000, STmodel = fit2))
sim3 <- data.frame(generateRF(n = 5000, STmodel = fit3))

## visual check
sim1 <- cbind(sim1, id = "1 - Gauss copula")
sim2 <- cbind(sim2, id = "2 - Student copula")
sim3 <- cbind(sim3, id = "3 - flipped Bárdossy copula")
dta <- rbind(sim1, sim2, sim3)

ggplot(data = dta) +
  geom_point(mapping = aes(x = X1, y = X2), alpha = 0.07) +
  facet_wrap(facets = ~id, nrow = 1)

## linear correlations
pc <- round(c(
  cor.G = cor(sim1[, 1:2])[1, 2], cor.S = cor(sim2[, 1:2])[1, 2],
  cor.B = cor(sim3[, 1:2])[1, 2]
), 3)
pc

## ----fig.cap = "Observed hourly precipitation record.", warning=FALSE, message=FALSE----
data("precip")
quickTSPlot(precip$value)

## ----eval=FALSE, warning=FALSE, message=FALSE---------------------------------
# precip_ggamma <- analyzeTS(TS = precip, season = "month", dist = "ggamma",
#                            acsID = "weibull", lag.max = 12)

## ----include=FALSE------------------------------------------------------------
precip_ggamma <- readRDS("figures/precip_ggamma.rds")

## ----fig.height = 9, fig.cap = "Fitted Generalized Gamma distribution (top) and Weibull ACS (bottom) against monthly precipitation data.", warning=FALSE, message=FALSE----
reportTS(aTS = precip_ggamma, method = "dist") 
reportTS(aTS = precip_ggamma, method = "acs") 

## ----eval=FALSE, warning=FALSE, message=FALSE---------------------------------
# precip_pareto <- analyzeTS(TS = precip, season = "month", dist = "paretoII", acsID = "fgn", lag.max = 12)

## ----include=FALSE------------------------------------------------------------
precip_pareto <- readRDS("figures/precip_pareto.rds")

## ----fig.height = 9, fig.cap = "Fitted Pareto II distribution (top) and FGN ACS (bottom) against monthly precipitation data: note the poor fit relative to the Generalized Gamma / Weibull combination.", warning=FALSE, message=FALSE----
reportTS(aTS = precip_pareto, method = "dist") 
reportTS(aTS = precip_pareto, method = "acs") 

## ----fig.cap = "Observed (top) and simulated (bottom) hourly precipitation.", warning=FALSE, message=FALSE----
sim_precip <- simulateTS(
  aTS = precip_ggamma, from = as.POSIXct(x = "1978-12-01 00:00:00"),
  to = as.POSIXct(x = "2008-12-01 00:00:00")
)
dta <- precip
dta[, id := "observed"]
sim_precip[, id := "simulated"]
dta <- rbind(dta, sim_precip)

ggplot(data = dta) +
  geom_line(mapping = aes(x = date, y = value)) +
  facet_wrap(facets = ~id, ncol = 1)

## ----eval=FALSE, warning=FALSE, message=FALSE---------------------------------
# data("disch")
# str <- analyzeTS(TS = disch, dist = "lnorm", norm = "N2", acsID = "paretoII",
#                  lag.max = 20, constrain = TRUE, season = "month")

## ----include=FALSE------------------------------------------------------------
data("disch")
str <- readRDS("figures/str.rds")

## ----fig.height = 9, fig.cap = "Fitted log-normal distribution and Pareto II ACS against monthly streamflow data for Nassawango Creek.", warning=FALSE, message=FALSE----
reportTS(aTS = str) 

## ----fig.height = 3.5, fig.cap = "Observed (top) and simulated (bottom) daily streamflow for Nassawango Creek.", warning=FALSE, message=FALSE----
sim_str <- simulateTS(aTS = str)

dta <- disch
dta[, id := "observed"]
sim_str[, id := "simulated"]
dta <- rbind(dta, sim_str)

ggplot(data = dta) +
  geom_line(mapping = aes(x = date, y = value)) +
  facet_wrap(facets = ~id, ncol = 1)

