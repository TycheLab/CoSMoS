#' CoSMoS: Complete Stochastic Modelling Solution
#'
#' CoSMoS is an R package that makes time series generation with desired properties easy. Just choose the characteristics of the time series you want to generate, and it will do the rest.
#'
#' The generated time series preserve any probability distribution and any linear autocorrelation structure. Users can generate as many and as long time series from processes such as precipitation, wind, temperature, relative humidity etc. It is based on a framework that unified, extended, and improved a modelling strategy that generates time series by transforming "parent" Gaussian time series having specific characteristics (Papalexiou, 2018).
#'
#' @author \strong{Coded by:} Filip Strnad \email{strnadf@fzp.czu.cz} and Francesco Serinaldi \email{francesco.serinaldi@ncl.ac.uk}
#' @author \strong{Conceptual design by:} Simon Michael Papalexiou \email{sm.papalexiou@usask.ca}
#' @author \strong{Tested and documented by:} Yannis Markonis \email{markonis@fzp.czu.cz}
#' @author \strong{Maintained by:} Kevin Shook \email{kevin.shook@usask.ca}
#'
#' @section Funding:
#'
#' The package was partly funded by the Global institute for Water Security (GIWS; \href{https://water.usask.ca/}{https://water.usask.ca/}) and the Global Water Futures (GWF; \href{https://gwf.usask.ca/}{https://gwf.usask.ca/}) program.
#'
#' @references Papalexiou, S.M. (2018). Unified theory for stochastic modelling of hydroclimatic processes: Preserving marginal distributions, correlation structures, and intermittency. Advances in Water Resources 115, 234-252, \doi{10.1016/j.advwatres.2018.02.013}
#' @references Papalexiou, S.M., Markonis, Y., Lombardo, F., AghaKouchak, A., Foufoula-Georgiou, E. (2018). Precise Temporal Disaggregation Preserving Marginals and Correlations (DiPMaC) for Stationary and Nonstationary Processes. Water Resources Research, 54(10), 7435-7458, \doi{10.1029/2018WR022726}
#' @references Papalexiou, S.M., Serinaldi, F. (2020). Random Fields Simplified: Preserving Marginal Distributions, Correlations, and Intermittency, With Applications From Rainfall to Humidity. Water Resources Research, 56(2), e2019WR026331, \doi{10.1029/2019WR026331}
#' @references Papalexiou, S.M., Serinaldi, F., Porcu, E. (2021). Advancing Space-Time Simulation of Random Fields: From Storms to Cyclones and Beyond. Water Resources Research, 57, e2020WR029466, \doi{10.1029/2020WR029466}
#' @docType package
#' @name CoSMoS-package
NULL
#> NULL
