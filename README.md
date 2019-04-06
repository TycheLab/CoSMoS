# CoSMoS - Complete Stochastic Modelling Solution
CoSMoS is an R package that makes time series generation with desired properties easy. Just choose the characteristics of the time series you want to generate, and it will do the rest.
The generated time series preserve any probability distribution and any linear autocorrelation structure. Users can generate as many and as long time series from processes such as precipitation, wind, temperature, relative humidity etc. It is based on a framework that unified, extended, and improved a modelling strategy that generates time series by transforming “parent” Gaussian time series having specific characteristics (Papalexiou, 2018).

## Install
To install the latest version of the package run:

```r
if (!require('devtools')) {install.packages('devtools'); library(devtools)} 

install_github('strnda/CoSMoS', upgrade = 'never')

library(CoSMoS)
```

# Funding
The package was partly funded by the Global institute for Water Security (GIWS; https://www.usask.ca/water/) and the Global Water Futures (GWF; https://gwf.usask.ca/) program.

# Authors
Coded and maintained by: Filip Strnad    
Conceptual design by: Simon Michael Papalexiou, and Filip Strnad    
Beta tested by: Filip Strnad, Yannis Markonis, and Simon Michael Papalexiou    

# References
Papalexiou, S.M., 2018. Unified theory for stochastic modelling of hydroclimatic processes: Preserving marginal distributions, correlation structures, and intermittency. Advances in Water Resources 115, 234-252. https://doi.org/10.1016/j.advwatres.2018.02.013
