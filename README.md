[![Travis build status](https://travis-ci.org/strnda/CoSMoS.svg?branch=master)](https://travis-ci.org/strnda/CoSMoS)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/strnda/CoSMoS?branch=master&svg=true)](https://ci.appveyor.com/project/strnda/CoSMoS)
[![CRAN_Release_Badge](http://www.r-pkg.org/badges/version-ago/CoSMoS)](https://CRAN.R-project.org/package=CoSMoS)
[![CRAN_Download_Badge](http://cranlogs.r-pkg.org/badges/grand-total/CoSMoS)](https://CRAN.R-project.org/package=CoSMoS)
[![license](https://img.shields.io/badge/license-GPL3-lightgrey.svg)](https://choosealicense.com/)


# CoSMoS - Complete Stochastic Modelling Solution
CoSMoS is an R package that makes time series generation with desired properties easy. Just choose the characteristics of the time series you want to generate, and it will do the rest.
The generated time series preserve any probability distribution and any linear autocorrelation structure. Users can generate as many and as long time series from processes such as precipitation, wind, temperature, relative humidity etc. It is based on a framework that unified, extended, and improved a modelling strategy that generates time series by transforming “parent” Gaussian time series having specific characteristics (Papalexiou, 2018).

## Install
To install the latest version of the package run:

```r
## copy-paste to get the latest version of CoSMoS

if (!require('devtools')) {install.packages('devtools'); library(devtools)} 

install_github('strnda/CoSMoS', upgrade = 'never')

library(CoSMoS)

?`CoSMoS-package`
```

# Funding
The package was partly funded by the Global institute for Water Security (GIWS; https://www.usask.ca/water/) and the Global Water Futures (GWF; https://gwf.usask.ca/) program.

# Authors
Coded and maintained by: Filip Strnad    
Conceptual design by: Simon Michael Papalexiou     
Tested and documented by: Yannis Markonis     

# References
Papalexiou, S.M., 2018. Unified theory for stochastic modelling of hydroclimatic processes: Preserving marginal distributions, correlation structures, and intermittency. Advances in Water Resources 115, 234-252. https://doi.org/10.1016/j.advwatres.2018.02.013
