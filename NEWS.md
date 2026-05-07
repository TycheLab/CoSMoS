CoSMoS v2.2.0 (Release date: 2026-05-06)
================

Internal refactoring — no breaking changes to the public API:
----------------
* Homogenised code style throughout (double-quoted strings, no trailing
  `return()`, 2-space indent, consistent `@return` and `@seealso` tags)
* Removed incorrect `@export` from internal helpers across 13 files (fixes
  spurious NAMESPACE entries for `AR1`, `ARp`, `YW`, `N`, `rMSE`, `MSE`,
  `ECDF`, `optimACS`, `lmom`, `erfc`, `inv.erfc`, `seasonalACF`,
  `seasonalAR`, `stratifySeasonData`, `acf*` helpers, and all `pop*`
  functions)
* Merged source files to reduce duplication:
  - `AR1.R` + `ARp.R` -> `ar-models.R`
  - `generateTS.R` -> `generate-ts.R`
  - `generateRF.r` + `generateMTS.r` + `generateRFFast.r` +
    `generateMTSFast.r` -> `generate-fields.R`
  - `ggamma.R` + `paretoII.R` + `burrXII.R` + `burrIII.R` + `gev.R` ->
    `distributions.R`
  - `actf.R` + `actfInv.R` + `acti.R` -> `actf.R`
  - `fitDist.R` + `fitACS.R` + `norm.R` + `errors.R` -> `fitting.R`
  - `analyzeTS.R` + `seasonalACF.R` + `seasonalAR.R` + `stratifyData.R` ->
    `seasonal.R`
  - All S3 plot methods -> `plot-methods.R`
  - Internal utilities -> `utils-internal.R`
* Fixed pre-existing bugs:
  - `generateTS` example used `p = TRUE` instead of `p = NULL` (BUG-06)
  - `regenerateTS` crashed after AR coefficient attribute rename (BUG-07)
  - `regenerateTS` applied AR coefficients in wrong order (BUG-08)
  - `qburrXII` operator-precedence ambiguity (BUG-09)
  - `getACSArg`/`getDistArg` failed when target functions were unexported,
    due to `base::exists()` searching the search path rather than the package
    namespace; fixed to use `asNamespace("CoSMoS")` lookup (BUG-10)
* Extracted shared `simulate_var()` helper from identical `generateRF` /
  `generateMTS` bodies; extracted `DHMgenSj` / `DHMgenSim` from both Fast
  functions; removed dead `stcsid` parameter from `generateRFFast` and
  `generateMTSFast`
* Fixed `simulateTS` to compute probability zero directly from stratified
  data instead of running the full `reportTS` reporting pipeline
* Fixed `rbind`-in-loop growth pattern in `seasonalAR`
* Fixed all `class(x)[1]` comparisons to use `inherits()` or `is.*()` per
  CRAN policy
* Removed deprecated `aes_string()` calls; replaced with `aes()` +
  `after_stat()` / `.data[[]]`
* Removed dead commented-out code throughout

CoSMoS v2.1.0 (Release date: 2021-05-25)
================

Major upgrade to package functionality:
----------------
* Adds Lagrangian Random fields functionality allowing simulation of storms,
  cyclones and other environmental fluxes

CoSMoS v2.0.0 (Release date: 2020-06-23)
================

Major upgrade to package functionality:
----------------
* added random vector and random field simulation

CoSMoS v1.1.3 (Release date: 2020-04-01)
================

Fix for R 4.0.0
----------------

CoSMoS v1.0.1 (Release date: 2019-05-08)
================

Major upgrade to package functionality:
----------------
* added timeseries simulation scheme
* added distribution and autocorrelation structure seasonal fitting

CoSMoS v0.4.1 (Release date: 2019-04-11)
================

First release
----------------
* ability to generate values from marginal distribution with target
  autocorrelation structure
