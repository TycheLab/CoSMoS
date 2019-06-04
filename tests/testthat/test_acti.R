test_that('test auto-correlation transformation integral and fit', {

  expect_equal(acti(1, 1, 'norm', list(), .2, 0),             0.0705,
               tolerance = 0.05)

  dist <- 'paretoII'
  distarg <- list(scale = 1, shape = .3)

  pnts <- actpnts(margdist = dist,
                  margarg = distarg,
                  p0 = 0)

  expect_equal(unlist(pnts[c(1, dim(pnts)[1]),]),             c(rhoz = c(0.10, 0.95), rhox = c(0.0549, 0.9134)),
               tolerance = 0.05)

  expect_equal(fitactf(pnts)$`actfcoef`,                      c(b = 1.6154, c = 1.2813),
               tolerance = 0.05)

  expect_error(pnts <- actpnts(margdist = 'gammma', ## misspelled dist
                               margarg = distarg,
                               p0 = 0))
})
