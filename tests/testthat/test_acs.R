test_that('test ACS generation', {

  t <- 2
  H <- .75
  scale <- 1
  shape <- .3
  shape1 <- .6
  shape2 <- .2

  expect_equal(acs('fgn', t, H),                         0.27,
               tolerance = 0.01)

  expect_equal(acs('burrXII', t, scale, shape1, shape2),   0.916,
               tolerance = 0.01)

  expect_equal(acs('paretoII', t, scale, shape),           0.209,
               tolerance = 0.01)

  expect_equal(acs('weibull', t, scale, shape),            0.292,
               tolerance = 0.01)

  expect_error(acs('wleibull', t, scale, shape)) ## misspelled distribution name
})
