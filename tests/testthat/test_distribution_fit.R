test_that('test distribution fitting', {

  expect_equal({
    x <- fitDist(rnorm(10000),
                 dist = 'norm',
                 n.points = 40,
                 norm = 'N2',
                 constrain = FALSE)
    attributes(x) <- NULL
    x
  },
  list(0, 1),
  tolerance = .5)
})
