test_that('test sample l moments', {

  expect_equal(lmom(rnorm(10000,
                          mean = 2,
                          sd = 1)),
               c(l1 = 2.01,
                 l2 = 0.55,
                 l3 = 0.00,
                 l4 = 0.12),
               tolerance = .2)

})
