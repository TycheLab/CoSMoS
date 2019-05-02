test_that('test moments function for normal and pareto type II distributions', {

  expect_equal(moments('norm', list(mean = 2,
                                    sd = 1)),
               list(m = c(m1 = 2, m2 = 5, m3 = 14, m4 = 43),
                    mu = c(mu1 = 2, mu2 = 1, mu3 = 0, mu4 =3),
                    coefficients = c(cv = .5, cs = 0, ck = 3)),
               tolerance = .01)

  expect_equal(moments(dist = 'paretoII',
                       distarg = list(scale = 1,
                                      shape = .3)),
               list(m = c(m1 = 1.4285, m2 = 7.1428, m3 = 214.2857, m4 = -4285.7142),
                    mu = c(mu1 = 1.4285, mu2 = 5.1020, mu3 = 189.5043, mu4 = -5435.2353),
                    coefficients = c(cv = 1.5811, cs = 16.4438, ck = -208.8)),
               tolerance = .01)
})

test_that('test sample.moments function for normal distributions', {

  expect_equal(sample.moments(rnorm(10000,
                                    mean = 2,
                                    sd = 1)),
               list(m = c(m1 = 2, m2 = 5, m3 = 14, m4 = 43),
                    mu = c(mu1 = 2, mu2 = 1, mu3 = 0, mu4 =3),
                    coefficients = c(cv = .5, cs = 0, ck = 3)),
               tolerance = 1)

})
