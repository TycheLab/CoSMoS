test_that('test random values deneration using AR1 and ARp models', {

  expect_equal({set.seed(666);
    AR1(n = 10,
        alpha = .8,
        mean = 0,
        sd = 1)},
    c(0.7533, 0.3895, 1.5285, -0.1072, 0.3692, -0.4883, -0.8721, -1.7730, -1.4436, 0.1350),
    tolerance = .01)


  dist <- 'paretoII'
  distarg <- list(scale = 1,
                  shape = .3)
  p0 <- .5
  order <- 5000

  pnts <- actpnts(margdist = dist,
                  margarg = distarg,
                  p0 = p0)
  fit <- fitactf(pnts)

  acsvalue <- acs(id = 'weibull',
                  t = 0:order,
                  scale = 10,
                  shape = .75)

  set.seed(666)
  val <- ARp(margdist = dist,
             margarg = distarg,
             acsvalue = acsvalue,
             actfpara = fit,
             n = 10000,
             p0 = p0)
  attributes(val) <- NULL

  expect_equal(val[1:10],
               c(0.1333, 0.4907, 0.8089, 0.2326, 0, 0, 0.0967, 1.1822, 1.9728, 1.2610),
    tolerance = .01)
})

test_that('test wrapper for acti and ARp functions', {
  set.seed(666)
  x <- generateTS(margdist = 'burrXII',
                  margarg = list(scale = 1,
                                 shape1 = .75,
                                 shape2 = .25),
                  acsvalue = acs(id = 'weibull',
                                 t = 0:30,
                                 scale = 10,
                                 shape = .75),
                  n = 1000, p = 30, p0 = .5, TSn = 3)

  expect_equal({val <- x[[1]];
    attributes(val) <- NULL;
    val[190:200]},
  c(0, 0, 0, 0, 0, 0, 0.1546, 0, 0, 0, 0),
  tolerance = .01)

  expect_equal({set.seed(666);
    val <- regenerateTS(x);
    unlist(val)[190:200]},
  c(0, 0, 0, 0, 0.0988, 0, 0, 0, 0, 0, 0.2773),
  tolerance = .01)
})
