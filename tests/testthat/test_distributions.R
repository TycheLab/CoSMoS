test_that('test burr type III dist', {

  scale <- 1
  shape1 <- .8
  shape2 <- .2

  expect_equal(dburrIII(1, scale, shape1, shape2),      0.39,
               tolerance = 0.01)

  expect_equal(pburrIII(1, scale, shape1, shape2),      0.878,
               tolerance = 0.01)

  set.seed(666)
  expect_equal(rburrIII(1, scale, shape1, shape2),      1.45,
               tolerance = 0.01)

  expect_equal(mburrIII(1, scale, shape1, shape2),      0.518,
               tolerance = 0.01)
})

test_that('test burr type XII dist', {

  scale <- 1
  shape1 <- .8
  shape2 <- .2

  expect_equal(dburrXII(1, scale, shape1, shape2),      0.267,
               tolerance = 0.01)

  expect_equal(pburrXII(1, scale, shape1, shape2),      0.68,
               tolerance = 0.01)

  set.seed(666)
  expect_equal(rburrXII(1, scale, shape1, shape2),      1.45,
               tolerance = 0.01)

  expect_equal(mburrXII(1, scale, shape1, shape2),      1.1,
               tolerance = 0.01)
})

test_that('test generalized gamma dist', {

  scale <- 1
  shape1 <- .8
  shape2 <- .2

  expect_equal(dggamma(1, scale, shape1, shape2),      0.0273,
               tolerance = 0.01)

  expect_equal(pggamma(1, scale, shape1, shape2),      0.019,
               tolerance = 0.01)

  set.seed(666)
  expect_equal(rggamma(1, scale, shape1, shape2),      4175,
               tolerance = 0.01)

  expect_equal(mggamma(1, scale, shape1, shape2),      6720,
               tolerance = 0.01)
})

test_that('test pareto type II dist', {

  scale <- 1
  shape <- .3

  expect_equal(dparetoII(1, scale, shape),      0.32,
               tolerance = 0.01)

  expect_equal(pparetoII(1, scale, shape),      0.58,
               tolerance = 0.01)

  set.seed(666)
  expect_equal(rparetoII(1, scale, shape),      1.88,
               tolerance = 0.01)

  expect_equal(mparetoII(1, scale, shape),      1.43,
               tolerance = 0.01)
})
