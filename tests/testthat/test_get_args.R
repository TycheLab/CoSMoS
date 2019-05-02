test_that('test get argument gunction', {

  expect_equal(getDistArg('norm'),
               c('mean', 'sd' ))

  expect_equal(getACSArg('fgn'),
               c('H' ))
})
