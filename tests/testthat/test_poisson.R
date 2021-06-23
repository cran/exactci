

context("Poisson CIs")


test_that("Poisson CIs x=1,conf.level<0.5",{
  # minlike
  expect_equal(
    as.vector(round(
      exactpoissonCI(x=1, tsmethod = "minlike", conf.level=0.366),4)),
      c(1.0051,2.4495)
    )
  # blaker
  expect_equal(
    as.vector(round(
      exactpoissonCI(x=1, tsmethod="blaker", conf.level=0.366),4)),
    c(0.6931,2.1559)
  )


})

test_that("Poisson x=100 conf.level=0.95",{
  # minlike
  expect_equal(
    as.vector(round(
      exactpoissonCI(x=100, tsmethod = "minlike"),4)),
    c(81.8391, 121.8372)
  )
  # blaker
  expect_equal(
    as.vector(round(
      exactpoissonCI(x=100, tsmethod="blaker"),4)),
    c(81.5387, 121.5313)
  )
})
