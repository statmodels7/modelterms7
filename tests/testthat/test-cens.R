# The censored-response constructor.

test_that("statuses follow from the values and the bounds", {
  r <- cens(c(0, 0.7, 2.4, 5, NA), lwr = 0, upr = c(5, 5, 5, 5, 4))
  expect_identical(r@status,
                   c("left", "observed", "observed", "right", "interval"))
  expect_identical(r@y[1L], 0)
  expect_identical(r@lwr, rep(0, 5))
})

test_that("scalar bounds recycle and defaults censor nothing", {
  r <- cens(c(-1, 0, 1))
  expect_identical(r@status, rep("observed", 3))
  expect_identical(r@lwr, rep(-Inf, 3))
  expect_identical(r@upr, rep(Inf, 3))
})

test_that("degenerate inputs are rejected", {
  expect_error(cens(1:3, lwr = c(0, 1)), "length")
  expect_error(cens(1:3, lwr = 2, upr = 2), "strictly below")
  expect_error(cens(c(1, NA), lwr = 0), "finite")
  expect_error(cens(1:3, lwr = NA), "NA")
})

test_that("print counts the statuses", {
  expect_output(print(cens(c(0, 1, 2), lwr = 0)),
                "3 observations: 2 observed, 1 left")
})
