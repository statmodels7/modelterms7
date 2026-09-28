test_that("term_charted names the coordinates that ride a chart with an edge", {
  expect_identical(term_charted(regime(2)), c("gap2", "alr1.1", "alr2.1"))
  expect_identical(term_charted(gas(p = 1, q = 1)), c("alpha1", "pacf1"))
  expect_identical(term_charted(linpar(~ x)), character(0))
  set.seed(1)
  dd <- data.frame(x = runif(50), g = factor(rep(1:5, 10)))
  tm <- term_build(nl(~ a * exp(-r * x) + c,
                      links = list(r = linkfunctions7::log_link(),
                                   c = linkfunctions7::logit_link())), dd)
  expect_identical(term_charted(tm), c(r = 2L, c = 3L))
  # a parameter with a subformula is not a scalar coordinate
  tm2 <- term_build(nl(~ a * exp(-r * x), a ~ 1 + ridge(~ g),
                       links = list(a = linkfunctions7::log_link(),
                                    r = linkfunctions7::log_link())), dd)
  expect_identical(term_charted(tm2), c(r = 7L))
})

test_that("a regime term starts its levels at the quantiles of the response", {
  set.seed(2)
  y <- c(rnorm(100, 55, 6), rnorm(100, 80, 6))
  s <- term_start(regime(2), target = y)
  q <- stats::quantile(y, c(0.25, 0.75), names = FALSE)
  expect_equal(s[["level1"]], q[1])
  expect_equal(exp(s[["gap2"]]), q[2] - q[1])
  expect_equal(unname(s[c("alr1.1", "alr2.1")]), c(0, 0))
  s3 <- term_start(regime(3), target = y)
  q3 <- stats::quantile(y, c(1, 3, 5) / 6, names = FALSE)
  expect_equal(exp(unname(s3[c("gap2", "gap3")])), diff(q3))
  # without a response, or with coinciding quantiles, the conventional start
  expect_true(all(term_start(regime(2)) == 0))
  expect_true(all(term_start(regime(2), target = rep(1, 20)) == 0))
})
