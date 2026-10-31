# term_block_deriv4() against a central difference of term_block_deriv3()

test_that("the fourth derivative of an nl() block matches the third differenced", {
  skip_if_not_installed("numDeriv")
  set.seed(5)
  dd <- data.frame(x = seq(0.2, 3, length.out = 30), g = rnorm(30))
  b <- term_build(nl(~ a * exp(-r * x) + k * x^2, a ~ 1 + g,
                     links = list(r = linkfunctions7::log_link()),
                     start = list(a = 2, r = 1.3, k = 0.1)), dd)
  m <- length(term_coef_names(b))
  cf <- c(1.8, 0.2, log(1.1), 0.05)
  v <- rnorm(m)
  u <- rnorm(m)
  s <- rnorm(m)
  r <- rnorm(m)
  d4 <- term_block_deriv4(b, coef = cf, v = v, u = u, s = s, r = r)
  num <- matrix(numDeriv::jacobian(function(t)
    as.vector(term_block_deriv3(b, coef = cf + t * r, v = v, u = u, s = s)),
    0), nrow(d4))
  expect_lt(max(abs(d4 - num)) / max(abs(num)), 1e-7)
  expect_equal(d4, term_block_deriv4(b, coef = cf, v = r, u = s, s = v,
                                     r = u), tolerance = 1e-12)
  bl <- term_build(linpar(~ x), dd)
  expect_true(all(term_block_deriv4(bl, v = 1, u = 1, s = 1, r = 1) == 0))
})
