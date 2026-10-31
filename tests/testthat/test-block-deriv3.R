# term_block_deriv3() against a central difference of term_block_deriv2()

test_that("the third derivative of an nl() block matches the second differenced", {
  skip_if_not_installed("numDeriv")
  set.seed(4)
  dd <- data.frame(x = seq(0.2, 3, length.out = 30), g = rnorm(30))
  dd$y <- 2 * exp(-1.3 * dd$x)
  b <- term_build(nl(~ a * exp(-r * x) + k * x^2, a ~ 1 + g,
                     start = list(a = 2, r = 1.3, k = 0.1)), dd)
  m <- length(term_coef_names(b))
  cf <- c(1.8, 0.2, 1.1, 0.05)
  v <- rnorm(m)
  u <- rnorm(m)
  s <- rnorm(m)
  d3 <- term_block_deriv3(b, coef = cf, v = v, u = u, s = s)
  num <- matrix(numDeriv::jacobian(function(t)
    as.vector(term_block_deriv2(b, coef = cf + t * s, v = v, u = u)), 0),
    nrow(d3))
  expect_lt(max(abs(d3 - num)) / max(abs(num)), 1e-8)
  # symmetric in its three directions
  expect_equal(d3, term_block_deriv3(b, coef = cf, v = s, u = v, s = u),
               tolerance = 1e-12)
  # a fixed block does not move
  bl <- term_build(linpar(~ x), dd)
  expect_true(all(term_block_deriv3(bl, v = 1, u = 1, s = 1) == 0))
})

test_that("a non-identity link enters the third derivative", {
  skip_if_not_installed("numDeriv")
  dd <- data.frame(x = seq(0.2, 3, length.out = 25))
  b <- term_build(nl(~ a * exp(-r * x), r ~ 1,
                     links = list(r = linkfunctions7::log_link()),
                     start = list(a = 2, r = 1.3)), dd)
  cf <- c(1.5, log(0.8))
  set.seed(2)
  v <- rnorm(2)
  u <- rnorm(2)
  s <- rnorm(2)
  d3 <- term_block_deriv3(b, coef = cf, v = v, u = u, s = s)
  num <- matrix(numDeriv::jacobian(function(t)
    as.vector(term_block_deriv2(b, coef = cf + t * s, v = v, u = u)), 0),
    nrow(d3))
  expect_lt(max(abs(d3 - num)) / max(abs(num)), 1e-8)
})
