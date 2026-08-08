# Nonlinear parametric terms: the Jacobian as the design block, the two
# ways of supplying the function, links, and parameter submodels.

set.seed(13)
n <- 60
dd <- data.frame(x = seq(0, 3, length.out = n),
                 g = factor(rep(c("a", "b"), length.out = n)))
dd$y <- 2 * exp(-1.3 * dd$x) + rnorm(n, sd = 0.05)

test_that("the block is the jacobian, and it is exact", {
  spec <- nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3))
  built <- term_build(spec, dd)
  expect_identical(term_coef_names(built), c("nl.a", "nl.r"))
  expect_identical(built@deriv_mode, "symbolic")

  b <- c(2, 1.3)
  # the contribution and its derivative, from the definition
  f <- function(v) v[1] * exp(-v[2] * dd$x)
  expect_equal(term_value(built), f(b), tolerance = 1e-12)
  expect_equal(unname(term_matrix(built)), numDeriv::jacobian(f, b),
               tolerance = 1e-7)
})

test_that("refreshing moves the block to the new parameters", {
  built <- term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3)), dd)
  b2 <- c(0.5, 2.5)
  moved <- term_refresh(built, b2)

  f <- function(v) v[1] * exp(-v[2] * dd$x)
  expect_equal(term_value(moved), f(b2), tolerance = 1e-12)
  expect_equal(unname(term_matrix(moved)), numDeriv::jacobian(f, b2),
               tolerance = 1e-7)
  # the term it came from is untouched
  expect_equal(term_value(built), f(c(2, 1.3)), tolerance = 1e-12)
  expect_error(term_refresh(built, 1), "length 2")

  # and an ordinary term ignores the refresh, its block being data alone
  lin <- term_build(linpar(~x), dd)
  expect_identical(term_matrix(term_refresh(lin, c(0, 0))),
                   term_matrix(lin))
})

test_that("a Gauss-Newton step on the linearization finds the truth", {
  # the whole point of the Jacobian contract: iterate value + J, and the
  # nonlinear least squares problem is solved by linear steps
  built <- term_build(nl(~ a * exp(-r * x), start = list(a = 1, r = 0.5)), dd)
  b <- c(linkfunctions7::linkfun(linkfunctions7::identity_link(), 1), 0.5)
  for (it in 1:25) {
    cur <- term_refresh(built, b)
    J <- term_matrix(cur)
    resid <- dd$y - term_value(cur)
    b <- b + solve(crossprod(J) + 1e-10 * diag(2), crossprod(J, resid))
  }
  expect_equal(as.numeric(b), c(2, 1.3), tolerance = 0.05)
})

test_that("a link holds a parameter in its own set", {
  spec <- nl(~ a * exp(-r * x),
             links = list(r = linkfunctions7::log_link()),
             start = list(a = 2, r = 1.3))
  built <- term_build(spec, dd)
  # the coefficient is now log(r), so the block carries the chain rule
  b <- c(2, log(1.3))
  f <- function(v) v[1] * exp(-exp(v[2]) * dd$x)
  expect_equal(term_value(built), f(b), tolerance = 1e-10)
  expect_equal(unname(term_matrix(built)), numDeriv::jacobian(f, b),
               tolerance = 1e-6)
  # and no coefficient can produce a negative rate
  expect_gt(min(exp(seq(-50, 50, length.out = 11))), 0)
})

test_that("a parameter can be developed with covariates", {
  spec <- nl(~ a * exp(-r * x), subformulas = list(a = ~g),
             start = list(r = 1.3))
  built <- term_build(spec, dd)
  expect_identical(term_coef_names(built),
                   c("nl.a.(Intercept)", "nl.a.gb", "nl.r"))

  Z <- stats::model.matrix(~g, dd)
  b <- c(2, -0.4, 1.3)
  f <- function(v) as.numeric(Z %*% v[1:2]) * exp(-v[3] * dd$x)
  moved <- term_refresh(built, b)
  expect_equal(term_value(moved), f(b), tolerance = 1e-10)
  expect_equal(unname(term_matrix(moved)), numDeriv::jacobian(f, b),
               tolerance = 1e-6)

  res <- check_term(spec, dd, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
})

test_that("a function is accepted and differenced, and says so", {
  fx <- function(x, theta) theta$a * exp(-theta$r * x)
  spec <- nl(fx, params = c("a", "r"), x = x, start = list(a = 2, r = 1.3))
  built <- term_build(spec, dd)
  expect_identical(built@deriv_mode, "numeric")

  b <- c(2, 1.3)
  f <- function(v) v[1] * exp(-v[2] * dd$x)
  expect_equal(term_value(built), f(b), tolerance = 1e-12)
  # differenced, so held to the accuracy of a central difference
  expect_equal(unname(term_matrix(built)), numDeriv::jacobian(f, b),
               tolerance = 1e-6)

  res <- check_term(spec, dd, verbose = FALSE)
  expect_true(all(res$status == "OK"))
})

test_that("an expression deriv() cannot read falls back to differences", {
  # besselJ has no symbolic derivative in R's table, so the formula route
  # must notice and difference instead of failing
  spec <- nl(~ a * besselJ(x + 1, 0) + r * x, start = list(a = 1, r = 1))
  built <- term_build(spec, dd)
  expect_identical(built@deriv_mode, "formula_fd")
  b <- c(1.5, 0.4)
  f <- function(v) v[1] * besselJ(dd$x + 1, 0) + v[2] * dd$x
  moved <- term_refresh(built, b)
  expect_equal(term_value(moved), f(b), tolerance = 1e-12)
  expect_equal(unname(term_matrix(moved)), numDeriv::jacobian(f, b),
               tolerance = 1e-6)
})

test_that("what only a formula can do is refused to a function", {
  fx <- function(x, theta) theta$a * x
  expect_error(nl(fx), "'params' must name the parameters")
  expect_error(nl(fx, params = "a", subformulas = list(a = ~g)),
               "only a formula says")
  expect_error(nl(y ~ a * x), "one-sided")
  expect_error(nl("not a formula"), "formula or a function")
  expect_error(nl(~ a * x, links = list(nope = 1)), NA)
  expect_error(term_build(nl(~ a * x, links = list(nope = 1)), dd),
               "not a parameter")
  expect_error(term_build(nl(~ x + 1), dd), "no parameters to estimate")
})

test_that("the interpreter routes it and print says how it differentiates", {
  out <- interpret_formula(y ~ nl(~ a * exp(-r * x)), dd)
  expect_named(out$terms, c("linpar", "nl(~a * exp(-r * x))"))
  built <- term_build(out$terms[["nl(~a * exp(-r * x))"]], dd)
  expect_output(print(built), "symbolic derivatives")
  expect_output(print(nl(~ a * x)), "specification")
})
