# The derivatives of a nonlinear term in its own parameters.

nld_data <- function(n = 40) {
  set.seed(5)
  d <- data.frame(x = seq(0.2, 3, length.out = n))
  d$y <- 2 * exp(-1.3 * d$x) + stats::rnorm(n, sd = 0.05)
  d
}

# the closed forms of f = a exp(-r x), transcribed by hand, which is what makes
# them a reference rather than the same arithmetic twice
nld_truth <- function(a, r, x) {
  e <- exp(-r * x)
  list(a = e, r = -a * x * e,
       a_a = 0 * e, a_r = -x * e, r_r = a * x^2 * e,
       a_a_a = 0 * e, a_a_r = 0 * e, a_r_r = x^2 * e, r_r_r = -a * x^3 * e,
       a_a_a_a = 0 * e, a_a_a_r = 0 * e, a_a_r_r = 0 * e,
       a_r_r_r = -x^3 * e, r_r_r_r = a * x^4 * e)
}

nld_worst <- function(got, truth) {
  max(vapply(names(got), function(nm) {
    t <- truth[[nm]]
    max(abs(as.numeric(got[[nm]]) - t)) / max(1e-8, max(abs(t)))
  }, numeric(1)))
}

nld_g <- function(theta, data) {
  e <- exp(-theta$r * data$x)
  list(a = e, r = -theta$a * data$x * e)
}
nld_h <- function(theta, data) {
  e <- exp(-theta$r * data$x)
  list(a_a = 0 * e, r_r = theta$a * data$x^2 * e, a_r = -data$x * e)
}

test_that("a formula is symbolic at every order", {
  d <- nld_data()
  b <- term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3)), d)
  tr <- nld_truth(2, 1.3, d$x)
  for (k in 1:4) expect_lt(nld_worst(nl_fderiv(b, order = k), tr), 1e-14)
})

test_that("supplied orders are used and the rest fall back", {
  # An OPAQUE function has nothing symbolic, so this is where the ladder shows.
  # Measured against the closed forms: with nothing supplied a single stencil
  # of order k on the value gives 8.4e-03, 4.62 and 1.96e+03 at orders 2, 3
  # and 4; with the gradient and the Hessian written out, 2.2e-12 and 8.9e-11.
  d <- nld_data()
  tr <- nld_truth(2, 1.3, d$x)
  fo <- function(x, theta) theta$a * exp(-theta$r * x$x)
  st <- list(a = 2, r = 1.3)

  bare <- term_build(nl(fo, params = c("a", "r"), start = st), d)
  expect_gt(nld_worst(nl_fderiv(bare, order = 3), tr), 1e-2)

  full <- term_build(nl(fo, params = c("a", "r"), start = st,
                        gradient = nld_g, hessian = nld_h), d)
  expect_lt(nld_worst(nl_fderiv(full, order = 1), tr), 1e-14)
  expect_lt(nld_worst(nl_fderiv(full, order = 2), tr), 1e-14)
  # the orders NOT supplied are one difference from the exact Hessian, which
  # is the whole point of letting the two be given separately
  expect_lt(nld_worst(nl_fderiv(full, order = 3), tr), 1e-9)
  expect_lt(nld_worst(nl_fderiv(full, order = 4), tr), 1e-8)
})

test_that("a supplied gradient is the design block", {
  d <- nld_data()
  st <- list(a = 2, r = 1.3)
  b1 <- term_build(nl(~ a * exp(-r * x), start = st), d)
  b2 <- term_build(nl(~ a * exp(-r * x), start = st, gradient = nld_g), d)
  expect_equal(as.matrix(term_matrix(b1)), as.matrix(term_matrix(b2)),
               tolerance = 1e-12)
})

test_that("the component names are normalized, not parsed", {
  # 'r_a' is 'a_r', and the order they are returned in does not matter: the
  # accepted spellings are BUILT from the parameter order, so a parameter whose
  # own name contains an underscore cannot be mis-read -- the trap this package
  # records for Hessian component names.
  d <- nld_data()
  swapped <- function(theta, data) {
    e <- exp(-theta$r * data$x)
    list(r_a = -data$x * e, a_a = 0 * e, r_r = theta$a * data$x^2 * e)
  }
  b <- term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3),
                     hessian = swapped), d)
  tr <- nld_truth(2, 1.3, d$x)
  expect_lt(nld_worst(nl_fderiv(b, order = 2), tr), 1e-14)
})

test_that("a wrong name is an error and not a silent fall-back", {
  d <- nld_data()
  st <- list(a = 2, r = 1.3)
  wrong <- function(theta, data) {
    list(a_a = data$x, a_q = data$x, r_r = data$x)
  }
  expect_error(term_build(nl(~ a * exp(-r * x), start = st, hessian = wrong), d),
               "not one of this term's")
  short <- function(theta, data) list(a_a = data$x)
  expect_error(term_build(nl(~ a * exp(-r * x), start = st, hessian = short), d),
               "did not return the component")
  unnamed <- function(theta, data) list(data$x, data$x, data$x)
  expect_error(term_build(nl(~ a * exp(-r * x), start = st,
                             hessian = unnamed), d),
               "NAMED list")
  expect_error(nl(~ a * exp(-r * x), gradient = 42), "must be a function")
  expect_error(nl(~ a * exp(-r * x), hessian = function(z) z),
               "must take \\(theta, data\\)")
})

test_that("the derivatives are in the parameters, and the links are ours", {
  # a log link on r means the coefficient is log(r); the components this
  # returns are in r itself, the chain rule onto the coefficient belonging to
  # the term
  d <- nld_data()
  b <- term_build(nl(~ a * exp(-r * x),
                     links = list(r = linkfunctions7::log_link()),
                     start = list(a = 2, r = 1.3), gradient = nld_g), d)
  tr <- nld_truth(2, 1.3, d$x)
  expect_lt(nld_worst(nl_fderiv(b, order = 1), tr), 1e-14)
  # and the block carries the chain rule, so its r column is scaled by r
  J <- as.matrix(term_matrix(b))
  expect_equal(J[, 2], tr$r * 1.3, tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("nl_fderiv reports what it cannot do", {
  d <- nld_data()
  b <- term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3)), d)
  expect_error(nl_fderiv(b, order = 5), "must be 1, 2, 3 or 4")
  expect_error(nl_fderiv(nl(~ a * exp(-r * x)), order = 1), "not built")
})
