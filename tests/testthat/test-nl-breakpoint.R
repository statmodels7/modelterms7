# A parameter of nl() developed over a break-point term: the sub-term's block
# moves with its own coefficients, so the parameter's design and predictor are
# read at the coefficients. Every order is checked against ONE numerical
# differentiation of the analytic order below it.

nlbp_data <- function(seed = 11, n = 120) {
  set.seed(seed)
  data.frame(t = seq(0, 10, length.out = n), x = stats::runif(n, 0, 3),
             id = factor(rep(1:4, each = n / 4)))
}
nlbp_specs <- function() {
  pr <- numericals7::smooth_probit()
  list(
    jump_smoothed = nl(~ a * exp(-r * x), r ~ jump(t, smoothed = pr),
                       start = list(a = 4, r = 0.5)),
    jseg_smoothed = nl(~ a * exp(-r * x), r ~ jseg(t, smoothed = pr),
                       start = list(a = 4, r = 0.5)),
    seg_sharp = nl(~ a * exp(-r * x), r ~ seg(t),
                   start = list(a = 4, r = 0.5)),
    developed = nl(~ a * exp(-r * x),
                   r ~ jump(t, psi ~ random(~1 | id), smoothed = pr),
                   a ~ id, start = list(a = 4, r = 0.5)),
    both = nl(~ a * exp(-r * x), r ~ jump(t, smoothed = pr),
              a ~ seg(t, smoothed = pr), start = list(a = 4, r = 0.5)))
}

test_that("a sharp jump or jseg is rejected by name, pointing at smoothed", {
  d <- nlbp_data()
  expect_error(term_build(nl(~ a * exp(-r * x), r ~ jump(t),
                             start = list(a = 4, r = 0.5)), d),
               "smoothed =")
  expect_error(term_build(nl(~ a * exp(-r * x), r ~ jseg(t),
                             start = list(a = 4, r = 0.5)), d),
               "smoothed =")
  # the admission is nl()'s alone: a break-point developing a break-point
  # term's own coefficient is still a fixed-design question
  expect_error(term_build(seg(x, psi ~ seg(t)), d), "fixed design")
})

test_that("the block and its derivatives are exact at every order", {
  skip_if_not_installed("numDeriv")
  d <- nlbp_data()
  n <- nrow(d)
  for (nm in names(nlbp_specs())) {
    b <- term_build(nlbp_specs()[[nm]], d)
    m <- length(b@coef_names)
    set.seed(2)
    cf <- b@blueprint$coef + stats::rnorm(m, sd = 0.05)
    v <- stats::rnorm(m)
    u <- stats::rnorm(m)
    A <- matrix(stats::rnorm(n * m), n, m)
    blk <- function(c) as.matrix(term_matrix(term_refresh(b, c)))
    # order one: the block is the derivative of the contribution
    J <- numDeriv::jacobian(function(c) term_value(b, c), cf)
    expect_equal(blk(cf), J, tolerance = 1e-7, ignore_attr = TRUE,
                 label = paste(nm, "block"))
    # order two along v, and its adjoint
    dX <- matrix(numDeriv::jacobian(function(s) as.numeric(blk(cf + s * v)),
                                    0), n, m)
    expect_equal(term_block_deriv(b, cf, v), dX, tolerance = 1e-7,
                 label = paste(nm, "deriv"))
    gc <- numDeriv::grad(function(c) sum(A * blk(c)), cf)
    expect_equal(term_block_contract(b, cf, A), gc, tolerance = 1e-7,
                 label = paste(nm, "contract"))
    # order three, against one difference of the analytic order two
    d2 <- matrix(numDeriv::jacobian(function(s)
      as.numeric(term_block_deriv(b, cf + s * u, v)), 0), n, m)
    b2 <- term_block_deriv2(b, cf, v, u)
    expect_equal(b2, d2, tolerance = 1e-6, label = paste(nm, "deriv2"))
    expect_equal(b2, term_block_deriv2(b, cf, u, v), tolerance = 1e-12)
  }
})

test_that("a nested break-point is read the same way on new rows", {
  d <- nlbp_data()
  b <- term_build(nlbp_specs()$developed, d)
  cf <- b@blueprint$coef + 0.03
  r <- term_refresh(b, cf)
  expect_equal(term_value(r, newdata = d), term_value(r), tolerance = 1e-12)
  expect_equal(as.matrix(term_predict(r, d)), as.matrix(term_matrix(r)),
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_length(term_value(r, newdata = d[1:5, ]), 5L)
})

test_that("the position is confined to the sub-term's interval", {
  d <- nlbp_data()
  b <- term_build(nlbp_specs()$jump_smoothed, d)
  pi <- grep("psi", b@coef_names)
  lim <- term_components(b)$r$subs[[2L]]@blueprint$lim
  far <- b@blueprint$coef
  at <- far
  far[pi] <- 1e6
  at[pi] <- lim[2L]
  expect_equal(term_value(b, far), term_value(b, at), tolerance = 1e-12)
})

test_that("the start search finds the break-point from the response", {
  # the smoothed jseg inside a log-rate, where a single start lands on a
  # local optimum more often than not
  set.seed(12)
  n <- 600
  d <- data.frame(t = seq(0, 10, length.out = n), x = stats::runif(n, 0, 3))
  lr <- log(0.5) + 0.02 * d$t + 0.25 * pmax(d$t - 6, 0) + 0.5 * (d$t > 6)
  d$y <- 5 * exp(-exp(lr) * d$x) + stats::rnorm(n, sd = 0.2)
  b <- term_build(nl(~ a * exp(-r * x),
                     r ~ jseg(t, smoothed = numericals7::smooth_probit(h = 0.05)),
                     links = list(r = linkfunctions7::log_link()),
                     start = list(a = 4)), d)
  cs <- term_coef_start(b, target = d$y)
  expect_lt(abs(cs[grep("psi", b@coef_names)] - 6), 0.2)
  # the refreshed sub-term reports the position through seg_psi()
  r <- term_refresh(b, cs)
  expect_equal(as.numeric(seg_psi(term_components(r)$r$subs[[2L]])),
               as.numeric(cs[grep("psi", b@coef_names)]))
})

test_that("a sharp seg's kink is carried to the term's coefficients", {
  d <- nlbp_data()
  b <- term_build(nlbp_specs()$seg_sharp, d)
  cf <- b@blueprint$coef
  pi <- grep("psi", b@coef_names)
  cf[pi] <- d$t[60]
  expect_identical(term_kinks(b, cf), pi)
  cf[pi] <- d$t[60] + 0.3 * diff(d$t[60:61])
  expect_identical(term_kinks(b, cf), integer(0))
})
