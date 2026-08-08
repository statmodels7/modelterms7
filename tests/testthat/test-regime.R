# Markov regimes: the forward recursion, its exact derivative, and the
# charts that keep the levels ordered and the rows stochastic.

set.seed(17)
n <- 50
dd <- data.frame(t = seq_len(n),
                 y = c(rnorm(25, 0), rnorm(25, 3)),
                 g = rep(c("a", "b"), each = 25))

gauss_ld <- function(y) function(e, i) stats::dnorm(y[i], e, log = TRUE)
gauss_sc <- function(y) function(e, i) y[i] - e

psi2 <- list(level1 = 0.1, gap2 = 2.8, alr1.1 = 2, alr2.1 = -2)

test_that("the parameters carry ordered levels and a stochastic chain", {
  expect_identical(term_params(regime(2)),
                   c("level1", "gap2", "alr1.1", "alr2.1"))
  # one free level, k - 1 positive gaps, and k(k - 1) log-ratios
  k <- 3L
  expect_identical(length(term_params(regime(k))), 1L + (k - 1L) + k * (k - 1L))
  lk <- term_links(regime(2))
  expect_identical(vapply(lk, function(l) l@link_name, character(1)),
                   c(level1 = "identity", gap2 = "log",
                     alr1.1 = "identity", alr2.1 = "identity"))
  expect_error(regime(1), "at least 2")
  expect_error(regime(2, label = ""), "non-empty")
})

test_that("the log-likelihood is the mixture over paths, summed by hand", {
  term <- term_build(regime(2, time = t), dd)
  out <- term_loglik(term, rep(0, n), dd$y, gauss_ld(dd$y), gauss_sc(dd$y),
                     psi2)

  # the same recursion written independently, unnormalized on a short
  # series where the product is still representable
  short <- 12
  P <- parameters7::param_value(parameters7::transition_matrix(2),
                                c(psi2$alr1.1, psi2$alr2.1))
  A <- diag(2) - P
  A[, 2] <- 1
  delta <- as.numeric(solve(t(A), c(0, 1)))
  mu <- c(psi2$level1, psi2$level1 + psi2$gap2)

  al <- delta * stats::dnorm(dd$y[1], mu)
  for (t in 2:short) {
    al <- as.numeric(al %*% P) * stats::dnorm(dd$y[t], mu)
  }
  expect_equal(sum(out$loglik[1:short]), log(sum(al)), tolerance = 1e-10)

  # and the normalization is what makes a long series representable at
  # all. The unnormalized product decays geometrically: measured on this
  # setup it is 4e-142 by the two hundredth observation and exactly zero
  # by the four hundred and seventy-first.
  set.seed(21)
  long <- c(stats::rnorm(300, 0), stats::rnorm(300, 3))
  al_long <- delta * stats::dnorm(long[1], mu)
  for (t in 2:600) {
    al_long <- as.numeric(al_long %*% P) * stats::dnorm(long[t], mu)
  }
  expect_true(all(al_long == 0))

  dlong <- data.frame(t = seq_along(long), y = long)
  tlong <- term_build(regime(2, time = t), dlong)
  out_long <- term_loglik(tlong, rep(0, 600), long, gauss_ld(long),
                          gauss_sc(long), psi2)
  expect_true(is.finite(sum(out_long$loglik)))
})

test_that("the jacobian of the recursion is exact", {
  for (k in c(2L, 3L)) {
    term <- term_build(regime(k, time = t), dd)
    nm <- term_params(term)
    v0 <- stats::setNames(
      c(0.1, rep(1.5, k - 1L), rep(0.3, k * (k - 1L))), nm)

    run <- function(v) {
      sum(term_loglik(term, rep(0, n), dd$y, gauss_ld(dd$y), gauss_sc(dd$y),
                      as.list(stats::setNames(v, nm)))$loglik)
    }
    got <- colSums(term_loglik(term, rep(0, n), dd$y, gauss_ld(dd$y),
                               gauss_sc(dd$y),
                               as.list(v0))$jacobian)
    expect_equal(unname(got), numDeriv::grad(run, v0), tolerance = 1e-6,
                 info = sprintf("k = %d", k))
  }
})

test_that("the static predictor enters every regime", {
  term <- term_build(regime(2, time = t), dd)
  base <- term_loglik(term, rep(0, n), dd$y, gauss_ld(dd$y), gauss_sc(dd$y),
                      psi2)
  # shifting the predictor by c is the same as shifting every level by c
  shifted <- term_loglik(term, rep(0.7, n), dd$y, gauss_ld(dd$y),
                         gauss_sc(dd$y), psi2)
  moved <- term_loglik(term, rep(0, n), dd$y, gauss_ld(dd$y), gauss_sc(dd$y),
                       utils::modifyList(psi2, list(level1 = 0.8)))
  expect_equal(shifted$loglik, moved$loglik, tolerance = 1e-10)
})

test_that("the fit recovers regimes it was given", {
  set.seed(4)
  m <- 400
  P <- matrix(c(0.95, 0.05, 0.08, 0.92), 2, 2, byrow = TRUE)
  s <- integer(m); s[1] <- 1L
  for (t in 2:m) s[t] <- sample(1:2, 1, prob = P[s[t - 1], ])
  yy <- stats::rnorm(m, c(0, 4)[s])
  d2 <- data.frame(t = seq_len(m), y = yy)
  term <- term_build(regime(2, time = t), d2)

  nll <- function(v) {
    -sum(term_loglik(term, rep(0, m), yy, gauss_ld(yy), gauss_sc(yy),
                     list(level1 = v[1], gap2 = exp(v[2]),
                          alr1.1 = v[3], alr2.1 = v[4]))$loglik)
  }
  fit <- stats::optim(c(0.5, log(2), 2, -2), nll, method = "BFGS")
  lev <- c(fit$par[1], fit$par[1] + exp(fit$par[2]))
  expect_equal(lev, c(0, 4), tolerance = 0.3)
  # and the estimated chain is persistent, as the one that generated it is
  Phat <- parameters7::param_value(parameters7::transition_matrix(2),
                                   fit$par[3:4])
  expect_gt(Phat[1, 1], 0.8)
  expect_gt(Phat[2, 2], 0.8)
})

test_that("groups run their own recursion", {
  term <- term_build(regime(2, by = g, time = t), dd)
  out <- term_loglik(term, rep(0, n), dd$y, gauss_ld(dd$y), gauss_sc(dd$y),
                     psi2)
  # a group's contributions depend only on its own rows
  d2 <- dd
  d2$y[dd$g == "b"] <- d2$y[dd$g == "b"] + 10
  out2 <- term_loglik(term, rep(0, n), d2$y, gauss_ld(d2$y), gauss_sc(d2$y),
                      psi2)
  expect_equal(out$loglik[dd$g == "a"], out2$loglik[dd$g == "a"])
  expect_false(isTRUE(all.equal(out$loglik[dd$g == "b"],
                                out2$loglik[dd$g == "b"])))
})

test_that("the degenerate cases are rejected and the term is routed", {
  term <- term_build(regime(2, time = t), dd)
  expect_error(
    term_loglik(term, rep(0, n), dd$y, gauss_ld(dd$y), gauss_sc(dd$y),
                utils::modifyList(psi2, list(gap2 = -1))),
    "ordered by construction")
  expect_error(term_loglik(regime(2), rep(0, n), dd$y, gauss_ld(dd$y),
                           gauss_sc(dd$y), psi2), "not been built")

  out <- interpret_formula(y ~ t + regime(2, time = t), dd)
  expect_named(out$terms, c("linpar", "regime(2, time = t)"))
  spec <- out$terms[["regime(2, time = t)"]]
  expect_true(S7::S7_inherits(spec, structural_term))
  expect_output(print(spec), "2 regimes")
  # the structural branch's other shape is refused, not faked
  expect_error(term_filter(term_build(spec, dd), rep(0, n), dd$y,
                           gauss_sc(dd$y), gauss_sc(dd$y), psi2),
               "does not implement term_filter")
})
