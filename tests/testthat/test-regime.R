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

test_that("the compiled recursion agrees with the R twin", {
  # Same operations in the same order, to a rounding: the kernel
  # accumulates sums the R form accumulates too, and a compiler free to
  # contract a multiply-add moves the last bit.
  set.seed(51)
  for (k in c(2L, 3L)) {
    n <- 200L; np <- (k - 1L) + k * (k - 1L) + 1L
    LF <- matrix(rnorm(n * k, -1, 0.5), n, k)
    SC <- matrix(rnorm(n * k), n, k)
    dmu <- matrix(0, k, np); dmu[, 1L] <- 1
    for (j in seq_len(k)) if (j > 1L) dmu[j, 1L + seq_len(j - 1L)] <- 1
    ch <- parameters7::transition_matrix(k)
    fv <- stats::setNames(runif(length(ch@free_names), -0.5, 0.5),
                          ch@free_names)
    P <- unclass(parameters7::param_value(ch, fv))
    d1 <- parameters7::param_d1(ch, fv)
    dP <- lapply(seq_len(np), function(i) matrix(0, k, k))
    off <- k - 1L
    for (i in seq_along(ch@free_names)) dP[[off + i]] <- d1[[ch@free_names[i]]]
    st <- modelterms7:::regime_stationary(P, dP)
    ddm <- do.call(rbind, st$ddelta)
    ord <- list(seq_len(120L), seq.int(121L, n))

    a <- modelterms7:::.regime_forward_r(ord, LF, SC, dmu, P, dP,
                                         st$delta, ddm)
    b <- modelterms7:::regime_forward_cpp(ord, LF, SC, dmu, P, dP,
                                          st$delta, ddm)
    expect_equal(b$loglik, a$loglik, tolerance = 1e-13, info = paste("k", k))
    expect_equal(b$jacobian, a$jacobian, tolerance = 1e-13,
                 info = paste("k", k))
  }
  # no rows is a legal group after a subset
  e <- modelterms7:::regime_forward_cpp(list(integer(0)),
                                        matrix(0, 0, 2), matrix(0, 0, 2),
                                        matrix(1, 2, 2), diag(2),
                                        list(diag(2), diag(2)),
                                        c(0.5, 0.5), matrix(0, 2, 2))
  expect_identical(length(e$loglik), 0L)
})

test_that("the recursion equals the sum over every state path", {
  # A hidden Markov likelihood IS a sum over paths, and at this length it
  # can be taken in full: a reference the recursion shares no arithmetic
  # with, and the one that would catch a wrong normalization.
  set.seed(52)
  K <- 2L; Tn <- 11L
  dd <- data.frame(y = rnorm(Tn, sd = 1.5), t = seq_len(Tn))
  br <- term_build(regime(K, time = t), dd)
  psi <- list(level1 = -0.6, gap2 = 2.1, alr1.1 = 0.8, alr2.1 = -0.5)
  eta <- seq(0, 0.5, length.out = Tn)
  got <- term_loglik(br, eta, dd$y,
                     logdens = function(e, i) dnorm(dd$y[i], e, 1, log = TRUE),
                     score = function(e, i) dd$y[i] - e, psi = psi)

  mu <- cumsum(c(psi$level1, psi$gap2))
  ch <- parameters7::transition_matrix(K)
  P <- unclass(parameters7::param_value(ch, unlist(psi[ch@free_names])))
  ev <- eigen(t(P))
  p0 <- Re(ev$vectors[, which.max(Re(ev$values))])
  p0 <- p0 / sum(p0)
  paths <- as.matrix(expand.grid(rep(list(seq_len(K)), Tn)))
  lp <- apply(paths, 1, function(sq) {
    l <- log(p0[sq[1]]) + dnorm(dd$y[1], eta[1] + mu[sq[1]], 1, log = TRUE)
    for (t in 2:Tn) {
      l <- l + log(P[sq[t - 1], sq[t]]) +
        dnorm(dd$y[t], eta[t] + mu[sq[t]], 1, log = TRUE)
    }
    l
  })
  expect_equal(sum(got$loglik), log(sum(exp(lp - max(lp)))) + max(lp),
               tolerance = 1e-10)
})

test_that("the closures are called with the whole index vector", {
  set.seed(53)
  dd <- data.frame(y = rnorm(30), t = 1:30)
  br <- term_build(regime(2, time = t), dd)
  psi <- list(level1 = 0, gap2 = 1, alr1.1 = 0.2, alr2.1 = -0.2)
  calls <- 0L
  ld <- function(e, i) { calls <<- calls + 1L; dnorm(dd$y[i], e, log = TRUE) }
  term_loglik(br, rep(0, 30), dd$y, logdens = ld,
              score = function(e, i) dd$y[i] - e, psi = psi)
  # one call per regime, not one per observation and regime
  expect_identical(calls, 2L)

  # and a closure that ignores the index is refused where it happens
  expect_error(
    term_loglik(br, rep(0, 30), dd$y,
                logdens = function(e, i) dnorm(dd$y[1], e[1], log = TRUE),
                score = function(e, i) dd$y[i] - e, psi = psi),
    "one value per observation")
})
