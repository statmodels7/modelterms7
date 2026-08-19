# The marginal break-point term: jump(x, psi ~ random(~1 | g),
# marginal = TRUE). The references are chosen so that no check shares its
# arithmetic with the implementation: the total against a brute-force
# quadrature of the marginal on a fine grid, the jacobian and the Hessian
# against numDeriv, and Fisher's identity against the posterior the term
# reports.

marg_data <- function(mI = 5L, nI = 10L, seed = 7) {
  set.seed(seed)
  id <- rep(seq_len(mI), each = nI)
  x <- as.numeric(replicate(mI, sort(runif(nI, 0, 10))))
  psi <- 5 + rnorm(mI, 0, 0.5)
  y <- 1 + 2 * (x >= psi[id]) + rnorm(mI * nI, 0, 0.4)
  data.frame(id = id, x = x, y = y, psi_true = psi[id])
}

marg_cb <- function(dd, s0 = 0.4) {
  list(ld = function(e, i) dnorm(dd$y[i], e, s0, log = TRUE),
       sc = function(e, i) (dd$y[i] - e) / s0^2)
}

test_that("the constructor takes the psi ~ random subformula and nothing else", {
  tm <- jump(x, psi ~ random(~1 | id), marginal = TRUE)
  expect_true(S7::S7_inherits(tm, MarginalBreakTerm))
  expect_true(S7::S7_inherits(tm, structural_term))
  expect_identical(term_params(tm), c("m1", "tau1", "delta1"))

  expect_error(jump(x, marginal = TRUE), "requires the break-point")
  expect_error(jump(x, psi ~ random(~1 | id), npsi = 4, marginal = TRUE),
               "up to three")
  expect_error(seg(x, psi ~ random(~1 | id), npsi = 2, marginal = TRUE),
               "one break-point")
  expect_error(jump(x, delta ~ id, psi ~ random(~1 | id), marginal = TRUE),
               "cannot be developed")
  expect_error(jump(x, psi ~ id, marginal = TRUE),
               "single random")
  expect_error(jump(x, psi ~ random(~ x | id), marginal = TRUE),
               "intercept-only")
  expect_error(jump(x, psi ~ random(~1 | id, distrib =
                                      distributions7::student_t1_distrib()),
                    marginal = TRUE),
               "fixed at zero")
  expect_error(jump(x, psi ~ random(~1 | id, distrib =
                                      distributions7::fixed(
                                        distributions7::student_t1_distrib(),
                                        mu = 0)),
                    npsi = 2, marginal = TRUE),
               "gaussian")
  expect_error(jump(x, psi ~ random(~1 | id), by = ~0 + id, marginal = TRUE),
               "does not combine")
  expect_error(seg(x, marginal = TRUE), "requires the break-point")
  expect_error(jseg(x, marginal = TRUE), "requires the break-point")

  expect_message(jump(x, psi ~ random(~1 | id), marginal = TRUE,
                      smoothed = penalties7::smooth_probit()),
                 "ignored")
  expect_message(jump(x, psi ~ random(~1 | id), marginal = TRUE, c0 = 0.1),
                 "ignored")
  expect_message(jump(x, psi ~ random(~1 | id), marginal = TRUE, n_boot = 5),
                 "exact profile")
})

test_that("the default marginal = FALSE leaves the sharp construction alone", {
  dd <- marg_data()
  a <- term_build(jump(x), dd)
  b <- term_build(jump(x, marginal = FALSE), dd)
  expect_identical(term_matrix(a), term_matrix(b))
  expect_identical(term_coef_names(a), term_coef_names(b))
})

test_that("the exact sum agrees with a brute-force quadrature of the marginal", {
  dd <- marg_data(4L, 8L)
  tm <- term_build(jump(x, psi ~ random(~1 | id), marginal = TRUE), dd)
  cb <- marg_cb(dd)
  eta <- rep(1, nrow(dd))
  psi <- list(m1 = 4.8, tau1 = 0.6, delta1 = 1.8)
  out <- term_loglik(tm, eta, dd$y, cb$ld, cb$sc, psi)

  brute <- 0
  for (g in unique(dd$id)) {
    r <- which(dd$id == g)
    ps <- seq(4.8 - 10 * 0.6, 4.8 + 10 * 0.6, length.out = 200001)
    lc <- vapply(ps, function(p)
      sum(dnorm(dd$y[r], eta[r] + 1.8 * (dd$x[r] >= p), 0.4, log = TRUE)),
      numeric(1))
    mx <- max(lc)
    brute <- brute +
      log(sum(exp(lc - mx) * dnorm(ps, 4.8, 0.6)) * (ps[2] - ps[1])) + mx
  }
  # the reference's own error is the grid's, first order at the jumps of
  # the conditional in psi: at this step it sits below 3e-6 relative
  expect_equal(sum(out$loglik), brute, tolerance = 3e-6)
})

test_that("the jacobian is exact against numDeriv, row by row", {
  skip_if_not_installed("numDeriv")
  dd <- marg_data(4L, 8L)
  tm <- term_build(jump(x, psi ~ random(~1 | id), marginal = TRUE), dd)
  cb <- marg_cb(dd)
  eta <- rep(1, nrow(dd))
  v0 <- c(4.8, 0.6, 1.8)
  out <- term_loglik(tm, eta, dd$y, cb$ld, cb$sc,
                     list(m1 = v0[1], tau1 = v0[2], delta1 = v0[3]))
  J <- numDeriv::jacobian(function(v)
    term_loglik(tm, eta, dd$y, cb$ld, cb$sc,
                list(m1 = v[1], tau1 = v[2], delta1 = v[3]))$loglik, v0)
  expect_lt(max(abs(J - out$jacobian)), 1e-6 * max(1, max(abs(J))))
})

test_that("the posterior carries Fisher's identity and the levels", {
  skip_if_not_installed("numDeriv")
  dd <- marg_data(4L, 8L)
  tm <- term_build(jump(x, psi ~ random(~1 | id), marginal = TRUE), dd)
  cb <- marg_cb(dd)
  eta <- rep(1, nrow(dd))
  psi <- list(m1 = 4.8, tau1 = 0.6, delta1 = 1.8)
  P <- term_posterior(tm, eta, dd$y, cb$ld, psi)
  expect_equal(rowSums(P), rep(1, nrow(dd)), tolerance = 1e-12)
  # the shifted probability grows with the covariate within a group, the
  # posterior over the intervals being cumulated along it
  for (g in unique(dd$id)) {
    r <- which(dd$id == g)
    expect_true(all(diff(P[r, 2L][order(dd$x[r])]) >= -1e-12))
  }
  expect_identical(term_levels(tm, psi), c(0, 1.8))

  # Fisher's identity in one predictor, against a difference of the total
  out <- term_loglik(tm, eta, dd$y, cb$ld, cb$sc, psi)
  s1 <- cb$sc(eta + psi$delta1, seq_len(nrow(dd)))
  expect_equal(sum(P[, 2L] * s1), sum(out$jacobian[, "delta1"]),
               tolerance = 1e-10)
  tt <- 7L
  g_num <- numDeriv::grad(function(e) {
    et <- eta
    et[tt] <- e
    sum(term_loglik(tm, et, dd$y, cb$ld, cb$sc, psi)$loglik)
  }, eta[tt])
  g_fi <- (1 - P[tt, 2L]) * cb$sc(eta[tt], tt) + P[tt, 2L] * s1[tt]
  expect_equal(g_fi, g_num, tolerance = 1e-7)
})

test_that("term_levels answers for a regime term too", {
  ps <- list(level1 = 0.5, gap2 = 3, alr1.1 = 2, alr2.1 = -2)
  expect_identical(term_levels(regime(2), ps), c(0.5, 3.5))
})

test_that("the observed Hessian agrees with numDeriv on the joint unknowns", {
  skip_if_not_installed("numDeriv")
  dd <- marg_data(4L, 8L)
  n <- nrow(dd)
  tm <- term_build(jump(x, psi ~ random(~1 | id), marginal = TRUE), dd)
  cb <- marg_cb(dd)
  eta <- rep(1, n)
  psi <- list(m1 = 4.8, tau1 = 0.6, delta1 = 1.8)
  seed <- list(cbind(rep(1, n), 0, 0, 0))
  oh <- term_hessian(tm, eta, dd$y, cb$ld,
                     grad = function(e, i)
                       matrix((dd$y[i] - e) / 0.4^2, ncol = 1L),
                     hess = function(e, i)
                       array(-1 / 0.4^2, c(length(i), 1L, 1L)),
                     psi = psi, seed = seed, cols = 2:4, level = 1L)
  f_u <- function(u) {
    sum(term_loglik(tm, rep(u[1], n), dd$y, cb$ld, cb$sc,
                    list(m1 = u[2], tau1 = exp(u[3]), delta1 = u[4]))$loglik)
  }
  u0 <- c(1, 4.8, log(0.6), 1.8)
  expect_equal(oh$gradient, numDeriv::grad(f_u, u0), tolerance = 1e-6)
  Hn <- numDeriv::hessian(f_u, u0)
  expect_lt(max(abs(Hn - oh$hessian)), 1e-5 * max(abs(Hn)))
  # the contributions are the same predictive decomposition term_loglik
  # returns
  out <- term_loglik(tm, eta, dd$y, cb$ld, cb$sc, psi)
  expect_equal(oh$loglik, out$loglik, tolerance = 1e-12)

  # a weight varying within a group has no reading on a marginal likelihood
  w <- rep(1, n)
  w[2L] <- 2
  expect_error(term_hessian(tm, eta, dd$y, cb$ld,
                            grad = function(e, i)
                              matrix((dd$y[i] - e) / 0.4^2, ncol = 1L),
                            hess = function(e, i)
                              array(-1 / 0.4^2, c(length(i), 1L, 1L)),
                            psi = psi, seed = seed, cols = 2:4, level = 1L,
                            weights = w),
               "constant within each group")
})

test_that("the latent posterior tracks the truth at the generating values", {
  dd <- marg_data(8L, 16L, seed = 11)
  tm <- term_build(jump(x, psi ~ random(~1 | id), marginal = TRUE), dd)
  cb <- marg_cb(dd)
  lat <- term_latent(tm, rep(1, nrow(dd)), dd$y, cb$ld,
                     list(m1 = 5, tau1 = 0.5, delta1 = 2))
  expect_identical(nrow(lat), 8L)
  expect_true(all(is.finite(lat$mean)) && all(lat$sd > 0))
  truth <- tapply(dd$psi_true, dd$id, function(v) v[1L])
  expect_gt(cor(lat$mean, as.numeric(truth)), 0.6)
})

test_that("the start reads the exact profile off the target", {
  dd <- marg_data(8L, 16L, seed = 13)
  tm <- term_build(jump(x, psi ~ random(~1 | id), marginal = TRUE), dd)
  st <- term_start(tm, target = dd$y)
  expect_identical(names(st), c("m1", "tau1", "delta1"))
  expect_lt(abs(st[["m1"]] - 5), 1)
  expect_gt(st[["delta1"]], 1)
  expect_true(is.finite(st[["tau1"]]))
  # without a target the covariate answers, with the change at zero
  st0 <- term_start(tm)
  expect_identical(st0[["delta1"]], 0)
  # a caller's psi seeds the position either way
  tm2 <- term_build(jump(x, psi ~ random(~1 | id), marginal = TRUE,
                         psi = 3.3), dd)
  expect_identical(term_start(tm2)[["m1"]], 3.3)
  expect_identical(term_start(tm2, target = dd$y)[["m1"]], 3.3)
})

# ---- several break-points (the product partition) ----

marg_data2 <- function(mI = 4L, nI = 12L, seed = 21) {
  set.seed(seed)
  id <- rep(seq_len(mI), each = nI)
  x <- as.numeric(replicate(mI, sort(runif(nI, 0, 10))))
  p1 <- 3 + rnorm(mI, 0, 0.4)
  p2 <- 7 + rnorm(mI, 0, 0.4)
  y <- 1 + 2 * (x >= p1[id]) - 1.5 * (x >= p2[id]) + rnorm(mI * nI, 0, 0.4)
  data.frame(id = id, x = x, y = y)
}

test_that("two latent break-points agree with a brute-force 2-D quadrature", {
  dd <- marg_data2(3L, 7L)
  tm <- term_build(jump(x, psi ~ random(~1 | id), npsi = 2, marginal = TRUE),
                   dd)
  expect_identical(term_params(tm),
                   c("m1", "tau1", "m2", "tau2", "delta1", "delta2"))
  cb <- marg_cb(dd)
  eta <- rep(1, nrow(dd))
  psi <- list(m1 = 3.1, tau1 = 0.5, m2 = 6.8, tau2 = 0.4,
              delta1 = 1.8, delta2 = -1.2)
  out <- term_loglik(tm, eta, dd$y, cb$ld, cb$sc, psi)

  # the independent route: the same exact sum written as plain loops, the
  # conditional evaluated from scratch at a representative point of every
  # cell of the product partition and the masses as direct cdf differences
  brute <- 0
  for (g in unique(dd$id)) {
    r <- which(dd$id == g)
    xs <- sort(dd$x[r])
    b <- c(-Inf, xs, Inf)
    parts <- numeric(0)
    for (j1 in seq_len(length(b) - 1L)) {
      for (j2 in seq_len(length(b) - 1L)) {
        lm1 <- log(pnorm((b[j1 + 1] - 3.1) / 0.5) -
                     pnorm((b[j1] - 3.1) / 0.5))
        lm2 <- log(pnorm((b[j2 + 1] - 6.8) / 0.4) -
                     pnorm((b[j2] - 6.8) / 0.4))
        if (!is.finite(lm1) || !is.finite(lm2)) next
        p1 <- if (j1 == 1L) xs[1] - 1 else if (j1 == length(b) - 1L)
          xs[length(xs)] + 1 else (b[j1] + b[j1 + 1]) / 2
        p2 <- if (j2 == 1L) xs[1] - 1 else if (j2 == length(b) - 1L)
          xs[length(xs)] + 1 else (b[j2] + b[j2 + 1]) / 2
        lc <- sum(dnorm(dd$y[r],
                        eta[r] + 1.8 * (dd$x[r] >= p1) -
                          1.2 * (dd$x[r] >= p2), 0.4, log = TRUE))
        parts <- c(parts, lm1 + lm2 + lc)
      }
    }
    mx <- max(parts)
    brute <- brute + mx + log(sum(exp(parts - mx)))
  }
  expect_equal(sum(out$loglik), brute, tolerance = 1e-8)

  skip_if_not_installed("numDeriv")
  v0 <- c(3.1, 0.5, 6.8, 0.4, 1.8, -1.2)
  J <- numDeriv::jacobian(function(v)
    term_loglik(tm, eta, dd$y, cb$ld, cb$sc,
                list(m1 = v[1], tau1 = v[2], m2 = v[3], tau2 = v[4],
                     delta1 = v[5], delta2 = v[6]))$loglik, v0)
  expect_lt(max(abs(J - out$jacobian)), 1e-6 * max(1, max(abs(J))))

  # the pattern posterior: rows sum to one, and the levels are the pattern
  # sums of the changes
  P <- term_posterior(tm, eta, dd$y, cb$ld, psi)
  expect_identical(ncol(P), 4L)
  expect_equal(rowSums(P), rep(1, nrow(dd)), tolerance = 1e-10)
  expect_equal(term_levels(tm, psi), c(0, 1.8, -1.2, 0.6))
})

test_that("the two-break-point Hessian agrees with numDeriv", {
  skip_if_not_installed("numDeriv")
  dd <- marg_data2(3L, 6L)
  n <- nrow(dd)
  tm <- term_build(jump(x, psi ~ random(~1 | id), npsi = 2, marginal = TRUE),
                   dd)
  cb <- marg_cb(dd)
  eta <- rep(1, n)
  psi <- list(m1 = 3.1, tau1 = 0.5, m2 = 6.8, tau2 = 0.4,
              delta1 = 1.8, delta2 = -1.2)
  seed <- list(cbind(rep(1, n), matrix(0, n, 6L)))
  oh <- term_hessian(tm, eta, dd$y, cb$ld,
                     grad = function(e, i)
                       matrix((dd$y[i] - e) / 0.4^2, ncol = 1L),
                     hess = function(e, i)
                       array(-1 / 0.4^2, c(length(i), 1L, 1L)),
                     psi = psi, seed = seed, cols = 2:7, level = 1L)
  f_u <- function(u) {
    sum(term_loglik(tm, rep(u[1], n), dd$y, cb$ld, cb$sc,
                    list(m1 = u[2], tau1 = exp(u[3]), m2 = u[4],
                         tau2 = exp(u[5]), delta1 = u[6],
                         delta2 = u[7]))$loglik)
  }
  u0 <- c(1, 3.1, log(0.5), 6.8, log(0.4), 1.8, -1.2)
  expect_equal(oh$gradient, numDeriv::grad(f_u, u0), tolerance = 1e-6)
  Hn <- numDeriv::hessian(f_u, u0)
  expect_lt(max(abs(Hn - oh$hessian)), 1e-4 * max(abs(Hn)))
})

test_that("the general step Hessian matches the propagated one at one break-point", {
  dd <- marg_data(4L, 8L)
  n <- nrow(dd)
  tm <- term_build(jump(x, psi ~ random(~1 | id), marginal = TRUE), dd)
  cb <- marg_cb(dd)
  eta <- rep(1, n)
  psi <- list(m1 = 4.8, tau1 = 0.6, delta1 = 1.8)
  seed <- list(cbind(rep(1, n), 0, 0, 0))
  gr <- function(e, i) matrix((dd$y[i] - e) / 0.4^2, ncol = 1L)
  he <- function(e, i) array(-1 / 0.4^2, c(length(i), 1L, 1L))
  a <- term_hessian(tm, eta, dd$y, cb$ld, gr, he, psi, seed = seed,
                    cols = 2:4, level = 1L)
  b <- .marg_jump_hessian(tm, eta, dd$y, cb$ld, gr, he, psi,
                          seed = lapply(seed, as.matrix), cols = 2:4,
                          level = 1L, w = rep(1, n))
  expect_equal(a$gradient, b$gradient, tolerance = 1e-10)
  expect_equal(a$hessian, b$hessian, tolerance = 1e-8)
})

# ---- the continuous kinds (quadrature per interval) ----

marg_seg_data <- function(mI = 4L, nI = 12L, seed = 31, kind = "seg") {
  set.seed(seed)
  id <- rep(seq_len(mI), each = nI)
  x <- as.numeric(replicate(mI, sort(runif(nI, 0, 10))))
  psi <- 5 + rnorm(mI, 0, 0.5)
  y <- 1 + 0.5 * x - 1.2 * pmax(x - psi[id], 0) +
    (if (kind == "jseg") 1.5 * (x >= psi[id]) else 0) +
    rnorm(mI * nI, 0, 0.4)
  data.frame(id = id, x = x, y = y, psi_true = psi[id])
}

test_that("the seg marginal agrees with a fine quadrature and numDeriv", {
  dd <- marg_seg_data(4L, 10L)
  tm <- term_build(seg(x, psi ~ random(~1 | id), marginal = TRUE), dd)
  expect_identical(term_params(tm), c("beta", "m1", "tau1", "gamma1"))
  cb <- marg_cb(dd)
  eta <- rep(1, nrow(dd))
  psi <- list(beta = 0.5, m1 = 4.8, tau1 = 0.6, gamma1 = -1.2)
  out <- term_loglik(tm, eta, dd$y, cb$ld, cb$sc, psi)

  # a fine trapezoid over the prior's support, sharing nothing with the
  # panels
  brute <- 0
  ps <- seq(4.8 - 9 * 0.6, 4.8 + 9 * 0.6, length.out = 40001)
  for (g in unique(dd$id)) {
    r <- which(dd$id == g)
    lc <- vapply(ps, function(p)
      sum(dnorm(dd$y[r],
                eta[r] + 0.5 * dd$x[r] - 1.2 * pmax(dd$x[r] - p, 0),
                0.4, log = TRUE)), numeric(1))
    lw <- lc + dnorm(ps, 4.8, 0.6, log = TRUE)
    mx <- max(lw)
    brute <- brute + mx + log(sum(exp(lw - mx)) * (ps[2] - ps[1]))
  }
  expect_equal(sum(out$loglik), brute, tolerance = 1e-6)

  skip_if_not_installed("numDeriv")
  v0 <- c(0.5, 4.8, 0.6, -1.2)
  J <- numDeriv::jacobian(function(v)
    term_loglik(tm, eta, dd$y, cb$ld, cb$sc,
                list(beta = v[1], m1 = v[2], tau1 = v[3],
                     gamma1 = v[4]))$loglik, v0)
  expect_lt(max(abs(J - out$jacobian)), 2e-5 * max(1, max(abs(J))))

  # the posterior over the nodes: rows sum to one and the latent summary
  # tracks the truth at the generating values
  P <- term_posterior(tm, eta, dd$y, cb$ld, psi)
  expect_equal(rowSums(P), rep(1, nrow(dd)), tolerance = 1e-10)
  lat <- term_latent(tm, eta, dd$y, cb$ld,
                     list(beta = 0.5, m1 = 5, tau1 = 0.5, gamma1 = -1.2))
  truth <- tapply(dd$psi_true, dd$id, function(v) v[1L])
  expect_gt(cor(lat$mean, as.numeric(truth)), 0.5)
})

test_that("the seg marginal jacobian carries the edge panels' node motion", {
  skip_if_not_installed("numDeriv")
  # m within a couple of prior sds of the data's lower end, where the lower
  # region carries real mass and a fixed-node gradient would be visibly
  # wrong
  set.seed(41)
  id <- rep(1:3, each = 8)
  x <- as.numeric(replicate(3, sort(runif(8, 0, 10))))
  y <- 1 - 0.8 * pmax(x - 1.6, 0) + rnorm(24, 0, 0.4)
  dd <- data.frame(id = id, x = x, y = y)
  tm <- term_build(seg(x, psi ~ random(~1 | id), marginal = TRUE,
                       linear = FALSE), dd)
  cb <- marg_cb(dd)
  eta <- rep(1, nrow(dd))
  v0 <- c(1.2, 0.7, -0.8)
  out <- term_loglik(tm, eta, dd$y, cb$ld, cb$sc,
                     list(m1 = v0[1], tau1 = v0[2], gamma1 = v0[3]))
  J <- numDeriv::jacobian(function(v)
    term_loglik(tm, eta, dd$y, cb$ld, cb$sc,
                list(m1 = v[1], tau1 = v[2], gamma1 = v[3]))$loglik, v0)
  expect_lt(max(abs(J - out$jacobian)), 2e-5 * max(1, max(abs(J))))
})

test_that("the seg and jseg Hessians agree with numDeriv", {
  skip_if_not_installed("numDeriv")
  for (kind in c("seg", "jseg")) {
    dd <- marg_seg_data(3L, 8L, seed = 51, kind = kind)
    n <- nrow(dd)
    ctor <- if (kind == "seg") seg else jseg
    tm <- term_build(ctor(x, psi ~ random(~1 | id), marginal = TRUE), dd)
    cb <- marg_cb(dd)
    eta <- rep(1, n)
    nmv <- term_params(tm)
    v0 <- c(beta = 0.5, m1 = 4.8, tau1 = 0.6, gamma1 = -1.2,
            delta1 = 1.5)[nmv]
    np <- length(nmv)
    seed <- list(cbind(rep(1, n), matrix(0, n, np)))
    oh <- term_hessian(tm, eta, dd$y, cb$ld,
                       grad = function(e, i)
                         matrix((dd$y[i] - e) / 0.4^2, ncol = 1L),
                       hess = function(e, i)
                         array(-1 / 0.4^2, c(length(i), 1L, 1L)),
                       psi = as.list(v0), seed = seed,
                       cols = 1L + seq_len(np), level = 1L)
    f_u <- function(u) {
      vv <- as.list(stats::setNames(u[-1L], nmv))
      vv$tau1 <- exp(vv$tau1)
      sum(term_loglik(tm, rep(u[1], n), dd$y, cb$ld, cb$sc, vv)$loglik)
    }
    z0 <- v0
    z0["tau1"] <- log(v0[["tau1"]])
    u0 <- c(1, unname(z0))
    expect_equal(oh$gradient, numDeriv::grad(f_u, u0), tolerance = 1e-5)
    Hn <- numDeriv::hessian(f_u, u0)
    expect_lt(max(abs(Hn - oh$hessian)), 5e-4 * max(abs(Hn)))
  }
})

# ---- an explicit prior (the cdf surface route) ----

test_that("a t prior's masses and derivatives ride the cdf surface", {
  skip_if_not_installed("numDeriv")
  dd <- marg_data(4L, 8L, seed = 61)
  pr <- distributions7::fixed(distributions7::student_t1_distrib(), mu = 0)
  tm <- term_build(jump(x, psi ~ random(~1 | id, distrib = pr),
                        marginal = TRUE), dd)
  expect_identical(term_params(tm), c("m1", "sigma", "nu", "delta1"))
  cb <- marg_cb(dd)
  eta <- rep(1, nrow(dd))
  psi <- list(m1 = 4.8, sigma = 0.6, nu = 5, delta1 = 1.8)
  out <- term_loglik(tm, eta, dd$y, cb$ld, cb$sc, psi)

  # the independent route: a fine quadrature of the same marginal with the
  # t prior's own density
  brute <- 0
  ps <- seq(4.8 - 40 * 0.6, 4.8 + 40 * 0.6, length.out = 400001)
  for (g in unique(dd$id)) {
    r <- which(dd$id == g)
    lc <- vapply(ps, function(p)
      sum(dnorm(dd$y[r], eta[r] + 1.8 * (dd$x[r] >= p), 0.4, log = TRUE)),
      numeric(1))
    lw <- lc + stats::dt((ps - 4.8) / 0.6, df = 5, log = TRUE) - log(0.6)
    mx <- max(lw)
    brute <- brute + mx + log(sum(exp(lw - mx)) * (ps[2] - ps[1]))
  }
  expect_equal(sum(out$loglik), brute, tolerance = 1e-5)

  # the jacobian in (m, sigma, delta) is closed through the surface; nu's
  # column carries the surface's own finite difference and is compared at
  # the tolerance that route documents
  v0 <- c(4.8, 0.6, 5, 1.8)
  J <- numDeriv::jacobian(function(v)
    term_loglik(tm, eta, dd$y, cb$ld, cb$sc,
                list(m1 = v[1], sigma = v[2], nu = v[3],
                     delta1 = v[4]))$loglik, v0)
  expect_lt(max(abs(J[, -3L] - out$jacobian[, -3L])),
            1e-6 * max(1, max(abs(J))))
  expect_lt(max(abs(J[, 3L] - out$jacobian[, 3L])), 1e-4)

  # the Hessian: the prior's rows ride one stencil on the analytic gradient
  n <- nrow(dd)
  seed <- list(cbind(rep(1, n), matrix(0, n, 4L)))
  oh <- term_hessian(tm, eta, dd$y, cb$ld,
                     grad = function(e, i)
                       matrix((dd$y[i] - e) / 0.4^2, ncol = 1L),
                     hess = function(e, i)
                       array(-1 / 0.4^2, c(length(i), 1L, 1L)),
                     psi = psi, seed = seed, cols = 2:5, level = 1L)
  lk_nu <- term_links(tm)$nu
  f_u <- function(u) {
    sum(term_loglik(tm, rep(u[1], n), dd$y, cb$ld, cb$sc,
                    list(m1 = u[2], sigma = exp(u[3]),
                         nu = linkfunctions7::linkinv(lk_nu, u[4]),
                         delta1 = u[5]))$loglik)
  }
  u0 <- c(1, 4.8, log(0.6), linkfunctions7::linkfun(lk_nu, 5), 1.8)
  Hn <- numDeriv::hessian(f_u, u0)
  expect_lt(max(abs(Hn - oh$hessian)), 5e-3 * max(abs(Hn)))
})

