# term_shift_jacobian(): the posterior-weighted shift of a term of the
# likelihood shape and its derivative in the caller's unknowns, against the
# shift that term_posterior() and term_levels() compute, differenced by
# numDeriv over the same unknowns.

shift_check <- function(tm, dd, psi0) {
  nm <- term_params(tm)
  links <- term_links(tm)
  z0 <- vapply(nm, function(j)
    linkfunctions7::linkfun(links[[j]], psi0[[j]]), numeric(1))
  n <- nrow(dd)
  X <- cbind(1, dd$x2)
  Z <- cbind(1, dd$w)
  u0 <- unname(c(0.8, 0.3, log(0.5), 0.2, z0))
  m <- length(u0)
  cols <- 4L + seq_along(nm)
  fun <- function(u) {
    u <- unname(u)
    em <- as.numeric(X %*% u[1:2])
    es <- as.numeric(Z %*% u[3:4])
    psi <- stats::setNames(lapply(seq_along(nm), function(j)
      linkfunctions7::linkinv(links[[nm[j]]], u[4L + j])), nm)
    ld <- function(e, i) dnorm(dd$y[i], e, exp(es[i]), log = TRUE)
    g <- term_posterior(tm, em, dd$y, ld, psi)
    mu <- term_levels(tm, psi, eta = em, y = dd$y, logdens = ld)
    if (is.matrix(mu)) rowSums(g * mu) else as.numeric(g %*% mu)
  }
  em <- as.numeric(X %*% u0[1:2])
  es <- as.numeric(Z %*% u0[3:4])
  ld <- function(e, i) dnorm(dd$y[i], e, exp(es[i]), log = TRUE)
  gr <- function(e, i) {
    s <- exp(es[i])
    cbind((dd$y[i] - e) / s^2, ((dd$y[i] - e) / s)^2 - 1)
  }
  seed <- list(cbind(X, matrix(0, n, m - 2L)),
               cbind(matrix(0, n, 2L), Z, matrix(0, n, m - 4L)))
  out <- term_shift_jacobian(tm, em, dd$y, ld, gr,
                             stats::setNames(as.list(psi0[nm]), nm),
                             seed, cols, 1L)
  list(out = out, shift = fun(u0), jac = numDeriv::jacobian(fun, u0))
}

marg_panel <- function(seed = 7) {
  set.seed(seed)
  id <- rep(1:5, each = 10)
  x <- as.numeric(replicate(5, sort(runif(10, 0, 10))))
  pp <- 5 + rnorm(5, 0, 0.5)
  dd <- data.frame(id = id, x = x, x2 = rnorm(50), w = rnorm(50))
  dd$y <- 1 + 2 * (x >= pp[id]) + 0.8 * pmax(x - pp[id], 0) +
    rnorm(50, 0, 0.4)
  dd
}

test_that("a regime term's posterior shift is differentiated exactly", {
  skip_if_not_installed("numDeriv")
  set.seed(1)
  dd <- data.frame(t = 1:60, x2 = rnorm(60), w = rnorm(60))
  dd$y <- c(rnorm(20, 0), rnorm(20, 2), rnorm(20, 4.5)) + 0.3 * dd$x2
  r <- shift_check(term_build(regime(3, time = t), dd), dd,
                   list(level1 = 0, gap2 = 2, gap3 = 2.2, alr1.1 = 2,
                        alr1.2 = -1, alr2.1 = -1.5, alr2.2 = 1.8,
                        alr3.1 = -2, alr3.2 = 0.5))
  expect_equal(r$out$shift, r$shift, tolerance = 1e-12)
  expect_lt(max(abs(r$out$jacobian - r$jac)) / max(abs(r$jac)), 1e-8)
})

test_that("a marginal step's posterior shift is differentiated exactly", {
  skip_if_not_installed("numDeriv")
  dd <- marg_panel()
  r <- shift_check(term_build(jump(x, psi ~ random(~1 | id),
                                   marginal = TRUE), dd), dd,
                   list(psi1.mean = 4.8, psi1.sd = 0.6, delta1 = 1.8))
  expect_equal(r$out$shift, r$shift, tolerance = 1e-12)
  expect_lt(max(abs(r$out$jacobian - r$jac)) / max(abs(r$jac)), 1e-8)
  pr <- distributions7::fixed(distributions7::student_t1_distrib(), mu = 0)
  r <- shift_check(term_build(jump(x, psi ~ random(~1 | id, distrib = pr),
                                   marginal = TRUE), dd), dd,
                   list(psi1.mean = 4.8, sigma = 0.6, nu = 5, delta1 = 1.8))
  expect_lt(max(abs(r$out$jacobian - r$jac)) / max(abs(r$jac)), 1e-8)
})

test_that("two marginal steps ride the side chain of four patterns", {
  skip_if_not_installed("numDeriv")
  set.seed(21)
  id <- rep(1:4, each = 12)
  x <- as.numeric(replicate(4, sort(runif(12, 0, 10))))
  p1 <- 3 + rnorm(4, 0, 0.4)
  p2 <- 7 + rnorm(4, 0, 0.4)
  dd <- data.frame(id = id, x = x, x2 = rnorm(48), w = rnorm(48))
  dd$y <- 1 + 2 * (x >= p1[id]) - 1.5 * (x >= p2[id]) + rnorm(48, 0, 0.4)
  r <- shift_check(term_build(jump(x, psi ~ random(~1 | id), npsi = 2,
                                   marginal = TRUE), dd), dd,
                   list(psi1.mean = 3, psi1.sd = 0.5, psi2.mean = 7,
                        psi2.sd = 0.5, delta1 = 2, delta2 = -1.5))
  expect_lt(max(abs(r$out$jacobian - r$jac)) / max(abs(r$jac)), 1e-8)
})

test_that("the continuous marginal kinds move the posterior over the nodes", {
  skip_if_not_installed("numDeriv")
  dd <- marg_panel()
  r <- shift_check(term_build(seg(x, psi ~ random(~1 | id),
                                  marginal = TRUE), dd), dd,
                   list(psi1.mean = 4.8, psi1.sd = 0.7, beta = 0.1,
                        gamma1 = 0.7))
  expect_equal(r$out$shift, r$shift, tolerance = 1e-12)
  expect_lt(max(abs(r$out$jacobian - r$jac)) / max(abs(r$jac)), 1e-7)
  r <- shift_check(term_build(jseg(x, psi ~ random(~1 | id),
                                   marginal = TRUE), dd), dd,
                   list(psi1.mean = 4.8, psi1.sd = 0.7, beta = 0.1,
                        gamma1 = 0.7, delta1 = 1.5))
  expect_lt(max(abs(r$out$jacobian - r$jac)) / max(abs(r$jac)), 1e-7)
})

test_that("a term of another shape and a malformed seed are rejected", {
  set.seed(1)
  dd <- data.frame(t = 1:20, y = rnorm(20))
  tm <- term_build(regime(2, time = t), dd)
  psi <- list(level1 = 0, gap2 = 1, alr1.1 = 1, alr2.1 = -1)
  ld <- function(e, i) dnorm(dd$y[i], e, log = TRUE)
  gr <- function(e, i) matrix(dd$y[i] - e, ncol = 1L)
  expect_error(term_shift_jacobian(tm, rep(0, 20), dd$y, ld, gr, psi,
                                   list(matrix(0, 19, 5)), 2:5, 1L),
               "must be 20 by 5")
  expect_error(term_shift_jacobian(tm, rep(0, 20), dd$y, ld, gr, psi,
                                   list(matrix(0, 20, 5)), 2:4, 1L),
               "'cols' must give the 4 columns")
  g <- term_build(gas(1, 1, time = t), dd)
  expect_error(term_shift_jacobian(g, rep(0, 20), dd$y, ld, gr, list(),
                                   list(matrix(0, 20, 5)), 2:5, 1L),
               "does not implement term_shift_jacobian")
})
