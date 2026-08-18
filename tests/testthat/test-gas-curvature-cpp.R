# The compiled second-order recursion of the submodel route against the R
# route, which stays the reference: same term, same point, the callbacks as
# lookups and the layer's blocks as data on one side, the plain callbacks on
# the other. The comparison carries a tolerance (the seg_block rule: a
# compiler is free to contract a multiply-add), while the SAME kernel across
# thread counts is compared with identical() -- the groups are decomposed by
# output element and no reduction is split.

set.seed(21)
n_c <- 80
dd_c <- data.frame(t = seq_len(n_c), y = rnorm(n_c),
                   z = as.numeric(scale(runif(n_c))),
                   g = factor(rep(c("a", "b"), each = n_c / 2)),
                   g8 = factor(rep(letters[1:8], each = n_c / 8)))

gauss_score_c <- function(y) function(e, i) y[i] - e
gauss_curv_c <- function(y) function(e, i) -1

# the layer's pieces for a synthetic two-parameter model, the filter on the
# first: the callback mirrors statmodels7's .structural_blocks() and the
# data list is the same arrays, so the two routes read one model
.mk_layer <- function(n, m) {
  set.seed(31)
  H <- cbind(0, 0.3 * sin(seq_len(n) / 7))
  D3 <- cbind(0.2 * cos(seq_len(n) / 9), 0.1 * sin(seq_len(n) / 5),
              0.05 * cos(seq_len(n) / 4), -0.15 * sin(seq_len(n) / 11))
  Vs <- list(matrix(0, n, m),
             matrix(rnorm(n * m, sd = 0.4), n, m))
  calls <- new.env()
  calls$k <- 0L
  cb <- function(e, i, D, act) {
    calls$k <- calls$k + 1L
    mk <- length(act)
    cross <- H[i, 2L] * Vs[[2L]][i, act]
    vr <- list(D, Vs[[2L]][i, act])
    M <- matrix(0, mk, mk)
    for (r in 1:2) {
      for (r2 in 1:2) {
        M <- M + D3[i, (r - 1L) * 2L + r2] * outer(vr[[r]], vr[[r2]])
      }
    }
    list(cross = cross, M = M)
  }
  list(cb = cb, calls = calls,
       data = list(H = H, D3 = D3, Vs = Vs, ap = 1L))
}

test_that("the compiled curvature is the R route's, and never calls back", {
  X <- cbind(1, as.numeric(scale(seq_len(n_c))))
  mb <- ncol(X)
  for (cfg in list(quote(gas(p = 1, q = 1, alpha1 ~ z, time = t)),
                   quote(gas(p = 2, q = 2, omega ~ z, alpha1 ~ g, time = t)),
                   quote(gas(p = 1, q = 1, omega ~ g, by = g, time = t)),
                   quote(gas(p = 1, q = 1, omega ~ z, by = g8, time = t)))) {
    term <- term_build(eval(cfg), dd_c)
    nm <- term_params(term)
    np <- length(nm)
    lk <- term_links(term)
    m <- mb + np
    seed <- cbind(X, matrix(0, n_c, np))
    set.seed(11)
    u0 <- c(0.3, -0.2, stats::runif(np, -0.5, 0.2))
    gw <- stats::rnorm(n_c)
    psi <- as.list(stats::setNames(vapply(seq_len(np), function(j)
      linkfunctions7::linkinv(lk[[nm[j]]], u0[mb + j]), numeric(1)), nm))
    eta0 <- as.numeric(X %*% u0[seq_len(mb)])

    # the lookups are the derivatives read at the predictor the recursion
    # reproduces, which is the layer's own construction
    sc <- gauss_score_c(dd_c$y)
    cu <- gauss_curv_c(dd_c$y)
    e_full <- term_filter(term, eta0, dd_c$y, sc, cu, psi)$eta
    s_at <- dd_c$y - e_full
    c_at <- rep(-1, n_c)
    sc2 <- function(e, i) s_at[i]
    cu2 <- function(e, i) c_at[i]

    lay <- .mk_layer(n_c, m)
    ref <- term_curvature(term, eta0, dd_c$y, sc2, cu2, psi, gw, seed,
                          lay$cb)
    calls_ref <- lay$calls$k
    lay$calls$k <- 0L
    got <- term_curvature(term, eta0, dd_c$y, sc2, cu2, psi, gw, seed,
                          lay$cb, score_values = s_at,
                          curvature_values = c_at,
                          blocks_data = lay$data, threads = 1L)
    lab <- paste(deparse(cfg), collapse = "")
    expect_gt(calls_ref, 0L)
    # the compiled route reads the data and never enters the callback
    expect_identical(lay$calls$k, 0L)
    expect_equal(got$jacobian, ref$jacobian, tolerance = 1e-12, info = lab)
    expect_equal(got$curvature, ref$curvature, tolerance = 1e-12, info = lab)
    expect_identical(got$curvature, t(got$curvature))

    # the same kernel across thread counts, to the bit: the last config has
    # eight groups, which is where the pool opens at all
    got2 <- term_curvature(term, eta0, dd_c$y, sc2, cu2, psi, gw, seed,
                           lay$cb, score_values = s_at,
                           curvature_values = c_at,
                           blocks_data = lay$data, threads = 2L)
    expect_identical(got$jacobian, got2$jacobian)
    expect_identical(got$curvature, got2$curvature)
  }
})

test_that("a developed autoregressive chart keeps the R route", {
  X <- cbind(1, as.numeric(scale(seq_len(n_c))))
  mb <- ncol(X)
  term <- term_build(gas(p = 1, q = 2, pacf1 ~ z, time = t), dd_c)
  nm <- term_params(term)
  np <- length(nm)
  lk <- term_links(term)
  m <- mb + np
  seed <- cbind(X, matrix(0, n_c, np))
  set.seed(12)
  u0 <- c(0.3, -0.2, stats::runif(np, -0.4, 0.1))
  gw <- stats::rnorm(n_c)
  psi <- as.list(stats::setNames(vapply(seq_len(np), function(j)
    linkfunctions7::linkinv(lk[[nm[j]]], u0[mb + j]), numeric(1)), nm))
  eta0 <- as.numeric(X %*% u0[seq_len(mb)])
  sc <- gauss_score_c(dd_c$y)
  cu <- gauss_curv_c(dd_c$y)
  e_full <- term_filter(term, eta0, dd_c$y, sc, cu, psi)$eta
  s_at <- dd_c$y - e_full
  sc2 <- function(e, i) s_at[i]
  cu2 <- function(e, i) -1

  lay <- .mk_layer(n_c, m)
  got <- term_curvature(term, eta0, dd_c$y, sc2, cu2, psi, gw, seed,
                        lay$cb, score_values = s_at,
                        curvature_values = rep(-1, n_c),
                        blocks_data = lay$data, threads = 2L)
  # varying_b: the kernel declines and the callback runs, so the result is
  # the R route's by construction
  expect_gt(lay$calls$k, 0L)
  ref <- term_curvature(term, eta0, dd_c$y, sc2, cu2, psi, gw, seed, lay$cb)
  expect_identical(got$jacobian, ref$jacobian)
  expect_identical(got$curvature, ref$curvature)
})
