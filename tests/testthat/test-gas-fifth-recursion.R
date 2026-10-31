# The fifth-order recursion against a Richardson difference of term_fourth()
# along the third direction. The layer is stubbed by l = -log cosh(e - y),
# whose six derivatives are all non-zero and bounded, as in the fourth
# order's own test one order up.

tanh_stub <- function(y) {
  tt <- function(e, i) tanh(e - y[i])
  list(
    sc = function(e, i) -tt(e, i),
    cu = function(e, i) -(1 - tt(e, i)^2),
    l3 = function(e, i) { t <- tt(e, i); 2 * t * (1 - t^2) },
    l4 = function(e, i) { t <- tt(e, i); (1 - t^2) * (2 - 6 * t^2) },
    l5 = function(e, i) { t <- tt(e, i); (1 - t^2) * (24 * t^3 - 16 * t) },
    l6 = function(e, i) {
      t <- tt(e, i)
      (1 - t^2) * (-2 * t * (24 * t^3 - 16 * t) + (1 - t^2) * (72 * t^2 - 16))
    })
}

# the one-equation blocks of term_fourth(), and the arrays of term_fifth()
blk4_of <- function(st, d1, d2v) function(e, i, D, act) {
  DD <- outer(D, D)
  dv <- sum(D * d1[act])
  dw <- sum(D * d2v[act])
  list(cross = numeric(length(D)), M = st$l3(e, i) * DD,
       dcurv = st$l3(e, i) * D, cppp = st$l3(e, i), Q = st$l4(e, i) * DD,
       N = list(st$l4(e, i) * dv * DD, st$l4(e, i) * dw * DD),
       P = st$l5(e, i) * dv * dw * DD)
}
arrays_of <- function(st, e, m) {
  idx <- seq_along(e)
  list(H = matrix(st$cu(e, idx)), D3 = matrix(st$l3(e, idx)),
       D4 = matrix(st$l4(e, idx)), D5 = matrix(st$l5(e, idx)),
       D6 = matrix(st$l6(e, idx)), Vs = list(matrix(0, length(e), m)),
       ap = 1L)
}

test_that("the fifth-order recursion is the derivative of the fourth", {
  set.seed(8)
  nn <- 40L
  d2 <- data.frame(t = seq_len(nn), y = stats::rnorm(nn))
  X <- cbind(1, as.numeric(scale(seq_len(nn))) * 0.3)
  mb <- ncol(X)
  st <- tanh_stub(d2$y)
  for (cfg in list(c(1, 1), c(2, 2), c(1, 5))) {
    term <- term_build(gas(p = cfg[1], q = cfg[2], time = t), d2)
    nm <- term_params(term)
    np <- length(nm)
    lk <- term_links(term)
    m <- mb + np
    seed <- cbind(X, matrix(0, nn, np))
    psi_of <- function(z) as.list(stats::setNames(vapply(seq_along(nm),
      function(j) linkfunctions7::linkinv(lk[[nm[j]]], z[j]), numeric(1)), nm))
    set.seed(9)
    u0 <- c(0.3, -0.2, stats::runif(np, -0.4, 0.2))
    gw <- stats::rnorm(nn)
    v <- stats::rnorm(m)
    w <- stats::rnorm(m)
    x <- stats::rnorm(m)
    eta0 <- function(u) as.numeric(X %*% u[seq_len(mb)])
    filt <- term_filter(term, eta0(u0), d2$y, st$sc, st$cu,
                        psi_of(u0[mb + seq_len(np)]))
    e <- filt$eta
    got <- term_fifth(term, eta0(u0), d2$y, st$sc,
                      psi_of(u0[mb + seq_len(np)]), gw, seed,
                      arrays_of(st, e, m), list(v, w, x))
    at <- function(u) term_fourth(term, eta0(u), d2$y, st$sc, st$cu,
      psi_of(u[mb + seq_len(np)]), gw, seed, blk4_of(st, v, w), list(v, w))
    dir_of <- function(h, f) (f(at(u0 + h * x)) - f(at(u0 - h * x))) / (2 * h)
    rich <- function(f) (4 * dir_of(1e-3, f) - dir_of(2e-3, f)) / 3
    lab <- sprintf("p = %d, q = %d", cfg[1], cfg[2])
    ref <- rich(function(r) r$curvature)
    expect_equal(got$curvature, ref, tolerance = 1e-6, info = lab)
    expect_equal(got$domega, rich(function(r) r$dpsi), tolerance = 1e-6,
                 info = lab)
    # the lower orders agree with term_fourth() at the point
    f4 <- at(u0)
    expect_equal(got$jacobian, f4$jacobian, tolerance = 1e-12, info = lab)
    expect_equal(got$dphi[[1L]], f4$dphi[[1L]], tolerance = 1e-12, info = lab)
    expect_equal(got$dphi[[2L]], f4$dphi[[2L]], tolerance = 1e-12, info = lab)
    expect_equal(got$dpsi[[1L]], f4$dpsi, tolerance = 1e-12, info = lab)
    # symmetric in the three directions, and as a matrix
    sw <- term_fifth(term, eta0(u0), d2$y, st$sc,
                     psi_of(u0[mb + seq_len(np)]), gw, seed,
                     arrays_of(st, e, m), list(x, v, w))
    expect_equal(sw$curvature, got$curvature, tolerance = 1e-12, info = lab)
    expect_equal(got$curvature, t(got$curvature), tolerance = 1e-14)
    # the negative control: a wrong sixth derivative is detected
    bad <- arrays_of(st, e, m)
    bad$D6 <- bad$D6 * 1.05
    wrong <- term_fifth(term, eta0(u0), d2$y, st$sc,
                        psi_of(u0[mb + seq_len(np)]), gw, seed, bad,
                        list(v, w, x))
    expect_false(isTRUE(all.equal(wrong$curvature, ref, tolerance = 1e-6)),
                 info = lab)
  }
})

test_that("term_fifth is zero for an additive term and refused where absent", {
  td <- data.frame(x = stats::rnorm(20), y = stats::rnorm(20))
  lp <- term_build(linpar(~x), td)
  out <- term_fifth(lp, rep(0, 20), td$y, function(e, i) 0, list(),
                    rep(1, 20), matrix(0, 20, 2), NULL,
                    list(c(1, 0), c(0, 1), c(1, 1)))
  expect_true(all(out$curvature == 0))
  expect_length(out$dpsi, 3L)
  rg <- term_build(regime(k = 2), td)
  expect_error(term_fifth(rg, rep(0, 20), td$y, function(e, i) 0, list(),
                          rep(1, 20), matrix(0, 20, 2), NULL,
                          list(c(1, 0), c(0, 1), c(1, 1))),
               "does not implement term_fifth")
  gt <- term_build(gas(p = 1, q = 1), td)
  expect_error(term_fifth(gt, rep(0, 20), td$y, function(e, i) 0,
                          list(omega = 0.1, kappa1 = 0.2, pacf1 = 0.3),
                          rep(1, 20), matrix(0, 20, 3), list(),
                          list(c(1, 0, 0), c(0, 1, 0))),
               "list of three")
})

test_that("the fifth order reads the other equation's predictor", {
  # Two predictors: l(e1, e2) = -exp(e2) log cosh(e1 - y) - e2^2 / 2, the
  # filter in e1 and e2 = Z gamma linear. Every derivative with e1 among its
  # indices is -exp(e2) L^(k)(e1 - y), with k the number of e1's, so the
  # cross derivatives of every order are non-zero.
  set.seed(18)
  nn <- 40L
  d2 <- data.frame(t = seq_len(nn), y = stats::rnorm(nn))
  X <- cbind(1, as.numeric(scale(seq_len(nn))) * 0.3)
  Z <- cbind(1, stats::runif(nn, -1, 1))
  st <- tanh_stub(d2$y)
  # L^(k) of log cosh at z = e1 - y, minus the stub's l_k
  Lk <- function(k, z) {
    t <- tanh(z)
    -switch(k, -t, -(1 - t^2), 2 * t * (1 - t^2), (1 - t^2) * (2 - 6 * t^2),
            (1 - t^2) * (24 * t^3 - 16 * t),
            (1 - t^2) * (-2 * t * (24 * t^3 - 16 * t) +
                           (1 - t^2) * (72 * t^2 - 16)))
  }
  term <- term_build(gas(p = 1, q = 2, time = t), d2)
  nm <- term_params(term)
  np <- length(nm)
  lk <- term_links(term)
  m <- 2L + 2L + np
  ib <- 1:2
  ig <- 3:4
  seed <- cbind(X, matrix(0, nn, 2L + np))
  Vz <- cbind(matrix(0, nn, 2L), Z, matrix(0, nn, np))
  psi_of <- function(z) as.list(stats::setNames(vapply(seq_along(nm),
    function(j) linkfunctions7::linkinv(lk[[nm[j]]], z[j]), numeric(1)), nm))
  set.seed(19)
  u0 <- c(0.3, -0.2, -0.3, 0.2, stats::runif(np, -0.4, 0.2))
  gw <- stats::rnorm(nn)
  v <- stats::rnorm(m)
  w <- stats::rnorm(m)
  x <- stats::rnorm(m)
  # ell_{1, r_1..r_k} at the predictors: -exp(e2) L^(1 + #1)(e1 - y)
  ell_arr <- function(k, e1, e2, yy = d2$y) {
    idx <- as.matrix(do.call(expand.grid, rep(list(1:2), k)))
    out <- matrix(0, length(e1), nrow(idx))
    for (c in seq_len(nrow(idx))) {
      out[, c] <- -exp(e2) * Lk(1 + sum(idx[c, ] == 1L), e1 - yy)
    }
    out
  }
  cb <- function(u) {
    e2 <- as.numeric(Z %*% u[ig])
    list(e2 = e2,
         sc = function(e, i) -exp(e2[i]) * Lk(1, e - d2$y[i]),
         cu = function(e, i) -exp(e2[i]) * Lk(2, e - d2$y[i]),
         ell = function(k, e, i) {
           a <- ell_arr(k, e, e2[i], d2$y[i])
           array(a, rep(2L, k))
         })
  }
  # the callback of term_fourth() for two equations, from the same arrays
  blk4 <- function(cbu, d1, d2v) function(e, i, D, act) {
    V <- rbind(D, Vz[i, act])
    l2 <- cbu$ell(1, e, i); l3 <- cbu$ell(2, e, i)
    l4 <- cbu$ell(3, e, i); l5 <- cbu$ell(4, e, i)
    dv <- as.numeric(V %*% d1[act])
    dw <- as.numeric(V %*% d2v[act])
    quad <- function(C) crossprod(V, C %*% V)
    list(cross = l2[2] * Vz[i, act],
         M = quad(l3),
         dcurv = as.numeric(crossprod(V, l3[1, ])),
         cppp = l3[1, 1],
         Q = quad(l4[1, , ]),
         N = list(quad(apply(l4, c(1, 2), function(z) sum(z * dv))),
                  quad(apply(l4, c(1, 2), function(z) sum(z * dw)))),
         P = quad(apply(l5, c(1, 2), function(z) sum(z * outer(dv, dw)))))
  }
  eta0 <- function(u) as.numeric(X %*% u[ib])
  at <- function(u) {
    cbu <- cb(u)
    term_fourth(term, eta0(u), d2$y, cbu$sc, cbu$cu,
                psi_of(u[4L + seq_len(np)]), gw, seed, blk4(cbu, v, w),
                list(v, w))
  }
  c0 <- cb(u0)
  e1 <- term_filter(term, eta0(u0), d2$y, c0$sc, c0$cu,
                    psi_of(u0[4L + seq_len(np)]))$eta
  arr <- list(H = ell_arr(1, e1, c0$e2), D3 = ell_arr(2, e1, c0$e2),
              D4 = ell_arr(3, e1, c0$e2), D5 = ell_arr(4, e1, c0$e2),
              D6 = ell_arr(5, e1, c0$e2),
              Vs = list(matrix(0, nn, m), Vz), ap = 1L)
  got <- term_fifth(term, eta0(u0), d2$y, c0$sc,
                    psi_of(u0[4L + seq_len(np)]), gw, seed, arr,
                    list(v, w, x))
  dir_of <- function(h, f) (f(at(u0 + h * x)) - f(at(u0 - h * x))) / (2 * h)
  rich <- function(f) (4 * dir_of(1e-3, f) - dir_of(2e-3, f)) / 3
  ref <- rich(function(r) r$curvature)
  expect_equal(got$curvature, ref, tolerance = 1e-6)
  expect_equal(got$domega, rich(function(r) r$dpsi), tolerance = 1e-6)
  f4 <- at(u0)
  expect_equal(got$dpsi[[1L]], f4$dpsi, tolerance = 1e-12)
  # the negative control: the other equation's sixth-order entries dropped
  bad <- arr
  bad$D6[, -1L] <- 0
  wrong <- term_fifth(term, eta0(u0), d2$y, c0$sc,
                      psi_of(u0[4L + seq_len(np)]), gw, seed, bad,
                      list(v, w, x))
  expect_false(isTRUE(all.equal(wrong$curvature, ref, tolerance = 1e-6)))
})
