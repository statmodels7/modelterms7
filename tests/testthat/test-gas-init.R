# gas(init = "static"): the first predictor of every group is the constant
# init_value, so the first level is c minus the static predictor of the
# group's first row. These tests check that start in every route of the
# scalar term: the filter and its R twin, the Jacobian, the adjoint,
# term_static_deriv(), and term_curvature(), term_third() and term_fourth()
# against numerical differentiation and against the compiled route.

.init_data <- function() {
  set.seed(5)
  n <- 90
  dd <- data.frame(t = seq_len(n), g = factor(rep(c("a", "b", "c"), each = 30)),
                   x = stats::rnorm(n))
  dd$y <- stats::rpois(n, exp(0.8 + 0.3 * sin(seq_len(n) / 6)))
  dd
}

test_that("a static start puts the first predictor of every group at the constant", {
  dd <- .init_data()
  n <- nrow(dd)
  X <- cbind(1, dd$x)
  mb <- ncol(X)
  sc <- function(e, i) dd$y[i] - exp(e)
  cu <- function(e, i) -exp(e)
  cval <- log(mean(dd$y))
  for (cfg in list(quote(gas(p = 1, q = 1, time = t, init = "static")),
                   quote(gas(p = 2, q = 2, time = t, init = "static")),
                   quote(gas(p = 1, q = 1, by = g, time = t, init = "static")),
                   quote(gas(p = 1, q = 3, by = g, time = t, init = "static")))) {
    lab <- paste(deparse(cfg), collapse = "")
    term <- S7::set_props(term_build(eval(cfg), dd), init_value = cval)
    nm <- term_params(term)
    np <- length(nm)
    lk <- term_links(term)
    m <- mb + np
    u0 <- c(0.2, -0.1, 0.05, rep(log(0.05), term@p), rep(0.3, term@q))
    psi_of <- function(u) as.list(stats::setNames(vapply(seq_len(np), function(j)
      linkfunctions7::linkinv(lk[[nm[j]]], u[mb + j]), numeric(1)), nm))
    eta_of <- function(u) {
      term_filter(term, as.numeric(X %*% u[seq_len(mb)]), dd$y, sc, cu,
                  psi_of(u))$eta
    }
    e0 <- as.numeric(X %*% u0[seq_len(mb)])
    out <- term_filter(term, e0, dd$y, sc, cu, psi_of(u0))
    first <- vapply(term@blueprint$order, `[`, 1L, FUN.VALUE = integer(1))
    expect_equal(out$eta[first], rep(cval, length(first)), tolerance = 1e-15,
                 info = lab)
    expect_true(all(out$jacobian[first, ] == 0), info = lab)

    # the R twin
    cf <- .gas_coefs(unlist(psi_of(u0))[.gas_base_params(term@p, term@q)],
                     term@p, term@q)
    tw <- gas_filter_r(e0, term@blueprint$order, term@p, term@q, cf$omega,
                       cf$a, cf$b, cf$db, cf$f0, cf$df0, cf$i_a, cf$np, sc, cu,
                       d0 = .gas_static_start(term, e0))
    expect_equal(tw$eta, out$eta, tolerance = 1e-14, info = lab)
    expect_equal(tw$jacobian, out$jacobian, tolerance = 1e-14, info = lab,
                 ignore_attr = TRUE)

    # the Jacobian in the unknowns and the Hessian of a Poisson
    # log-likelihood, from the second-order recursion
    Jn <- numDeriv::jacobian(eta_of, u0)
    seed <- cbind(X, matrix(0, n, np))
    e_at <- out$eta
    s_at <- dd$y - exp(e_at)
    c_at <- -exp(e_at)
    sc2 <- function(e, i) s_at[i]
    cu2 <- function(e, i) c_at[i]
    cb <- function(e, i, D, act) list(cross = numeric(length(act)),
                                      M = -exp(e_at[i]) * outer(D, D))
    cr <- term_curvature(term, e0, dd$y, sc2, cu2, psi_of(u0), s_at, seed, cb)
    expect_equal(cr$jacobian, Jn, tolerance = 1e-8, info = lab)
    ll <- function(u) { e <- eta_of(u); sum(dd$y * e - exp(e)) }
    H <- crossprod(cr$jacobian, c_at * cr$jacobian) + cr$curvature
    expect_equal(H, numDeriv::hessian(ll, u0), tolerance = 1e-8, info = lab)

    # the adjoint, the scores included
    gfix <- stats::rnorm(n)
    obj_e <- function(e) sum(gfix * term_filter(term, e, dd$y, sc, cu,
                                                psi_of(u0))$eta)
    adj <- term_adjoint(term, e0, dd$y, sc, cu, psi_of(u0), gfix)
    expect_equal(adj$deta, numDeriv::grad(obj_e, e0), tolerance = 1e-7,
                 info = lab)

    # the derivative of the predictor in the equation's coefficients
    sd <- term_static_deriv(term, out$curv, X, psi_of(u0))
    expect_equal(sd, Jn[, seq_len(mb)], tolerance = 1e-8, info = lab,
                 ignore_attr = TRUE)

    # the compiled second order against the R route
    lay <- list(H = matrix(0, n, 1), D3 = matrix(-exp(e_at), n, 1),
                Vs = list(matrix(0, n, m)), ap = 1L)
    gc <- term_curvature(term, e0, dd$y, sc2, cu2, psi_of(u0), s_at, seed, cb,
                         score_values = s_at, curvature_values = c_at,
                         blocks_data = lay)
    expect_equal(gc$curvature, cr$curvature, tolerance = 1e-12, info = lab)
    expect_equal(gc$jacobian, cr$jacobian, tolerance = 1e-12, info = lab)
  }
})

test_that("the third and fourth orders of a static start are the exact directional derivatives", {
  # the stub of test-gas.R: l = -log cosh(e - y), five derivatives non-zero
  # and bounded
  set.seed(8)
  nn <- 50L
  d2 <- data.frame(t = seq_len(nn), y = stats::rnorm(nn),
                   g = factor(rep(c("a", "b"), each = 25)))
  X <- cbind(1, as.numeric(scale(seq_len(nn))) * 0.3)
  mb <- ncol(X)
  tt <- function(e, i) tanh(e - d2$y[i])
  sc <- function(e, i) -tt(e, i)
  cu <- function(e, i) -(1 - tt(e, i)^2)
  l3 <- function(e, i) { t <- tt(e, i); 2 * t * (1 - t^2) }
  l4 <- function(e, i) { t <- tt(e, i); (1 - t^2) * (2 - 6 * t^2) }
  l5 <- function(e, i) { t <- tt(e, i); (1 - t^2) * (24 * t^3 - 16 * t) }
  d2$x <- sin(seq_len(nn) / 4)
  for (cfg in list(quote(gas(p = 1, q = 1, time = t, init = "static")),
                   quote(gas(p = 2, q = 2, time = t, init = "static")),
                   quote(gas(p = 1, q = 3, by = g, time = t, init = "static")),
                   quote(gas(p = 1, q = 1, kappa1 ~ x, time = t, init = "static")),
                   quote(gas(p = 1, q = 2, omega ~ g, by = g, time = t,
                             init = "static")))) {
    lab <- paste(deparse(cfg), collapse = "")
    term <- S7::set_props(term_build(eval(cfg), d2), init_value = 0.2)
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
    eta0 <- function(u) as.numeric(X %*% u[seq_len(mb)])
    blk2 <- function(e, i, D, act) list(cross = numeric(length(D)),
                                        M = l3(e, i) * outer(D, D))
    blk3 <- function(dir) function(e, i, D, act) list(
      cross = numeric(length(D)), M = l3(e, i) * outer(D, D),
      dcurv = l3(e, i) * D,
      N = l4(e, i) * sum(D * dir[act]) * outer(D, D))
    blk4 <- function(d1, d2v) function(e, i, D, act) {
      DD <- outer(D, D)
      dv <- sum(D * d1[act])
      dw <- sum(D * d2v[act])
      list(cross = numeric(length(D)), M = l3(e, i) * DD,
           dcurv = l3(e, i) * D, cppp = l3(e, i), Q = l4(e, i) * DD,
           N = list(l4(e, i) * dv * DD, l4(e, i) * dw * DD),
           P = l5(e, i) * dv * dw * DD)
    }
    psi0 <- psi_of(u0[mb + seq_len(np)])
    rich <- function(at, what, dir) {
      d <- function(h) (at(u0 + h * dir)[[what]] - at(u0 - h * dir)[[what]]) / (2 * h)
      (4 * d(1e-3) - d(2e-3)) / 3
    }
    # third order: the derivative of the second along v
    got3 <- term_third(term, eta0(u0), d2$y, sc, cu, psi0, gw, seed, blk3(v), v)
    at2 <- function(u) term_curvature(term, eta0(u), d2$y, sc, cu,
      psi_of(u[mb + seq_len(np)]), gw, seed, blk2)
    expect_equal(got3$curvature, rich(at2, "curvature", v), tolerance = 1e-6,
                 info = lab)
    expect_equal(got3$dphi, rich(at2, "jacobian", v), tolerance = 1e-6,
                 info = lab)
    # fourth order: the derivative of the third along w
    got4 <- term_fourth(term, eta0(u0), d2$y, sc, cu, psi0, gw, seed,
                        blk4(v, w), list(v, w))
    at3 <- function(u) term_third(term, eta0(u), d2$y, sc, cu,
      psi_of(u[mb + seq_len(np)]), gw, seed, blk3(v), v)
    expect_equal(got4$curvature, rich(at3, "curvature", w), tolerance = 1e-6,
                 info = lab)
    expect_equal(got4$dpsi, rich(at3, "dphi", w), tolerance = 1e-6, info = lab)

    # the compiled route at orders 3 and 4 against the R route, the model's
    # pieces read at the filtered predictor
    e_at <- term_filter(term, eta0(u0), d2$y, sc, cu, psi0)$eta
    s_at <- vapply(seq_len(nn), function(i) sc(e_at[i], i), numeric(1))
    c_at <- vapply(seq_len(nn), function(i) cu(e_at[i], i), numeric(1))
    f3 <- vapply(seq_len(nn), function(i) l3(e_at[i], i), numeric(1))
    f4 <- vapply(seq_len(nn), function(i) l4(e_at[i], i), numeric(1))
    f5 <- vapply(seq_len(nn), function(i) l5(e_at[i], i), numeric(1))
    sc2 <- function(e, i) s_at[i]
    cu2 <- function(e, i) c_at[i]
    lay <- list(H = matrix(0, nn, 1), D3 = matrix(f3), D4 = matrix(f4),
                D5 = matrix(f5), Vs = list(matrix(0, nn, m)), ap = 1L)
    b3 <- function(e, i, D, act) list(
      cross = numeric(length(D)), M = f3[i] * outer(D, D), dcurv = f3[i] * D,
      N = f4[i] * sum(D * v[act]) * outer(D, D))
    b4 <- function(e, i, D, act) {
      DD <- outer(D, D)
      dv <- sum(D * v[act])
      dw <- sum(D * w[act])
      list(cross = numeric(length(D)), M = f3[i] * DD, dcurv = f3[i] * D,
           cppp = f3[i], Q = f4[i] * DD,
           N = list(f4[i] * dv * DD, f4[i] * dw * DD), P = f5[i] * dv * dw * DD)
    }
    r3 <- term_third(term, eta0(u0), d2$y, sc2, cu2, psi0, gw, seed, b3, v)
    c3 <- term_third(term, eta0(u0), d2$y, sc2, cu2, psi0, gw, seed, b3, v,
                     score_values = s_at, curvature_values = c_at,
                     blocks_data = lay)
    expect_equal(c3$curvature, r3$curvature, tolerance = 1e-12, info = lab)
    expect_equal(c3$dphi, r3$dphi, tolerance = 1e-12, info = lab)
    r4 <- term_fourth(term, eta0(u0), d2$y, sc2, cu2, psi0, gw, seed, b4,
                      list(v, w))
    c4 <- term_fourth(term, eta0(u0), d2$y, sc2, cu2, psi0, gw, seed, b4,
                      list(v, w), score_values = s_at, curvature_values = c_at,
                      blocks_data = lay)
    expect_equal(c4$curvature, r4$curvature, tolerance = 1e-12, info = lab)
    expect_equal(c4$dpsi, r4$dpsi, tolerance = 1e-12, info = lab)
  }
})

test_that("a static start with subformulas has the derivatives of its predictor", {
  dd <- .init_data()
  n <- nrow(dd)
  X <- cbind(1, dd$x)
  mb <- ncol(X)
  sc <- function(e, i) dd$y[i] - exp(e)
  cu <- function(e, i) -exp(e)
  for (cfg in list(quote(gas(p = 1, q = 1, kappa1 ~ x, time = t, init = "static")),
                   quote(gas(p = 1, q = 2, omega ~ g, pacf1 ~ x, by = g, time = t,
                             init = "static")))) {
    lab <- paste(deparse(cfg), collapse = "")
    term <- term_build(eval(cfg), dd)
    ng <- length(term@blueprint$order)
    cval <- log(mean(dd$y)) + 0.05 * seq_len(ng)
    term <- S7::set_props(term, init_value = cval)
    nm <- term_params(term)
    np <- length(nm)
    lk <- term_links(term)
    m <- mb + np
    set.seed(3)
    u0 <- c(0.2, -0.1, stats::runif(np, -0.2, 0.2))
    u0[mb + grep("^kappa1", nm)[1L]] <- log(0.05)
    psi_of <- function(u) as.list(stats::setNames(vapply(seq_len(np), function(j)
      linkfunctions7::linkinv(lk[[nm[j]]], u[mb + j]), numeric(1)), nm))
    eta_of <- function(u) {
      term_filter(term, as.numeric(X %*% u[seq_len(mb)]), dd$y, sc, cu,
                  psi_of(u))$eta
    }
    e0 <- as.numeric(X %*% u0[seq_len(mb)])
    out <- term_filter(term, e0, dd$y, sc, cu, psi_of(u0))
    first <- vapply(term@blueprint$order, `[`, 1L, FUN.VALUE = integer(1))
    expect_equal(out$eta[first], cval, tolerance = 1e-15, info = lab)
    tw <- .gas_filter_sub_r(term, e0, dd$y, sc, cu, unlist(psi_of(u0)))
    expect_equal(tw$eta, out$eta, tolerance = 1e-13, info = lab)
    expect_equal(tw$jacobian, out$jacobian, tolerance = 1e-12, info = lab,
                 ignore_attr = TRUE)
    Jn <- numDeriv::jacobian(eta_of, u0)
    seed <- cbind(X, matrix(0, n, np))
    e_at <- out$eta
    s_at <- dd$y - exp(e_at)
    c_at <- -exp(e_at)
    sc2 <- function(e, i) s_at[i]
    cu2 <- function(e, i) c_at[i]
    cb <- function(e, i, D, act) list(cross = numeric(length(act)),
                                      M = -exp(e_at[i]) * outer(D, D))
    cr <- term_curvature(term, e0, dd$y, sc2, cu2, psi_of(u0), s_at, seed, cb)
    expect_equal(cr$jacobian, Jn, tolerance = 1e-7, info = lab,
                 ignore_attr = TRUE)
    ll <- function(u) { e <- eta_of(u); sum(dd$y * e - exp(e)) }
    H <- crossprod(cr$jacobian, c_at * cr$jacobian) + cr$curvature
    expect_equal(H, numDeriv::hessian(ll, u0), tolerance = 1e-7, info = lab,
                 ignore_attr = TRUE)
    gfix <- stats::rnorm(n)
    obj_e <- function(e) sum(gfix * term_filter(term, e, dd$y, sc, cu,
                                                psi_of(u0))$eta)
    adj <- term_adjoint(term, e0, dd$y, sc, cu, psi_of(u0), gfix)
    expect_equal(adj$deta, numDeriv::grad(obj_e, e0), tolerance = 1e-7,
                 info = lab)
    lay <- list(H = matrix(0, n, 1), D3 = matrix(-exp(e_at), n, 1),
                Vs = list(matrix(0, n, m)), ap = 1L)
    gc <- term_curvature(term, e0, dd$y, sc2, cu2, psi_of(u0), s_at, seed, cb,
                         score_values = s_at, curvature_values = c_at,
                         blocks_data = lay)
    expect_equal(gc$curvature, cr$curvature, tolerance = 1e-12, info = lab)
    expect_equal(gc$jacobian, cr$jacobian, tolerance = 1e-12, info = lab)
  }
})

test_that("a static start needs its constant", {
  dd <- .init_data()
  term <- term_build(gas(p = 1, q = 1, time = t, init = "static"), dd)
  expect_error(term_filter(term, rep(0, nrow(dd)), dd$y,
                           function(e, i) dd$y[i] - exp(e),
                           function(e, i) -exp(e),
                           list(omega = 0, kappa1 = 0.05, pacf1 = 0.3)),
               "init_value")
  expect_identical(gas()@init, "level")
})
