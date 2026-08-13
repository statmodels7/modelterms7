# The submodel route of gas(): a parameter developed with covariates, the
# recursion read per observation, and every derivative against numDeriv or
# against the special case it generalizes.

set.seed(15)
n <- 60
dd <- data.frame(t = seq_len(n), y = rnorm(n),
                 z = as.numeric(scale(runif(n))),
                 g = factor(rep(c("a", "b"), each = n / 2)))

gauss_score <- function(y) function(e, i) y[i] - e
gauss_curv <- function(y) function(e, i) -1

test_that("the syntax is validated where it is written", {
  expect_error(gas(p = 1, q = 1, gamma ~ z), "the parameters are")
  expect_error(gas(p = 1, q = 1, alpha1 ~ z, alpha1 ~ g),
               "two subformulas")
  expect_error(gas(p = 1, q = 1, ~z), "two-sided formulas")
  expect_error(gas(p = 1, q = 1, by = ~ridge(~g), alpha1 ~ z),
               "mixing it with per-parameter")
  # by as a formula is the shorthand: every parameter gets the subformula
  sh <- gas(p = 1, q = 1, by = ~z)
  expect_identical(names(sh@submodels), c("omega", "alpha1", "pacf1"))
})

test_that("a development of ~1 is the scalar filter in other coordinates", {
  term_sc <- term_build(gas(p = 1, q = 1, time = t), dd)
  psi <- list(omega = 0.2, alpha1 = 0.3, pacf1 = 0.5)
  sc <- term_filter(term_sc, rep(0, n), dd$y, gauss_score(dd$y),
                    gauss_curv(dd$y), psi)

  # omega ~ 1 on the identity chart: the coordinate IS omega
  t_om <- term_build(gas(p = 1, q = 1, omega ~ 1, time = t), dd)
  expect_identical(term_params(t_om),
                   c("omega.(Intercept)", "alpha1", "pacf1"))
  ps <- list("omega.(Intercept)" = 0.2, alpha1 = 0.3, pacf1 = 0.5)
  got <- term_filter(t_om, rep(0, n), dd$y, gauss_score(dd$y),
                     gauss_curv(dd$y), ps)
  expect_equal(got$eta, sc$eta, tolerance = 1e-12)
  expect_equal(unname(got$jacobian), unname(sc$jacobian), tolerance = 1e-10)

  # alpha1 ~ 1 rides the log chart inside, so its coordinate is
  # log(alpha1) and the jacobian column carries the chain factor
  t_al <- term_build(gas(p = 1, q = 1, alpha1 ~ 1, time = t), dd)
  ps2 <- list(omega = 0.2, "alpha1.(Intercept)" = log(0.3), pacf1 = 0.5)
  got2 <- term_filter(t_al, rep(0, n), dd$y, gauss_score(dd$y),
                      gauss_curv(dd$y), ps2)
  expect_equal(got2$eta, sc$eta, tolerance = 1e-12)
  expect_equal(unname(got2$jacobian[, 2L]), unname(sc$jacobian[, 2L]) * 0.3,
               tolerance = 1e-10)
})

test_that("a time-varying loading reproduces a recursion written by hand", {
  term <- term_build(gas(p = 1, q = 1, alpha1 ~ z, time = t), dd)
  g0 <- log(0.2)
  g1 <- 0.4
  om <- 0.1
  rho <- 0.5
  ps <- list(omega = om, "alpha1.(Intercept)" = g0, "alpha1.z" = g1,
             pacf1 = rho)
  got <- term_filter(term, rep(0, n), dd$y, gauss_score(dd$y),
                     gauss_curv(dd$y), ps)

  # the same recursion, written independently: the loading of each
  # observation from its own covariate, positive by the chart
  a_t <- exp(g0 + g1 * dd$z)
  f <- numeric(n)
  s <- numeric(n)
  f0 <- om / (1 - rho)
  for (tt in seq_len(n)) {
    f[tt] <- om + a_t[tt] * (if (tt > 1) s[tt - 1] else 0) +
      rho * (if (tt > 1) f[tt - 1] else f0)
    s[tt] <- dd$y[tt] - f[tt]
  }
  expect_equal(got$eta, f, tolerance = 1e-12)
  expect_true(all(a_t > 0))
})

test_that("the jacobian of a developed filter is exact", {
  for (cfg in list(list(f = quote(gas(p = 1, q = 1, alpha1 ~ z, time = t))),
                   list(f = quote(gas(p = 1, q = 1, omega ~ g,
                                      pacf1 ~ z, time = t))),
                   list(f = quote(gas(p = 2, q = 2, alpha2 ~ z,
                                      time = t))))) {
    term <- term_build(eval(cfg$f), dd)
    nm <- term_params(term)
    np <- length(nm)
    lk <- term_links(term)
    set.seed(7)
    u0 <- stats::runif(np, -0.6, 0.3)
    psi_of <- function(v) as.list(stats::setNames(vapply(seq_len(np),
      function(j) linkfunctions7::linkinv(lk[[nm[j]]], v[j]), numeric(1)),
      nm))
    got <- term_filter(term, rep(0, n), dd$y, gauss_score(dd$y),
                       gauss_curv(dd$y), psi_of(u0))
    # the numerical reference differentiates in the same coordinates the
    # filter reports, the links being identity for every development
    e_of <- function(v) term_filter(term, rep(0, n), dd$y,
                                    gauss_score(dd$y), gauss_curv(dd$y),
                                    psi_of(v))$eta
    num <- numDeriv::jacobian(e_of, u0)
    # chain the psi-scale jacobian onto the coordinates, as the layer does
    chain <- vapply(seq_len(np), function(j)
      linkfunctions7::dlinkinv(lk[[nm[j]]], u0[j]), numeric(1))
    expect_equal(unname(got$jacobian * rep(chain, each = n)), num,
                 tolerance = 1e-6, info = deparse(cfg$f))
  }
})

test_that("a per-group development is the population-and-departures model", {
  # omega ~ random(~1|g) with by = g: each group's level is the intercept
  # plus its own departure, on the parameter's unconstrained scale, and
  # each group is filtered independently. The reference is the per-group
  # recursion written by hand, sharing no code with the route.
  t_sub <- term_build(gas(p = 1, q = 1, omega ~ random(~1 | g), by = g,
                          time = t), dd)
  nm_sub <- term_params(t_sub)
  expect_length(nm_sub, 5L)

  om <- 0.2; al <- 0.3; rho <- 0.5; d_a <- 0.15; d_b <- -0.1
  ps_sub <- stats::setNames(list(om, d_a, d_b, al, rho), nm_sub)
  b <- term_filter(t_sub, rep(0, n), dd$y, gauss_score(dd$y),
                   gauss_curv(dd$y), ps_sub)
  hand <- numeric(n)
  for (l in c("a", "b")) {
    rows <- which(dd$g == l)
    om_l <- om + if (l == "a") d_a else d_b
    f <- numeric(length(rows)); s <- numeric(length(rows))
    f0 <- om_l / (1 - rho)
    for (tt in seq_along(rows)) {
      f[tt] <- om_l + al * (if (tt > 1) s[tt - 1] else 0) +
        rho * (if (tt > 1) f[tt - 1] else f0)
      s[tt] <- dd$y[rows[tt]] - f[tt]
    }
    hand[rows] <- f
  }
  expect_equal(b$eta, hand, tolerance = 1e-12)

  # on the log chart the same model: a loading developed per group is
  # exp(gamma0 + delta) = alpha_pop * exp(delta), positive whatever the
  # departure is -- the property the development exists for
  t_sl <- term_build(gas(p = 1, q = 1, alpha1 ~ random(~1 | g), by = g,
                         time = t), dd)
  ps_sl <- stats::setNames(list(om, log(al), d_a, d_b, rho),
                           term_params(t_sl))
  b2 <- term_filter(t_sl, rep(0, n), dd$y, gauss_score(dd$y),
                    gauss_curv(dd$y), ps_sl)
  hand2 <- numeric(n)
  for (l in c("a", "b")) {
    rows <- which(dd$g == l)
    al_l <- exp(log(al) + if (l == "a") d_a else d_b)
    expect_gt(al_l, 0)
    f <- numeric(length(rows)); s <- numeric(length(rows))
    f0 <- om / (1 - rho)
    for (tt in seq_along(rows)) {
      f[tt] <- om + al_l * (if (tt > 1) s[tt - 1] else 0) +
        rho * (if (tt > 1) f[tt - 1] else f0)
      s[tt] <- dd$y[rows[tt]] - f[tt]
    }
    hand2[rows] <- f
  }
  expect_equal(b2$eta, hand2, tolerance = 1e-12)
})

test_that("the compiled recursion is the R route's twin", {
  # .gas_filter_sub_r is the same loop in R, sharing no code with the
  # kernel; the comparison carries a tolerance rather than asking identity,
  # a compiler being free to contract a multiply-add (the seg_block lesson)
  for (cfg in list(quote(gas(p = 1, q = 1, alpha1 ~ z, time = t)),
                   quote(gas(p = 1, q = 1, omega ~ random(~1 | g),
                             pacf1 ~ z, by = g, time = t)),
                   quote(gas(p = 2, q = 2, alpha2 ~ z, omega ~ g,
                             time = t)))) {
    term <- term_build(eval(cfg), dd)
    nm <- term_params(term)
    np <- length(nm)
    lk <- term_links(term)
    set.seed(21)
    u0 <- stats::runif(np, -0.5, 0.2)
    psi <- stats::setNames(vapply(seq_len(np), function(j)
      linkfunctions7::linkinv(lk[[nm[j]]], u0[j]), numeric(1)), nm)
    got <- .gas_filter_sub(term, rep(0.05, n), dd$y, gauss_score(dd$y),
                           gauss_curv(dd$y), psi)
    ref <- .gas_filter_sub_r(term, rep(0.05, n), dd$y, gauss_score(dd$y),
                             gauss_curv(dd$y), psi)
    lab <- deparse(cfg)
    expect_equal(got$eta, ref$eta, tolerance = 1e-13, info = lab)
    expect_equal(got$jacobian, ref$jacobian, tolerance = 1e-13, info = lab)
  }
})

test_that("the adjoint of a developed filter matches the definition", {
  term <- term_build(gas(p = 1, q = 1, alpha1 ~ z, time = t), dd)
  ps <- list(omega = 0.1, "alpha1.(Intercept)" = log(0.25),
             "alpha1.z" = 0.3, pacf1 = 0.4)
  set.seed(9)
  gw <- rnorm(n)
  out <- term_adjoint(term, rep(0, n), dd$y, gauss_score(dd$y),
                      gauss_curv(dd$y), ps, g = gw)
  # dL/d eta0 with L = sum(gw * eta): the static predictor enters every
  # later score, which is what the reverse pass accounts for
  L_of <- function(e0) sum(gw * term_filter(term, e0, dd$y,
                                            gauss_score(dd$y),
                                            gauss_curv(dd$y), ps)$eta)
  expect_equal(out$deta, numDeriv::grad(L_of, rep(0, n)), tolerance = 1e-6)
})

test_that("the curvature of a developed filter is exact", {
  X <- cbind(1, as.numeric(scale(seq_len(n))))
  mb <- ncol(X)
  sc <- gauss_score(dd$y)
  cu <- gauss_curv(dd$y)
  for (cfg in list(quote(gas(p = 1, q = 1, alpha1 ~ z, time = t)),
                   quote(gas(p = 1, q = 1, omega ~ g, time = t)),
                   quote(gas(p = 1, q = 2, pacf1 ~ z, time = t)))) {
    term <- term_build(eval(cfg), dd)
    nm <- term_params(term)
    np <- length(nm)
    lk <- term_links(term)
    m <- mb + np
    seed <- cbind(X, matrix(0, n, np))
    blocks <- function(e, i, D, act) list(cross = numeric(m),
                                          M = matrix(0, m, m))
    psi_of <- function(z) as.list(stats::setNames(vapply(seq_len(np),
      function(j) linkfunctions7::linkinv(lk[[nm[j]]], z[j]), numeric(1)),
      nm))
    set.seed(11)
    u0 <- c(0.3, -0.2, stats::runif(np, -0.5, 0.2))
    gw <- stats::rnorm(n)
    got <- term_curvature(term, as.numeric(X %*% u0[seq_len(mb)]), dd$y,
                          sc, cu, psi_of(u0[mb + seq_len(np)]), gw, seed,
                          blocks)
    e_of <- function(u) {
      term_filter(term, as.numeric(X %*% u[seq_len(mb)]), dd$y, sc, cu,
                  psi_of(u[mb + seq_len(np)]))$eta
    }
    lab <- deparse(cfg)
    expect_equal(got$jacobian, numDeriv::jacobian(e_of, u0),
                 tolerance = 1e-6, info = lab)
    gsum <- function(u) sum(gw * e_of(u))
    expect_equal(got$curvature, numDeriv::hessian(gsum, u0),
                 tolerance = 1e-6, info = lab)
    expect_identical(got$curvature, t(got$curvature))
  }
})

test_that("the contract views answer for a developed term", {
  term <- term_build(gas(p = 1, q = 1, omega ~ random(~1 | g),
                         alpha1 ~ z, by = g, time = t), dd)
  # penalties: the random intercept's, keyed parameter::subterm, over the
  # indicator coordinates in the term's own numbering
  ent <- term_penalties(term)
  expect_length(ent, 1L)
  expect_true(startsWith(ent[[1L]]$name, "omega::"))
  expect_true(all(startsWith(term_params(term)[ent[[1L]]$index],
                             "omega.random")))
  # the start projects each parameter's own target onto its design: omega
  # at zero, the loading at 0.1 through its log chart
  z0 <- term_start(term)
  expect_equal(unname(z0[["alpha1.(Intercept)"]]), log(0.1),
               tolerance = 1e-8)
  expect_equal(unname(z0[["omega.(Intercept)"]]), 0, tolerance = 1e-8)
  # the level's constant coordinate is the intercept of its development
  expect_identical(term_level_param(term), "omega.(Intercept)")
  # every coordinate of a development is unconstrained
  lk <- term_links(term)
  expect_true(all(vapply(lk[startsWith(names(lk), "omega.")],
                         function(l) l@link_name == "identity", logical(1))))
  expect_output(print(term), "developed")
})
