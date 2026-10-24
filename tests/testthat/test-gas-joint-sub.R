# A joint filter whose terms develop their parameters with subformulas: the
# coefficients of the recursion vary by observation. The references are the
# single term's submodel route, the joint filter without the subformula, and
# differences of the filter.

set.seed(23)
ns <- 80
ds <- data.frame(t = rep(seq_len(ns / 2), 2), y = rnorm(ns, 1, 1.5),
                 g = rep(c("a", "b"), each = ns / 2),
                 x = rnorm(ns))

gs_score <- function(y) function(e, i) {
  s2 <- exp(2 * e[2])
  r <- y[i] - e[1]
  c(r / s2, r^2 / s2 - 1)
}
gs_curv <- function(y) function(e, i) {
  s2 <- exp(2 * e[2])
  r <- y[i] - e[1]
  matrix(c(-1 / s2, -2 * r / s2, -2 * r / s2, -2 * r^2 / s2), 2, 2)
}

# values for every entry of term_params(): a developed coordinate is
# unconstrained, a scalar parameter is on its own scale
gs_psi <- function(nm) {
  v <- stats::setNames(numeric(length(nm)), nm)
  for (k in nm) {
    v[[k]] <- if (grepl("omega", k)) 0.3
      else if (grepl("kappa1\\.\\(Intercept\\)", k)) log(0.2)
      else if (grepl("kappa1\\.", k)) 0.15
      else if (grepl("_", k)) 0.05
      else if (grepl("kappa", k)) 0.2
      else if (grepl("pacf", k)) 0.5
      else 0
  }
  as.list(v)
}

test_that("one developed term filtered jointly is the term's own submodel route", {
  tm <- term_build(gas(p = 1, q = 1, kappa1 ~ g, by = g, time = t), ds)
  jt <- gas_joint(list(mu = tm))
  psi <- gs_psi(term_params(tm))
  sc1 <- function(e, i) (ds$y[i] - e) / 2
  cv1 <- function(e, i) -1 / 2
  eta0 <- rep(0.1, ns)
  a <- term_filter(tm, eta0, ds$y, sc1, cv1, psi)
  b <- term_filter(jt, matrix(eta0, ns, 1), ds$y,
                   function(e, i) sc1(e, i), function(e, i) matrix(cv1(e, i)),
                   list(mu = psi))
  expect_equal(as.numeric(b$eta), as.numeric(a$eta), tolerance = 1e-12)
  expect_equal(unname(b$jacobian$mu), unname(a$jacobian), tolerance = 1e-12)
})

test_that("a subformula of one is the joint filter without it", {
  plain <- gas_joint(list(
    mu = term_build(gas(p = 1, q = 1, by = g, time = t, p_cross = c(sigma = 1)), ds),
    sigma = term_build(gas(p = 1, q = 1, by = g, time = t, p_cross = c(mu = 1)), ds)))
  dev <- gas_joint(list(
    mu = term_build(gas(p = 1, q = 1, kappa1 ~ 1, by = g, time = t,
                        p_cross = c(sigma = 1)), ds),
    sigma = term_build(gas(p = 1, q = 1, by = g, time = t, p_cross = c(mu = 1)), ds)))
  pp <- gs_psi(term_params(plain))
  pd <- gs_psi(term_params(dev))
  eta0 <- matrix(c(rep(1, ns), rep(0.3, ns)), ns, 2)
  a <- term_filter(plain, eta0, ds$y, gs_score(ds$y), gs_curv(ds$y), pp)
  b <- term_filter(dev, eta0, ds$y, gs_score(ds$y), gs_curv(ds$y), pd)
  expect_equal(b$eta, a$eta, tolerance = 1e-12)
  # the developed coordinate is log kappa: its column is kappa times the
  # scalar loading's column
  for (k in c("mu", "sigma")) {
    Ja <- a$jacobian[[k]]
    Jb <- b$jacobian[[k]]
    expect_equal(Jb[, "mu.kappa1.(Intercept)"], 0.2 * Ja[, "mu.kappa1"],
                 tolerance = 1e-12)
    others <- setdiff(colnames(Ja), "mu.kappa1")
    expect_equal(unname(Jb[, others]), unname(Ja[, others]), tolerance = 1e-12)
  }
})

test_that("the Jacobian of a developed joint filter is its derivative", {
  jt <- gas_joint(list(
    mu = term_build(gas(p = 1, q = 1, kappa1 ~ x, omega ~ g, by = g,
                        time = t, p_cross = c(sigma = 1)), ds),
    sigma = term_build(gas(p = 1, q = 1, by = g, time = t, p_cross = c(mu = 1)), ds)))
  nm <- term_params(jt)
  psi <- unlist(gs_psi(nm))
  eta0 <- matrix(c(rep(1, ns), rep(0.3, ns)), ns, 2)
  run <- function(v) term_filter(jt, eta0, ds$y, gs_score(ds$y),
                                 gs_curv(ds$y), as.list(v))
  base <- run(psi)
  for (j in seq_along(nm)) {
    h <- 1e-5
    up <- psi; up[j] <- up[j] + h
    dn <- psi; dn[j] <- dn[j] - h
    d <- (run(up)$eta - run(dn)$eta) / (2 * h)
    expect_equal(base$jacobian$mu[, j], d[, 1], tolerance = 1e-6,
                 label = paste("mu", nm[j]))
    expect_equal(base$jacobian$sigma[, j], d[, 2], tolerance = 1e-6,
                 label = paste("sigma", nm[j]))
  }
})

test_that("a cross loading takes a subformula only where p_cross declares it", {
  expect_error(gas(p = 1, q = 1, kappa1_mu ~ x), "does not\\s+declare")
  expect_error(gas(p = 1, q = 1, p_cross = c(mu = 1), kappa2_mu ~ x),
               "does not\\s+declare")
  expect_error(gas(p = 1, q = 1, p_cross = c(mu = 2), kappa2_mu ~ 0),
               "leaving its lag out")
  expect_silent(gas(p = 1, q = 1, p_cross = c(mu = 1), kappa1_mu ~ x,
                    links = list(kappa1_mu = linkfunctions7::identity_link())))
  mk <- function(pc, f = NULL) {
    if (is.null(f)) term_build(gas(p = 1, q = 1, by = g, time = t,
                                   p_cross = pc), ds)
    else term_build(gas(p = 1, q = 1, by = g, time = t, p_cross = pc, f), ds)
  }
  expect_error(gas_joint(list(mu = mk(c(nu = 1), kappa1_nu ~ x),
                              sigma = mk(NULL))), "has no gas\\(\\) term")
  jt <- gas_joint(list(mu = mk(c(sigma = 1), kappa1_sigma ~ x),
                       sigma = mk(NULL)))
  nm <- term_params(jt)
  expect_identical(nm, c("mu.omega", "mu.kappa1", "mu.pacf1",
                         "mu.kappa1_sigma.(Intercept)", "mu.kappa1_sigma.x",
                         "sigma.omega", "sigma.kappa1", "sigma.pacf1"))
  expect_identical(term_links(jt)[["mu.kappa1_sigma.x"]]@link_name, "identity")
})

test_that("a cross loading developed by a subformula of one is the scalar loading", {
  mk <- function(dev) list(
    mu = term_build(if (dev) gas(p = 1, q = 1, by = g, time = t, p_cross = c(sigma = 1),
                                 kappa1_sigma ~ 1)
                    else gas(p = 1, q = 1, by = g, time = t, p_cross = c(sigma = 1)), ds),
    sigma = term_build(gas(p = 1, q = 1, by = g, time = t, p_cross = c(mu = 1)), ds))
  plain <- gas_joint(mk(FALSE))
  dev <- gas_joint(mk(TRUE))
  pp <- gs_psi(term_params(plain))
  pd <- gs_psi(term_params(dev))
  # the identity link: the coordinate is the loading itself
  pd[["mu.kappa1_sigma.(Intercept)"]] <- pp[["mu.kappa1_sigma"]]
  eta0 <- matrix(c(rep(1, ns), rep(0.3, ns)), ns, 2)
  a <- term_filter(plain, eta0, ds$y, gs_score(ds$y), gs_curv(ds$y), pp)
  b <- term_filter(dev, eta0, ds$y, gs_score(ds$y), gs_curv(ds$y), pd)
  expect_equal(b$eta, a$eta, tolerance = 1e-12)
  for (k in c("mu", "sigma")) {
    Ja <- a$jacobian[[k]]
    Jb <- b$jacobian[[k]]
    colnames(Jb)[colnames(Jb) == "mu.kappa1_sigma.(Intercept)"] <- "mu.kappa1_sigma"
    expect_equal(unname(Jb[, colnames(Ja)]), unname(Ja), tolerance = 1e-12)
  }
})

test_that("the Jacobian of a filter with a developed cross loading is its derivative", {
  jt <- gas_joint(list(
    mu = term_build(gas(p = 1, q = 1, by = g, time = t, p_cross = c(sigma = 1),
                        kappa1_sigma ~ x), ds),
    sigma = term_build(gas(p = 1, q = 1, by = g, time = t, p_cross = c(mu = 1),
                           kappa1_mu ~ g), ds)))
  nm <- term_params(jt)
  psi <- unlist(gs_psi(nm))
  psi[grepl("_sigma.", nm, fixed = TRUE) | grepl("_mu.", nm, fixed = TRUE)] <-
    0.04
  eta0 <- matrix(c(rep(1, ns), rep(0.3, ns)), ns, 2)
  run <- function(v) term_filter(jt, eta0, ds$y, gs_score(ds$y),
                                 gs_curv(ds$y), as.list(v))
  base <- run(psi)
  for (j in seq_along(nm)) {
    h <- 1e-5
    up <- psi; up[j] <- up[j] + h
    dn <- psi; dn[j] <- dn[j] - h
    d <- (run(up)$eta - run(dn)$eta) / (2 * h)
    expect_equal(base$jacobian$mu[, j], d[, 1], tolerance = 1e-6,
                 label = paste("mu", nm[j]))
    expect_equal(base$jacobian$sigma[, j], d[, 2], tolerance = 1e-6,
                 label = paste("sigma", nm[j]))
  }
})

test_that("a one-way graph with a cross lag beyond p is the plain filter at a subformula of one", {
  # mu loads sigma's score at the lags 1..3 (its own p is 1), sigma loads
  # nothing: the subformula route against the route without subformulas
  mk <- function(f = NULL) list(
    mu = term_build(if (is.null(f)) gas(p = 1, q = 1, by = g, time = t,
                                        p_cross = c(sigma = 3))
                    else gas(p = 1, q = 1, by = g, time = t,
                             p_cross = c(sigma = 3), f), ds),
    sigma = term_build(gas(p = 1, q = 2, by = g, time = t), ds))
  plain <- gas_joint(mk())
  dev <- gas_joint(mk(kappa3_sigma ~ 1))
  pp <- gs_psi(term_params(plain))
  pp[["mu.kappa2_sigma"]] <- -0.03
  pp[["mu.kappa3_sigma"]] <- 0.02
  pd <- gs_psi(term_params(dev))
  pd[["mu.kappa2_sigma"]] <- -0.03
  pd[["mu.kappa3_sigma.(Intercept)"]] <- 0.02
  eta0 <- matrix(c(rep(1, ns), rep(0.3, ns)), ns, 2)
  a <- term_filter(plain, eta0, ds$y, gs_score(ds$y), gs_curv(ds$y), pp)
  b <- term_filter(dev, eta0, ds$y, gs_score(ds$y), gs_curv(ds$y), pd)
  expect_equal(b$eta, a$eta, tolerance = 1e-12)
  for (k in c("mu", "sigma")) {
    Jb <- b$jacobian[[k]]
    colnames(Jb)[colnames(Jb) == "mu.kappa3_sigma.(Intercept)"] <-
      "mu.kappa3_sigma"
    expect_equal(unname(Jb[, colnames(a$jacobian[[k]])]),
                 unname(a$jacobian[[k]]), tolerance = 1e-12)
  }
  gw <- matrix(stats::rnorm(2 * ns), ns, 2)
  ra <- term_adjoint(plain, eta0, ds$y, gs_score(ds$y), gs_curv(ds$y), pp, gw)
  rb <- term_adjoint(dev, eta0, ds$y, gs_score(ds$y), gs_curv(ds$y), pd, gw)
  expect_equal(rb$deta, ra$deta, tolerance = 1e-12)
  expect_equal(rb$dscore, ra$dscore, tolerance = 1e-12)
})

test_that("the Jacobian of a one-way graph with a developed cross lag beyond p is its derivative", {
  jt <- gas_joint(list(
    mu = term_build(gas(p = 1, q = 1, by = g, time = t,
                        p_cross = c(sigma = 3), kappa3_sigma ~ x), ds),
    sigma = term_build(gas(p = 1, q = 1, by = g, time = t), ds)))
  nm <- term_params(jt)
  expect_identical(nm, c("mu.omega", "mu.kappa1", "mu.pacf1",
                         "mu.kappa1_sigma", "mu.kappa2_sigma",
                         "mu.kappa3_sigma.(Intercept)", "mu.kappa3_sigma.x",
                         "sigma.omega", "sigma.kappa1", "sigma.pacf1"))
  psi <- unlist(gs_psi(nm))
  psi[grepl("_sigma", nm)] <- c(0.04, -0.03, 0.02, 0.01)
  eta0 <- matrix(c(rep(1, ns), rep(0.3, ns)), ns, 2)
  run <- function(v) term_filter(jt, eta0, ds$y, gs_score(ds$y),
                                 gs_curv(ds$y), as.list(v))
  base <- run(psi)
  for (j in seq_along(nm)) {
    h <- 1e-5
    up <- psi; up[j] <- up[j] + h
    dn <- psi; dn[j] <- dn[j] - h
    d <- (run(up)$eta - run(dn)$eta) / (2 * h)
    expect_equal(base$jacobian$mu[, j], d[, 1], tolerance = 1e-6,
                 label = paste("mu", nm[j]))
    expect_equal(base$jacobian$sigma[, j], d[, 2], tolerance = 1e-6,
                 label = paste("sigma", nm[j]))
  }
})
