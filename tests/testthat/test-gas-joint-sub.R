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
    mu = term_build(gas(p = 1, q = 1, by = g, time = t, cross = "a"), ds),
    sigma = term_build(gas(p = 1, q = 1, by = g, time = t, cross = "a"), ds)))
  dev <- gas_joint(list(
    mu = term_build(gas(p = 1, q = 1, kappa1 ~ 1, by = g, time = t,
                        cross = "a"), ds),
    sigma = term_build(gas(p = 1, q = 1, by = g, time = t, cross = "a"), ds)))
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
                        time = t, cross = "a"), ds),
    sigma = term_build(gas(p = 1, q = 1, by = g, time = t, cross = "a"), ds)))
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
