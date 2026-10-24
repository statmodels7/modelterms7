# Several gas() terms filtered as one recursion: the joint filter, its
# Jacobian, its reverse pass, and what gas_joint() rejects.

set.seed(21)
nj <- 80
dj <- data.frame(t = seq_len(nj), y = rnorm(nj, 1, 1.5),
                 g = rep(c("a", "b"), each = nj / 2))

# a gaussian in (mu, log sigma), the driving vector being the score as it is
gj_score <- function(y) function(e, i) {
  s2 <- exp(2 * e[2])
  r <- y[i] - e[1]
  c(r / s2, r^2 / s2 - 1)
}
gj_curv <- function(y) function(e, i) {
  s2 <- exp(2 * e[2])
  r <- y[i] - e[1]
  matrix(c(-1 / s2, -2 * r / s2, -2 * r / s2, -2 * r^2 / s2), 2, 2)
}

gj_terms <- function(p, q, cross = NULL, by = FALSE, data = dj) {
  mk <- function() {
    if (by) gas(p = p, q = q, by = g, time = t, cross = cross)
    else gas(p = p, q = q, time = t, cross = cross)
  }
  list(mu = term_build(mk(), data), sigma = term_build(mk(), data))
}

gj_psi <- function(jt, seed = 1) {
  set.seed(seed)
  lapply(stats::setNames(jt@params, jt@params), function(k) {
    nm <- term_params(jt@terms[[k]])
    v <- stats::setNames(numeric(length(nm)), nm)
    for (j in nm) {
      v[[j]] <- if (j == "omega") stats::runif(1, -0.1, 0.1) else
        if (grepl("_", j)) stats::runif(1, -0.05, 0.05) else
        if (startsWith(j, "kappa")) stats::runif(1, 0.02, 0.08) else
        stats::runif(1, 0.3, 0.7)
    }
    as.list(v)
  })
}

test_that("one term filtered jointly is the scalar filter", {
  for (cfg in list(c(1, 1), c(2, 2), c(1, 0))) {
    for (by in c(FALSE, TRUE)) {
      tm <- if (by) term_build(gas(p = cfg[1], q = cfg[2], by = g, time = t),
                               dj)
            else term_build(gas(p = cfg[1], q = cfg[2], time = t), dj)
      jt <- gas_joint(list(mu = tm))
      psi <- gj_psi(jt)
      eta0 <- 0.2 + 0.1 * sin(seq_len(nj))
      sc <- function(e, i) dj$y[i] - e
      cu <- function(e, i) -1
      a <- term_filter(tm, eta0, dj$y, sc, cu, psi$mu)
      b <- term_filter(jt, matrix(eta0), dj$y, sc, cu, psi)
      lab <- sprintf("p = %d, q = %d, by = %s", cfg[1], cfg[2], by)
      expect_equal(as.numeric(b$eta), a$eta, tolerance = 1e-14, info = lab)
      expect_equal(unname(b$jacobian$mu), unname(a$jacobian),
                   tolerance = 1e-13, info = lab)
      gw <- stats::rnorm(nj)
      aa <- term_adjoint(tm, eta0, dj$y, sc, cu, psi$mu, gw)
      ab <- term_adjoint(jt, matrix(eta0), dj$y, sc, cu, psi, matrix(gw))
      expect_equal(as.numeric(ab$deta), aa$deta, tolerance = 1e-13, info = lab)
      expect_equal(as.numeric(ab$dscore), aa$dscore, tolerance = 1e-13,
                   info = lab)
    }
  }
})

test_that("the joint filter reproduces the recursion written out by hand", {
  jt <- gas_joint(gj_terms(1, 1, cross = "a"))
  psi <- gj_psi(jt, 3)
  eta0 <- cbind(rep(0.8, nj), rep(0.3, nj))
  out <- term_filter(jt, eta0, dj$y, gj_score(dj$y), gj_curv(dj$y), psi)
  pm <- psi$mu
  ps <- psi$sigma
  f <- c(pm$omega / (1 - pm$pacf1), ps$omega / (1 - ps$pacf1))
  sc <- gj_score(dj$y)
  ref <- matrix(0, nj, 2)
  for (t in seq_len(nj)) {
    if (t > 1L) {
      u <- sc(ref[t - 1L, ], t - 1L)
      f <- c(pm$omega + pm$kappa1 * u[1] + pm$kappa1_sigma * u[2] +
               pm$pacf1 * f[1],
             ps$omega + ps$kappa1 * u[2] + ps$kappa1_mu * u[1] +
               ps$pacf1 * f[2])
    }
    ref[t, ] <- eta0[t, ] + f
  }
  expect_equal(unname(out$eta), ref, tolerance = 1e-14)
})

test_that("the joint Jacobian is the derivative of the filter", {
  for (cfg in list(list(1, 1, NULL, FALSE), list(1, 1, "a", FALSE),
                   list(2, 2, "a", TRUE), list(2, 1, NULL, TRUE))) {
    jt <- gas_joint(gj_terms(cfg[[1]], cfg[[2]], cross = cfg[[3]],
                             by = cfg[[4]]))
    psi <- gj_psi(jt, 5)
    eta0 <- cbind(0.9 + 0.1 * cos(seq_len(nj)), rep(0.35, nj))
    run <- function(ps) term_filter(jt, eta0, dj$y, gj_score(dj$y),
                                    gj_curv(dj$y), ps)$eta
    got <- term_filter(jt, eta0, dj$y, gj_score(dj$y), gj_curv(dj$y), psi)
    nm <- term_params(jt)
    lab <- sprintf("p = %d, q = %d, cross = %s, by = %s", cfg[[1]], cfg[[2]],
                   format(cfg[[3]]), cfg[[4]])
    h <- 1e-6
    for (k in jt@params) {
      for (j in names(psi[[k]])) {
        up <- dn <- psi
        up[[k]][[j]] <- up[[k]][[j]] + h
        dn[[k]][[j]] <- dn[[k]][[j]] - h
        fd <- (run(up) - run(dn)) / (2 * h)
        col <- match(paste0(k, ".", j), nm)
        for (r in jt@params) {
          expect_equal(got$jacobian[[r]][, col], fd[, r], tolerance = 1e-7,
                       info = paste(lab, k, j, r))
        }
      }
    }
  }
})

test_that("the joint reverse pass is the derivative in the static predictors and the scores", {
  jt <- gas_joint(gj_terms(2, 1, cross = "a", by = TRUE))
  psi <- gj_psi(jt, 7)
  eta0 <- cbind(0.9 + 0.1 * cos(seq_len(nj)), rep(0.35, nj))
  set.seed(4)
  gw <- matrix(stats::rnorm(2 * nj), nj, 2)
  obj <- function(e0, sc = gj_score(dj$y)) {
    sum(gw * term_filter(jt, e0, dj$y, sc, gj_curv(dj$y), psi)$eta)
  }
  ad <- term_adjoint(jt, eta0, dj$y, gj_score(dj$y), gj_curv(dj$y), psi, gw)
  h <- 1e-6
  for (r in c(1, 2, 17, 40, 41, 66, 80)) {
    for (k in 1:2) {
      up <- dn <- eta0
      up[r, k] <- up[r, k] + h
      dn[r, k] <- dn[r, k] - h
      expect_equal(unname(ad$deta[r, k]), (obj(up) - obj(dn)) / (2 * h),
                   tolerance = 1e-7, info = sprintf("deta row %d, %d", r, k))
      # a perturbation of the driving score at one row, the curvature being
      # unchanged, isolates the derivative in that score
      bump <- function(dl) function(e, i) {
        s <- gj_score(dj$y)(e, i)
        if (i == r) s[k] <- s[k] + dl
        s
      }
      # the objective is of order ten, so the difference carries a rounding
      # of about 1e-9 against entries of order 1e-2
      expect_equal(unname(ad$dscore[r, k]),
                   (obj(eta0, bump(h)) - obj(eta0, bump(-h))) / (2 * h),
                   tolerance = 1e-6, info = sprintf("dscore row %d, %d", r, k))
    }
  }
})

test_that("cross labels add loadings to the terms of a group only", {
  tm <- list(mu = term_build(gas(p = 2, q = 1, time = t, cross = "a"), dj),
             sigma = term_build(gas(p = 1, q = 1, time = t, cross = "a"), dj),
             nu = term_build(gas(p = 1, q = 1, time = t), dj))
  jt <- gas_joint(tm)
  expect_identical(term_params(jt@terms$mu),
                   c("omega", "kappa1", "kappa2", "pacf1", "kappa1_sigma",
                     "kappa2_sigma"))
  expect_identical(term_params(jt@terms$sigma),
                   c("omega", "kappa1", "pacf1", "kappa1_mu"))
  expect_identical(term_params(jt@terms$nu), c("omega", "kappa1", "pacf1"))
  lk <- term_links(jt@terms$mu)
  expect_identical(lk$kappa1_sigma@link_name, "identity")
  expect_identical(term_start(jt@terms$mu)[["kappa1_sigma"]], 0)
  expect_identical(jt@offsets, c(0L, 6L, 10L))
})

test_that("gas_joint rejects terms that cannot be filtered together", {
  a <- term_build(gas(p = 1, q = 1, time = t), dj)
  expect_error(gas_joint(list(mu = a, sigma = term_build(
    gas(p = 1, q = 1, time = t, scaling = 1), dj))), "same 'scaling'")
  expect_error(gas_joint(list(
    mu = term_build(gas(p = 1, q = 1, time = t, scaling = 0.3), dj),
    sigma = term_build(gas(p = 1, q = 1, time = t, scaling = 0.3), dj))),
    "must be 0, 1/2 or 1")
  expect_silent(gas_joint(list(
    mu = term_build(gas(p = 1, q = 1, time = t, scaling = 0.3), dj))))
  expect_error(gas_joint(list(mu = a, sigma = term_build(
    gas(p = 1, q = 1, time = t, cross = "a"), dj))), "one gas\\(\\) term only")
  expect_error(gas_joint(list(mu = a, sigma = term_build(
    gas(p = 1, q = 1, by = g, time = t), dj))), "walk different series")
  expect_error(gas_joint(list(mu = a, mu = a)), "distinct parameters")
  expect_error(gas(cross = c("a", "b")), "single non-empty")
})

test_that("the joint continuation is the filter with the scores at zero", {
  n0 <- 50
  nh <- 6
  dfull <- data.frame(t = seq_len(n0 + nh), y = c(dj$y[1:n0], rep(0, nh)),
                      g = "a")
  for (cfg in list(list(1, 1, NULL), list(2, 2, "a"))) {
    mk <- function(dat) gas_joint(list(
      mu = term_build(gas(p = cfg[[1]], q = cfg[[2]], time = t,
                          cross = cfg[[3]]), dat),
      sigma = term_build(gas(p = cfg[[1]], q = cfg[[2]], time = t,
                             cross = cfg[[3]]), dat)))
    jt <- mk(dfull[1:n0, ])
    jf <- mk(dfull)
    psi <- gj_psi(jt, 11)
    e0 <- cbind(rep(0.9, n0 + nh), rep(0.35, n0 + nh))
    sc <- function(e, i) if (i > n0) c(0, 0) else gj_score(dfull$y)(e, i)
    cu <- function(e, i) if (i > n0) matrix(0, 2, 2) else
      gj_curv(dfull$y)(e, i)
    full <- term_filter(jf, e0, dfull$y, sc, cu, psi)
    obs <- term_filter(jt, e0[1:n0, ], dfull$y[1:n0], sc, cu, psi)
    f_past <- obs$eta - e0[1:n0, ]
    s_past <- t(vapply(seq_len(n0), function(i) sc(obs$eta[i, ], i),
                       numeric(2)))
    got <- term_continue(jt, psi, f_past, s_past, dfull[n0 + seq_len(nh), ])
    lab <- paste(format(cfg), collapse = " ")
    expect_equal(unname(got), unname(full$eta[n0 + seq_len(nh), ] -
                                       e0[n0 + seq_len(nh), ]),
                 tolerance = 1e-13, info = lab)
    # the derivative of the continuation in the joint vector: the past
    # levels' derivative is the filter's Jacobian, the past scores' is the
    # curvature times it
    M <- length(term_params(jt))
    dfp <- lapply(1:2, function(k) obs$jacobian[[k]])
    dsp <- lapply(1:2, function(l) {
      out <- matrix(0, n0, M)
      for (i in seq_len(n0)) {
        C <- matrix(obs$curv[i, , ], 2, 2)
        for (k in 1:2) out[i, ] <- out[i, ] + C[l, k] * obs$jacobian[[k]][i, ]
      }
      out
    })
    gd <- term_continue(jt, psi, f_past, s_past, dfull[n0 + seq_len(nh), ],
                        deriv = list(df_past = dfp, ds_past = dsp,
                                     dpsi = diag(M)))
    nm <- term_params(jt)
    h <- 1e-6
    for (k in jt@params) for (j in names(psi[[k]])) {
      up <- dn <- psi
      up[[k]][[j]] <- up[[k]][[j]] + h
      dn[[k]][[j]] <- dn[[k]][[j]] - h
      fd <- (term_filter(jf, e0, dfull$y, sc, cu, up)$eta -
               term_filter(jf, e0, dfull$y, sc, cu, dn)$eta) / (2 * h)
      col <- match(paste0(k, ".", j), nm)
      for (r in 1:2) {
        expect_equal(attr(gd, "jacobian")[[r]][, col],
                     fd[n0 + seq_len(nh), r], tolerance = 1e-7,
                     info = paste(lab, k, j, r))
      }
    }
  }
})

test_that("a joint filter is a structural term with prefixed parameters", {
  jt <- gas_joint(gj_terms(1, 2, cross = "a"))
  expect_true(S7::S7_inherits(jt, structural_term))
  nm <- term_params(jt)
  expect_identical(nm, c("mu.omega", "mu.kappa1", "mu.pacf1", "mu.pacf2",
                         "mu.kappa1_sigma", "sigma.omega", "sigma.kappa1",
                         "sigma.pacf1", "sigma.pacf2", "sigma.kappa1_mu"))
  expect_identical(names(term_links(jt)), nm)
  expect_identical(names(term_start(jt)), nm)
  expect_identical(unname(term_start(jt)[c("mu.kappa1", "mu.kappa1_sigma")]),
                   c(log(0.1), 0))
  expect_identical(term_level_param(jt), c("mu.omega", "sigma.omega"))
  expect_identical(term_penalties(jt), list())
  z <- as.list(stats::setNames(seq_along(nm) / 20, nm))
  r <- term_readable(jt, z)
  rm <- term_readable(jt@terms$mu, as.list(stats::setNames(
    unlist(z[1:5]), term_params(jt@terms$mu))))
  expect_identical(r$name[1:5], paste0("mu.", rm$name))
  expect_identical(r$value[1:5], rm$value)
  expect_identical(unname(r$jacobian[1:5, 1:5]), unname(rm$jacobian))
  expect_true(all(r$jacobian[1:5, 6:10] == 0))
  # building again keeps the groups and the label
  jb <- term_build(jt, dj[1:60, ])
  expect_identical(jb@blueprint$n, 60L)
  expect_identical(term_params(jb), nm)
  expect_output(print(jt), "mu: p = 1, q = 2, loads sigma")
})

test_that("a developed member keeps its layout and its cross loadings in the joint vector", {
  jt <- gas_joint(list(
    mu = term_build(gas(p = 1, q = 1, omega ~ g, time = t, cross = "a"), dj),
    sigma = term_build(gas(p = 1, q = 1, time = t, cross = "a"), dj)))
  nm <- term_params(jt)
  expect_identical(nm[seq_len(5L)], c("mu.omega.(Intercept)", "mu.omega.gb",
                                      "mu.kappa1", "mu.pacf1",
                                      "mu.kappa1_sigma"))
  expect_identical(nm[-seq_len(5L)], c("sigma.omega", "sigma.kappa1",
                                       "sigma.pacf1", "sigma.kappa1_mu"))
  cps <- term_components(jt)
  expect_identical(cps[["mu.omega"]]$index, 1:2)
  expect_identical(cps[["mu.kappa1_sigma"]]$index, 5L)
  expect_identical(cps[["sigma.kappa1_mu"]]$index, 9L)
})
