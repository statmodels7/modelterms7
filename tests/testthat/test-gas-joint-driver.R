# A driver: a parameter without a gas() term whose score the levels of the
# other terms load through p_cross. It is a score-driven term with no
# parameters and a level of zero, so it must equal a gas(p = 1, q = 0) term
# on the same parameter whose omega and kappa1 (identity link) are zero.

gd_dl <- function(y, e, idx) {
  j <- sum(idx == 1L)
  k <- sum(idx == 2L)
  r <- y - e[1]
  dj <- switch(min(j, 3L) + 1L, r^2, -2 * r, 2, 0)
  out <- -0.5 * dj * (-2)^k * exp(-2 * e[2])
  if (j == 0L && k == 1L) out <- out - 1
  out
}
gd_blocks <- function(y, E) {
  n <- nrow(E)
  u <- matrix(0, n, 2)
  U1 <- array(0, c(n, 2, 2))
  U2 <- array(0, c(n, 2, 2, 2))
  U3 <- array(0, c(n, 2, 2, 2, 2))
  U4 <- array(0, c(n, 2, 2, 2, 2, 2))
  for (i in seq_len(n)) for (l in 1:2) {
    u[i, l] <- gd_dl(y[i], E[i, ], l)
    for (a in 1:2) {
      U1[i, l, a] <- gd_dl(y[i], E[i, ], c(l, a))
      for (b in 1:2) {
        U2[i, l, a, b] <- gd_dl(y[i], E[i, ], c(l, a, b))
        for (c in 1:2) {
          U3[i, l, a, b, c] <- gd_dl(y[i], E[i, ], c(l, a, b, c))
          for (d in 1:2)
            U4[i, l, a, b, c, d] <- gd_dl(y[i], E[i, ], c(l, a, b, c, d))
        }
      }
    }
  }
  list(filtered = 1:2, u = u, U1 = U1, U2 = U2, U3 = U3, U4 = U4)
}

set.seed(31)
nd <- 30
dd <- data.frame(t = seq_len(nd), y = rnorm(nd, 1, 1.4),
                 g = rep(c("a", "b"), each = nd / 2), x = sin(seq_len(nd) / 3))
Xd <- cbind(1, as.numeric(scale(seq_len(nd))) * 0.3)

# mu's term from mk_mu; with phantom, sigma carries gas(p = 1, q = 0) whose
# coordinates are held at zero by the caller, otherwise sigma is a driver
gd_model <- function(mk_mu, phantom) {
  mu <- term_build(mk_mu(), dd)
  jt <- if (phantom) {
    gas_joint(list(mu = mu, sigma = term_build(gas(
      p = 1, q = 0, by = g, time = t,
      links = list(kappa1 = linkfunctions7::identity_link())), dd)))
  } else {
    gas_joint(list(mu = mu), drivers = "sigma")
  }
  M <- length(term_params(jt))
  mb <- 3L
  seed <- list(mu = cbind(Xd, 0, matrix(0, nd, M)),
               sigma = cbind(matrix(0, nd, 2), 1, matrix(0, nd, M)))
  psi_of <- function(u) {
    z <- u[mb + seq_len(M)]
    out <- list()
    at <- 0L
    for (k in jt@params) {
      nm <- term_params(jt@terms[[k]])
      lk <- term_links(jt@terms[[k]])
      out[[k]] <- as.list(stats::setNames(vapply(seq_along(nm), function(j)
        linkfunctions7::linkinv(lk[[nm[j]]], z[at + j]), numeric(1)), nm))
      at <- at + length(nm)
    }
    out
  }
  static <- function(u) cbind(as.numeric(Xd %*% u[1:2]), rep(u[3], nd))
  sc <- function(e, i) vapply(1:2, function(l) gd_dl(dd$y[i], e, l), 0)
  cu <- function(e, i) matrix(vapply(1:2, function(a) vapply(1:2, function(l)
    gd_dl(dd$y[i], e, c(l, a)), 0), numeric(2)), 2, 2)
  filt <- function(u) term_filter(jt, static(u), dd$y, sc, cu, psi_of(u))
  core <- function(u, gw, dirs = NULL) {
    E <- static(u)
    E[, 1:2] <- filt(u)$eta
    .gas_joint_core(jt, psi_of(u), gw, seed, gd_blocks(dd$y, E),
                    direction = dirs)
  }
  list(jt = jt, m = mb + M, M = M, mb = mb, filt = filt, core = core,
       static = static, sc = sc, cu = cu, psi_of = psi_of)
}

gd_cases <- list(
  plain = function() gas(p = 1, q = 1, by = g, time = t,
                         p_cross = c(sigma = 2)),
  developed = function() gas(p = 1, q = 1, kappa1 ~ x, by = g, time = t,
                             p_cross = c(sigma = 2)),
  developed_cross = function() gas(p = 1, q = 1, by = g, time = t,
                                   p_cross = c(sigma = 2), kappa2_sigma ~ x))

test_that("a driver is a gas() term with its coordinates at zero, at every order", {
  for (cn in names(gd_cases)) {
    A <- gd_model(gd_cases[[cn]], FALSE)
    B <- gd_model(gd_cases[[cn]], TRUE)
    nmA <- term_params(A$jt)
    common <- match(nmA, term_params(B$jt))
    expect_false(anyNA(common), info = cn)
    set.seed(4)
    zA <- c(0.9, 0.1, 0.25, stats::runif(A$M, -0.05, 0.05))
    zA[A$mb + grep("kappa1$|kappa1[.][(]Intercept[)]", nmA)] <- log(0.05)
    zA[A$mb + grep("pacf", nmA)] <- 0.6
    lift <- function(x) {
      out <- numeric(B$m)
      out[1:3] <- x[1:3]
      out[B$mb + common] <- x[A$mb + seq_len(A$M)]
      out
    }
    set.seed(7)
    gw <- matrix(stats::rnorm(nd * 2), nd, 2)
    v <- stats::rnorm(A$m)
    w <- stats::rnorm(A$m)
    keep <- c(1:3, B$mb + common)
    for (route in c(FALSE, TRUE)) {
      old <- options(modelterms7.joint_core_r = route)
      lab <- sprintf("%s, R route %s", cn, route)
      fA <- A$filt(zA)
      fB <- B$filt(lift(zA))
      expect_equal(fA$eta, fB$eta, tolerance = 1e-12, info = lab)
      expect_equal(unname(fA$jacobian$mu),
                   unname(fB$jacobian$mu[, common]), tolerance = 1e-12,
                   info = lab)
      expect_true(all(fA$jacobian$sigma == 0), info = lab)
      a2 <- A$core(zA, gw)
      b2 <- B$core(lift(zA), gw)
      expect_equal(a2$curvature, b2$curvature[keep, keep], tolerance = 1e-12,
                   info = lab)
      a3 <- A$core(zA, gw, v)
      b3 <- B$core(lift(zA), gw, lift(v))
      expect_equal(a3$curvature, b3$curvature[keep, keep], tolerance = 1e-12,
                   info = lab)
      expect_equal(unname(a3$dphi$mu[[1L]]),
                   unname(b3$dphi$mu[[1L]][, keep]), tolerance = 1e-12,
                   info = lab)
      a4 <- A$core(zA, gw, list(v, w))
      b4 <- B$core(lift(zA), gw, list(lift(v), lift(w)))
      expect_equal(a4$curvature, b4$curvature[keep, keep], tolerance = 1e-12,
                   info = lab)
      expect_equal(unname(a4$dpsi$mu), unname(b4$dpsi$mu[, keep]),
                   tolerance = 1e-12, info = lab)
      options(old)
    }
    # the reverse pass reaches the driver's static predictor through its
    # score, as the phantom's does
    eta0 <- A$static(zA)
    ga <- term_adjoint(A$jt, eta0, dd$y, A$sc, A$cu, A$psi_of(zA), gw)
    gb <- term_adjoint(B$jt, eta0, dd$y, B$sc, B$cu, B$psi_of(lift(zA)), gw)
    expect_equal(ga$deta, gb$deta, tolerance = 1e-12, info = cn)
    expect_equal(ga$dscore, gb$dscore, tolerance = 1e-12, info = cn)
  }
})

test_that("a driver's continuation is the phantom's", {
  n0 <- 24
  past <- dd[seq_len(n0), ]
  new <- dd[n0 + 1:4, ]
  A <- gas_joint(list(mu = term_build(gas(p = 1, q = 1, time = t,
                                          p_cross = c(sigma = 3)), past)),
                 drivers = "sigma")
  B <- gas_joint(list(mu = term_build(gas(p = 1, q = 1, time = t,
                                          p_cross = c(sigma = 3)), past),
                      sigma = term_build(gas(p = 1, q = 0, time = t,
                        links = list(kappa1 = linkfunctions7::identity_link())),
                        past)))
  pm <- list(omega = 0.1, kappa1 = 0.2, pacf1 = 0.5, kappa1_sigma = 0.04,
             kappa2_sigma = -0.03, kappa3_sigma = 0.02)
  set.seed(9)
  fp <- cbind(stats::rnorm(n0, 0, 0.1), 0)
  sp <- cbind(stats::rnorm(n0), stats::rnorm(n0))
  a <- term_continue(A, list(mu = pm), fp, sp, new)
  b <- term_continue(B, list(mu = pm, sigma = list(omega = 0, kappa1 = 0)),
                     fp, sp, new)
  expect_equal(unname(a), unname(b), tolerance = 1e-14)
  expect_true(all(a[, "sigma"] == 0))
})

test_that("a driver has no parameters and is printed as static", {
  jt <- gas_joint(list(mu = term_build(gas(p = 1, q = 1, time = t,
                                           p_cross = c(sigma = 1)), dd)),
                  drivers = "sigma")
  expect_identical(jt@params, c("mu", "sigma"))
  expect_true(S7::S7_inherits(jt@terms$sigma, GasDriver))
  expect_identical(term_params(jt), c("mu.omega", "mu.kappa1", "mu.pacf1",
                                      "mu.kappa1_sigma"))
  expect_identical(names(term_links(jt)), term_params(jt))
  expect_identical(names(term_start(jt)), term_params(jt))
  expect_identical(term_level_param(jt), "mu.omega")
  expect_identical(names(term_components(jt)),
                   c("mu.omega", "mu.kappa1", "mu.pacf1", "mu.kappa1_sigma"))
  r <- term_readable(jt, as.list(stats::setNames(rep(0.1, 4),
                                                 term_params(jt))))
  expect_identical(r$name, c("mu.omega", "mu.kappa1", "mu.xi1",
                             "mu.kappa1_sigma"))
  expect_output(print(jt), "sigma: static, its score drives the filter")
  expect_output(print(jt), "1 filtered parameters")
  # building again on other rows keeps the driver
  jb <- term_build(jt, dd[1:20, ])
  expect_identical(jb@params, c("mu", "sigma"))
  expect_true(S7::S7_inherits(jb@terms$sigma, GasDriver))
  expect_identical(jb@blueprint$n, 20L)
})

test_that("gas_joint rejects drivers that it cannot use", {
  a <- term_build(gas(p = 1, q = 1, time = t, p_cross = c(sigma = 1)), dd)
  s <- term_build(gas(p = 1, q = 1, time = t), dd)
  expect_error(gas_joint(list(mu = a, sigma = s), drivers = "sigma"),
               "carries a gas\\(\\) term and is named in 'drivers'")
  expect_error(gas_joint(list(mu = a)), "is not among 'drivers'")
  expect_error(gas_joint(list(mu = a), drivers = c("sigma", "sigma")),
               "distinct parameter names")
  expect_error(gas_joint(list(mu = term_build(gas(p = 1, q = 1, time = t,
                                                  scaling = 0.3,
                                                  p_cross = c(sigma = 1)),
                                              dd)), drivers = "sigma"),
               "must be 0, 1/2 or 1")
})
