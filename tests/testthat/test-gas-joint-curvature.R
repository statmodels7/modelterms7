# The second, third and fourth derivatives of a joint filter: against the
# scalar core at K = 1, and against one central difference of the analytic
# order below at K = 2.

# a gaussian in (mu, s = log sigma): l = -s - (y - mu)^2 exp(-2 s) / 2, whose
# derivative with j indices on mu and k on s is closed
gjc_dl <- function(y, e, idx) {
  j <- sum(idx == 1L)
  k <- sum(idx == 2L)
  r <- y - e[1]
  dj <- switch(min(j, 3L) + 1L, r^2, -2 * r, 2, 0)
  out <- -0.5 * dj * (-2)^k * exp(-2 * e[2])
  if (j == 0L && k == 1L) out <- out - 1
  out
}

# the driving data at predictors E (n x P): u and U1..U4 for the filtered
# parameters `fp`, with the scaling at zero so that U is the family's own
gjc_blocks <- function(y, E, fp, order = 4L) {
  n <- nrow(E)
  P <- ncol(E)
  K <- length(fp)
  u <- matrix(0, n, K)
  U1 <- array(0, c(n, K, P))
  U2 <- array(0, c(n, K, P, P))
  U3 <- array(0, c(n, K, P, P, P))
  U4 <- array(0, c(n, K, P, P, P, P))
  for (i in seq_len(n)) {
    for (l in seq_len(K)) {
      L <- fp[l]
      u[i, l] <- gjc_dl(y[i], E[i, ], L)
      for (a in seq_len(P)) {
        U1[i, l, a] <- gjc_dl(y[i], E[i, ], c(L, a))
        for (b in seq_len(P)) {
          U2[i, l, a, b] <- gjc_dl(y[i], E[i, ], c(L, a, b))
          if (order < 3L) next
          for (c in seq_len(P)) {
            U3[i, l, a, b, c] <- gjc_dl(y[i], E[i, ], c(L, a, b, c))
            if (order < 4L) next
            for (d in seq_len(P)) {
              U4[i, l, a, b, c, d] <- gjc_dl(y[i], E[i, ], c(L, a, b, c, d))
            }
          }
        }
      }
    }
  }
  list(filtered = fp, u = u, U1 = U1, U2 = U2, U3 = U3, U4 = U4)
}

set.seed(31)
nc <- 40
dc <- data.frame(t = seq_len(nc), y = rnorm(nc, 1, 1.4),
                 g = rep(c("a", "b"), each = nc / 2))
Xc <- cbind(1, as.numeric(scale(seq_len(nc))) * 0.3)

# a model of the gaussian with filters in the equations `fp` (1 = mu,
# 2 = log sigma), mu carrying the columns of Xc and log sigma an intercept;
# the unknowns are those coefficients followed by the filter's parameters on
# their unconstrained scale
gjc_model <- function(fp, p, q, cross = NULL, by = FALSE) {
  mk <- function() {
    if (by) gas(p = p, q = q, by = g, time = t, cross = cross)
    else gas(p = p, q = q, time = t, cross = cross)
  }
  nmp <- c("mu", "sigma")[fp]
  jt <- gas_joint(stats::setNames(lapply(nmp, function(k)
    term_build(mk(), dc)), nmp))
  M <- length(term_params(jt))
  mb <- ncol(Xc) + 1L
  m <- mb + M
  seed <- list(mu = cbind(Xc, 0, matrix(0, nc, M)),
               sigma = cbind(matrix(0, nc, ncol(Xc)), 1, matrix(0, nc, M)))
  links <- unlist(lapply(jt@params, function(k) term_links(jt@terms[[k]])),
                  recursive = FALSE)
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
  static <- function(u) cbind(as.numeric(Xc %*% u[1:2]), rep(u[3], nc))
  sc <- function(e, i) vapply(fp, function(L) gjc_dl(dc$y[i], e, L), 0)
  cu <- function(e, i) {
    matrix(vapply(fp, function(a) vapply(fp, function(L)
      gjc_dl(dc$y[i], e, c(L, a)), 0), numeric(length(fp))),
      length(fp), length(fp))
  }
  # the filter only sees its own equations; the others are static
  preds <- function(u) {
    E <- static(u)
    ef <- function(e, i) { full <- E[i, ]; full[fp] <- e; full }
    out <- term_filter(jt, E[, fp, drop = FALSE], dc$y,
                       function(e, i) sc(ef(e, i), i),
                       function(e, i) cu(ef(e, i), i), psi_of(u))
    E[, fp] <- out$eta
    E
  }
  core <- function(u, gw, dirs = NULL) {
    E <- preds(u)
    .gas_joint_core(jt, psi_of(u), gw, seed, gjc_blocks(dc$y, E, fp),
                    direction = dirs)
  }
  list(jt = jt, m = m, M = M, mb = mb, preds = preds, core = core)
}

gjc_u0 <- function(md, seed = 2) {
  set.seed(seed)
  z <- numeric(md$M)
  nm <- term_params(md$jt)
  z[grepl("omega$", nm)] <- stats::runif(sum(grepl("omega$", nm)), -0.1, 0.1)
  z[grepl("kappa[0-9]+$", nm)] <- log(stats::runif(sum(grepl("kappa[0-9]+$", nm)),
                                                   0.03, 0.08))
  z[grepl("_", nm)] <- stats::runif(sum(grepl("_", nm)), -0.04, 0.04)
  z[grepl("pacf", nm)] <- stats::runif(sum(grepl("pacf", nm)), 0.3, 0.8)
  c(0.9, 0.1, 0.25, z)
}

test_that("the joint core at K = 1 is the scalar core", {
  # the stub of test-gas.R, l = -log cosh(e - y), five non-zero derivatives
  set.seed(8)
  nn <- 50L
  d2 <- data.frame(t = seq_len(nn), y = stats::rnorm(nn))
  X <- cbind(1, as.numeric(scale(seq_len(nn))) * 0.3)
  mb <- ncol(X)
  tt <- function(e, i) tanh(e - d2$y[i])
  l1 <- function(e, i) -tt(e, i)
  l2 <- function(e, i) -(1 - tt(e, i)^2)
  l3 <- function(e, i) { t <- tt(e, i); 2 * t * (1 - t^2) }
  l4 <- function(e, i) { t <- tt(e, i); (1 - t^2) * (2 - 6 * t^2) }
  l5 <- function(e, i) { t <- tt(e, i); (1 - t^2) * (24 * t^3 - 16 * t) }
  for (cfg in list(c(1, 1), c(2, 2), c(1, 4))) {
    term <- term_build(gas(p = cfg[1], q = cfg[2], time = t), d2)
    jt <- gas_joint(list(mu = term))
    nm <- term_params(term)
    np <- length(nm)
    lk <- term_links(term)
    m <- mb + np
    seed <- cbind(X, matrix(0, nn, np))
    set.seed(9)
    u0 <- c(0.3, -0.2, stats::runif(np, -0.4, 0.2))
    psi <- as.list(stats::setNames(vapply(seq_along(nm), function(j)
      linkfunctions7::linkinv(lk[[nm[j]]], u0[mb + j]), numeric(1)), nm))
    gw <- stats::rnorm(nn)
    v <- stats::rnorm(m)
    w <- stats::rnorm(m)
    eta0 <- as.numeric(X %*% u0[seq_len(mb)])
    e <- term_filter(term, eta0, d2$y, l1, l2, psi)$eta
    ix <- seq_len(nn)
    bl <- list(filtered = 1L, u = matrix(l1(e, ix)),
               U1 = array(l2(e, ix), c(nn, 1, 1)),
               U2 = array(l3(e, ix), c(nn, 1, 1, 1)),
               U3 = array(l4(e, ix), c(nn, 1, 1, 1, 1)),
               U4 = array(l5(e, ix), c(nn, 1, 1, 1, 1, 1)))
    blk4 <- function(e, i, D, act) {
      DD <- outer(D, D)
      dv <- sum(D * v[act])
      dw <- sum(D * w[act])
      list(cross = numeric(length(D)), M = l3(e, i) * DD,
           dcurv = l3(e, i) * D, cppp = l3(e, i), Q = l4(e, i) * DD,
           N = list(l4(e, i) * dv * DD, l4(e, i) * dw * DD),
           P = l5(e, i) * dv * dw * DD)
    }
    blk3 <- function(e, i, D, act) list(
      cross = numeric(length(D)), M = l3(e, i) * outer(D, D),
      dcurv = l3(e, i) * D,
      N = l4(e, i) * sum(D * v[act]) * outer(D, D))
    blk2 <- function(e, i, D, act) list(cross = numeric(length(D)),
                                        M = l3(e, i) * outer(D, D))
    lab <- sprintf("p = %d, q = %d", cfg[1], cfg[2])
    s2 <- term_curvature(term, eta0, d2$y, l1, l2, psi, gw, seed, blk2)
    j2 <- term_curvature(jt, NULL, NULL, NULL, NULL, list(mu = psi),
                         matrix(gw), list(mu = seed), bl)
    expect_equal(j2$curvature, s2$curvature, tolerance = 1e-12, info = lab)
    expect_equal(j2$jacobian$mu, s2$jacobian, tolerance = 1e-13, info = lab)
    s3 <- term_third(term, eta0, d2$y, l1, l2, psi, gw, seed, blk3, v)
    j3 <- term_third(jt, NULL, NULL, NULL, NULL, list(mu = psi), matrix(gw),
                     list(mu = seed), bl, v)
    expect_equal(j3$curvature, s3$curvature, tolerance = 1e-12, info = lab)
    expect_equal(j3$dphi$mu[[1L]], s3$dphi, tolerance = 1e-12, info = lab)
    s4 <- term_fourth(term, eta0, d2$y, l1, l2, psi, gw, seed, blk4,
                      list(v, w))
    j4 <- term_fourth(jt, NULL, NULL, NULL, NULL, list(mu = psi), matrix(gw),
                      list(mu = seed), bl, list(v, w))
    expect_equal(j4$curvature, s4$curvature, tolerance = 1e-11, info = lab)
    expect_equal(j4$dpsi$mu, s4$dpsi, tolerance = 1e-12, info = lab)
  }
})

test_that("the joint core's Jacobian is the derivative of the predictors", {
  md <- gjc_model(1:2, 1, 1, cross = "a")
  u0 <- gjc_u0(md)
  got <- md$core(u0, matrix(0, nc, 2))
  h <- 1e-6
  for (j in seq_len(md$m)) {
    e <- numeric(md$m)
    e[j] <- h
    fd <- (md$preds(u0 + e) - md$preds(u0 - e)) / (2 * h)
    expect_equal(got$jacobian$mu[, j], fd[, 1], tolerance = 1e-7,
                 info = paste("mu", j))
    expect_equal(got$jacobian$sigma[, j], fd[, 2], tolerance = 1e-7,
                 info = paste("sigma", j))
  }
})

test_that("the joint curvature is the derivative of the joint Jacobian", {
  # one filtered equation beside a static one, both filtered, and a group of
  # cross loadings over a panel
  for (cfg in list(list(1L, 1, 1, NULL, FALSE), list(2L, 2, 1, NULL, FALSE),
                   list(1:2, 1, 1, NULL, FALSE), list(1:2, 2, 2, "a", TRUE))) {
    md <- gjc_model(cfg[[1]], cfg[[2]], cfg[[3]], cross = cfg[[4]],
                    by = cfg[[5]])
    K <- length(cfg[[1]])
    u0 <- gjc_u0(md)
    set.seed(5)
    gw <- matrix(stats::rnorm(nc * K), nc, K)
    got <- md$core(u0, gw)
    jac_of <- function(u) {
      J <- md$core(u, gw)$jacobian
      Reduce(`+`, lapply(seq_len(K), function(k) crossprod(gw[, k], J[[k]])))
    }
    h <- 1e-5
    fd <- t(vapply(seq_len(md$m), function(j) {
      e <- numeric(md$m)
      e[j] <- h
      as.numeric(jac_of(u0 + e) - jac_of(u0 - e)) / (2 * h)
    }, numeric(md$m)))
    lab <- paste(format(cfg), collapse = " ")
    expect_equal(got$curvature, fd, tolerance = 1e-7, info = lab)
  }
})

test_that("the joint third and fourth orders are the derivatives of the order below", {
  for (cfg in list(list(1L, 1, 1, NULL, FALSE), list(1:2, 1, 1, "a", FALSE),
                   list(1:2, 2, 2, "a", TRUE))) {
    md <- gjc_model(cfg[[1]], cfg[[2]], cfg[[3]], cross = cfg[[4]],
                    by = cfg[[5]])
    K <- length(cfg[[1]])
    u0 <- gjc_u0(md, 3)
    set.seed(6)
    gw <- matrix(stats::rnorm(nc * K), nc, K)
    v <- stats::rnorm(md$m)
    w <- stats::rnorm(md$m)
    lab <- paste(format(cfg), collapse = " ")
    dir_of <- function(f, x, h) (f(u0 + h * x) - f(u0 - h * x)) / (2 * h)
    rich <- function(f, x) (4 * dir_of(f, x, 1e-3) - dir_of(f, x, 2e-3)) / 3

    g3 <- md$core(u0, gw, v)
    expect_equal(g3$curvature, rich(function(u) md$core(u, gw)$curvature, v),
                 tolerance = 1e-6, info = lab)
    for (k in seq_len(K)) {
      expect_equal(g3$dphi[[k]][[1L]],
                   rich(function(u) md$core(u, gw)$jacobian[[k]], v),
                   tolerance = 1e-6, info = lab)
    }
    g4 <- md$core(u0, gw, list(v, w))
    expect_equal(g4$curvature,
                 rich(function(u) md$core(u, gw, v)$curvature, w),
                 tolerance = 1e-6, info = lab)
    for (k in seq_len(K)) {
      expect_equal(g4$dpsi[[k]],
                   rich(function(u) md$core(u, gw, v)$dphi[[k]][[1L]], w),
                   tolerance = 1e-6, info = lab)
    }
    # the lower orders come back as the lower calls give them
    expect_equal(g4$jacobian, g3$jacobian, tolerance = 1e-15)
    sw <- md$core(u0, gw, list(w, v))
    expect_equal(sw$curvature, g4$curvature, tolerance = 1e-11, info = lab)
  }
})

test_that("the compiled joint curvature is the R route's", {
  # gas_joint_curv_cpp() against the R recursion it replaces at the second
  # order, over the lags, cross loadings and a panel of groups
  for (cfg in list(list(fp = 1:2, p = 1, q = 1, cross = NULL, by = FALSE),
                   list(fp = 1:2, p = 2, q = 1, cross = "a", by = FALSE),
                   list(fp = 1:2, p = 1, q = 2, cross = "a", by = TRUE),
                   list(fp = 1, p = 1, q = 1, cross = NULL, by = TRUE))) {
    md <- gjc_model(cfg$fp, cfg$p, cfg$q, cross = cfg$cross, by = cfg$by)
    u0 <- gjc_u0(md)
    set.seed(5)
    gw <- matrix(stats::rnorm(nc * length(cfg$fp)), nc, length(cfg$fp))
    compiled <- md$core(u0, gw)
    old <- options(modelterms7.joint_core_r = TRUE)
    r <- md$core(u0, gw)
    options(old)
    lab <- sprintf("p = %d, q = %d, cross = %s, by = %s", cfg$p, cfg$q,
                   format(cfg$cross), cfg$by)
    expect_identical(names(compiled), names(r), info = lab)
    expect_equal(compiled$curvature, r$curvature, tolerance = 1e-12, info = lab)
    for (k in names(r$jacobian)) {
      expect_equal(compiled$jacobian[[k]], r$jacobian[[k]], tolerance = 1e-13,
                   info = paste(lab, k))
    }
  }
})
