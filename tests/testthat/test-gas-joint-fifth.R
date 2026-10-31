# The fifth order of a joint filter: at K = 1 it is the scalar term's, and
# with two members (cross loadings, groups, subformulas) it is a Richardson
# difference of the joint fourth order along the third direction.

test_that("the joint fifth order at K = 1 is the scalar term's", {
  set.seed(8)
  nn <- 40L
  d2 <- data.frame(t = seq_len(nn), y = stats::rnorm(nn),
                   x = sin(seq_len(nn) / 4))
  X <- cbind(1, as.numeric(scale(seq_len(nn))) * 0.3)
  mb <- ncol(X)
  tt <- function(e, i) tanh(e - d2$y[i])
  l1 <- function(e, i) -tt(e, i)
  l2 <- function(e, i) -(1 - tt(e, i)^2)
  ld <- list(
    function(t) -(1 - t^2), function(t) 2 * t * (1 - t^2),
    function(t) (1 - t^2) * (2 - 6 * t^2),
    function(t) (1 - t^2) * (24 * t^3 - 16 * t),
    function(t) (1 - t^2) * (-2 * t * (24 * t^3 - 16 * t) +
                               (1 - t^2) * (72 * t^2 - 16)))
  specs <- list(quote(gas(p = 1, q = 2, time = t)),
                quote(gas(p = 1, q = 1, kappa1 ~ x, time = t)))
  for (sp in specs) {
    term <- term_build(eval(sp), d2)
    jt <- gas_joint(list(mu = term))
    nm <- term_params(term)
    np <- length(nm)
    lk <- term_links(term)
    base <- .gas_base_params(term@p, term@q)
    m <- mb + np
    seed <- cbind(X, matrix(0, nn, np))
    set.seed(9)
    u0 <- c(0.3, -0.2, stats::runif(np, -0.4, 0.2))
    if ("kappa1.(Intercept)" %in% nm) {
      u0[mb + which(nm == "kappa1.(Intercept)")] <- log(0.1)
    }
    psi <- as.list(stats::setNames(vapply(seq_along(nm), function(j)
      if (nm[j] %in% base) linkfunctions7::linkinv(lk[[nm[j]]], u0[mb + j])
      else u0[mb + j], numeric(1)), nm))
    gw <- stats::rnorm(nn)
    dirs <- lapply(1:3, function(i) stats::rnorm(m))
    eta0 <- as.numeric(X %*% u0[seq_len(mb)])
    e <- term_filter(term, eta0, d2$y, l1, l2, psi)$eta
    tz <- tanh(e - d2$y)
    arr <- list(H = matrix(ld[[1]](tz)), D3 = matrix(ld[[2]](tz)),
                D4 = matrix(ld[[3]](tz)), D5 = matrix(ld[[4]](tz)),
                D6 = matrix(ld[[5]](tz)), Vs = list(matrix(0, nn, m)),
                ap = 1L)
    bl <- list(filtered = 1L, u = matrix(-tz),
               U1 = array(ld[[1]](tz), c(nn, 1, 1)),
               U2 = array(ld[[2]](tz), c(nn, 1, 1, 1)),
               U3 = array(ld[[3]](tz), c(nn, 1, 1, 1, 1)),
               U4 = array(ld[[4]](tz), c(nn, 1, 1, 1, 1, 1)),
               U5 = array(ld[[5]](tz), c(nn, 1, 1, 1, 1, 1, 1)))
    s5 <- term_fifth(term, eta0, d2$y, l1, psi, gw, seed, arr, dirs)
    j5 <- term_fifth(jt, NULL, NULL, NULL, list(mu = psi), matrix(gw),
                     list(mu = seed), bl, dirs)
    lab <- deparse(sp)
    expect_equal(j5$curvature, s5$curvature, tolerance = 1e-12, info = lab)
    expect_equal(j5$domega$mu, s5$domega, tolerance = 1e-12, info = lab)
    expect_equal(j5$dpsi$mu, s5$dpsi, tolerance = 1e-12, info = lab)
  }
})

# The two-member model of test-gas-joint-curvature.R: a gaussian with mean
# and log-scale filtered, the derivatives of the driving scores in closed
# form to any order.
gjf_dl <- function(y, e, idx) {
  j <- sum(idx == 1L)
  k <- sum(idx == 2L)
  r <- y - e[1]
  dj <- switch(min(j, 3L) + 1L, r^2, -2 * r, 2, 0)
  out <- -0.5 * dj * (-2)^k * exp(-2 * e[2])
  if (j == 0L && k == 1L) out <- out - 1
  out
}

gjf_blocks <- function(y, E, fp, order = 5L) {
  n <- nrow(E)
  P <- ncol(E)
  K <- length(fp)
  out <- list(filtered = fp, u = matrix(0, n, K))
  for (o in seq_len(order)) {
    idx <- as.matrix(do.call(expand.grid, rep(list(seq_len(P)), o)))
    A <- array(0, c(n, K, rep(P, o)))
    for (i in seq_len(n)) {
      for (l in seq_len(K)) {
        vals <- apply(idx, 1L, function(ix) gjf_dl(y[i], E[i, ], c(fp[l], ix)))
        sl <- c(list(i, l), rep(list(seq_len(P)), o))
        A <- do.call(`[<-`, c(list(A), sl, list(value = array(vals, rep(P, o)))))
      }
    }
    out[[paste0("U", o)]] <- A
  }
  for (i in seq_len(n)) {
    for (l in seq_len(K)) out$u[i, l] <- gjf_dl(y[i], E[i, ], fp[l])
  }
  out
}

set.seed(31)
nf <- 30
df5 <- data.frame(t = seq_len(nf), y = rnorm(nf, 1, 1.4),
                  g = rep(c("a", "b"), each = nf / 2))
df5$x <- sin(seq_len(nf) / 3)
Xf <- cbind(1, as.numeric(scale(seq_len(nf))) * 0.3)

gjf_model <- function(fp, p, q, cross = NULL, by = FALSE, dev = FALSE,
                      devx = NULL) {
  oth <- c(mu = "sigma", sigma = "mu")
  mk <- function(k) {
    cross <- if (identical(cross, "sym")) stats::setNames(p, oth[[k]])
             else if (identical(cross, "upper") && k == "mu") c(sigma = p + 2)
    if (!is.null(devx) && k == "mu") {
      return(gas(p = p, q = q, devx, by = g, time = t, p_cross = cross))
    }
    if (dev && k == "mu") {
      if (by) gas(p = p, q = q, kappa1 ~ x, by = g, time = t, p_cross = cross)
      else gas(p = p, q = q, kappa1 ~ x, time = t, p_cross = cross)
    } else if (by) gas(p = p, q = q, by = g, time = t, p_cross = cross)
    else gas(p = p, q = q, time = t, p_cross = cross)
  }
  nmp <- c("mu", "sigma")[fp]
  jt <- gas_joint(stats::setNames(lapply(nmp, function(k)
    term_build(mk(k), df5)), nmp))
  M <- length(term_params(jt))
  mb <- ncol(Xf) + 1L
  m <- mb + M
  seed <- list(mu = cbind(Xf, 0, matrix(0, nf, M)),
               sigma = cbind(matrix(0, nf, ncol(Xf)), 1, matrix(0, nf, M)))
  psi_of <- function(u) {
    z <- u[mb + seq_len(M)]
    out <- list()
    at <- 0L
    for (k in jt@params) {
      tm <- jt@terms[[k]]
      nm <- term_params(tm)
      lay <- .gas_sub_layout(tm)
      vals <- numeric(length(nm))
      for (j in names(lay$idx)) {
        ix <- lay$idx[[j]]
        vals[ix] <- if (is.null(tm@build_info$sub[[j]])) {
          linkfunctions7::linkinv(.gas_param_link(tm, j), z[at + ix])
        } else z[at + ix]
      }
      out[[k]] <- as.list(stats::setNames(vals, nm))
      at <- at + length(nm)
    }
    out
  }
  static <- function(u) cbind(as.numeric(Xf %*% u[1:2]), rep(u[3], nf))
  sc <- function(e, i) vapply(fp, function(L) gjf_dl(df5$y[i], e, L), 0)
  cu <- function(e, i) {
    matrix(vapply(fp, function(a) vapply(fp, function(L)
      gjf_dl(df5$y[i], e, c(L, a)), 0), numeric(length(fp))),
      length(fp), length(fp))
  }
  preds <- function(u) {
    E <- static(u)
    ef <- function(e, i) { full <- E[i, ]; full[fp] <- e; full }
    out <- term_filter(jt, E[, fp, drop = FALSE], df5$y,
                       function(e, i) sc(ef(e, i), i),
                       function(e, i) cu(ef(e, i), i), psi_of(u))
    E[, fp] <- out$eta
    E
  }
  core <- function(u, gw, dirs = NULL) {
    E <- preds(u)
    .gas_joint_core(jt, psi_of(u), gw, seed, gjf_blocks(df5$y, E, fp, 4L),
                    direction = dirs, eta = static(u)[, fp, drop = FALSE])
  }
  fifth <- function(u, gw, dirs, blocks = NULL) {
    E <- preds(u)
    if (is.null(blocks)) blocks <- gjf_blocks(df5$y, E, fp)
    term_fifth(jt, static(u)[, fp, drop = FALSE], df5$y, NULL, psi_of(u),
               gw, seed, blocks, dirs)
  }
  list(jt = jt, m = m, M = M, mb = mb, preds = preds, core = core,
       fifth = fifth)
}

gjf_u0 <- function(md, seed = 2) {
  set.seed(seed)
  z <- numeric(md$M)
  nm <- term_params(md$jt)
  z[grepl("omega$", nm)] <- stats::runif(sum(grepl("omega$", nm)), -0.1, 0.1)
  z[grepl("kappa[0-9]+$", nm)] <- log(stats::runif(sum(grepl("kappa[0-9]+$", nm)),
                                                   0.03, 0.08))
  z[grepl("_", nm)] <- stats::runif(sum(grepl("_", nm)), -0.04, 0.04)
  z[grepl("pacf", nm)] <- stats::runif(sum(grepl("pacf", nm)), 0.3, 0.8)
  z[endsWith(nm, "kappa1.(Intercept)")] <- log(0.05)
  z[endsWith(nm, "kappa1.x")] <- 0.2
  c(0.9, 0.1, 0.25, z)
}

test_that("the joint fifth order is the derivative of the joint fourth", {
  cfgs <- list(list(1:2, 1, 1, "sym", FALSE, FALSE, NULL),
               list(1:2, 2, 1, "upper", TRUE, FALSE, NULL),
               list(1:2, 1, 2, "sym", TRUE, TRUE, NULL),
               list(1:2, 1, 1, "sym", TRUE, FALSE, kappa1_sigma ~ x))
  for (cfg in cfgs) {
    md <- gjf_model(cfg[[1]], cfg[[2]], cfg[[3]], cross = cfg[[4]],
                    by = cfg[[5]], dev = cfg[[6]], devx = cfg[[7]])
    K <- length(cfg[[1]])
    u0 <- gjf_u0(md, 3)
    set.seed(6)
    gw <- matrix(stats::rnorm(nf * K), nf, K)
    v <- stats::rnorm(md$m)
    w <- stats::rnorm(md$m)
    x <- stats::rnorm(md$m)
    lab <- paste(format(cfg), collapse = " ")
    dir_of <- function(f, h) (f(u0 + h * x) - f(u0 - h * x)) / (2 * h)
    rich <- function(f) (4 * dir_of(f, 1e-3) - dir_of(f, 2e-3)) / 3
    g5 <- md$fifth(u0, gw, list(v, w, x))
    ref <- rich(function(u) md$core(u, gw, list(v, w))$curvature)
    expect_equal(g5$curvature, ref, tolerance = 1e-6, info = lab)
    for (k in seq_len(K)) {
      expect_equal(g5$domega[[k]],
                   rich(function(u) md$core(u, gw, list(v, w))$dpsi[[k]]),
                   tolerance = 1e-6, info = lab)
    }
    g4 <- md$core(u0, gw, list(v, w))
    for (k in seq_len(K)) {
      expect_equal(g5$dpsi[[k]][[1L]], g4$dpsi[[k]], tolerance = 1e-11,
                   info = lab)
      expect_equal(g5$jacobian[[k]], g4$jacobian[[k]], tolerance = 1e-12,
                   info = lab)
    }
    # the negative control: a perturbed fifth-order array is detected
    bl <- gjf_blocks(df5$y, md$preds(u0), cfg[[1]])
    bl$U5 <- bl$U5 * 1.05
    bad <- md$fifth(u0, gw, list(v, w, x), blocks = bl)
    expect_false(isTRUE(all.equal(bad$curvature, ref, tolerance = 1e-6)),
                 info = lab)
  }
})
