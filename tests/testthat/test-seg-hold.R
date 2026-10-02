test_that("seg_polish_exact() reaches the minimum of the profile over every interval", {
  design <- function(kind, x, q) {
    if (kind == "jseg") cbind(1, x, as.numeric(x > q), pmax(x - q, 0))
    else cbind(1, as.numeric(x > q))
  }
  rssq <- function(kind, x, y, q) sum(stats::lm.fit(design(kind, x, q), y)$residuals^2)
  for (kind in c("jseg", "jump")) for (s in 1:4) {
    set.seed(s)
    n <- 150
    x <- runif(n, 0, 10)
    y <- if (kind == "jseg") 1 + 0.3 * x + ifelse(x > 6, 1.5 - 0.6 * (x - 6), 0) +
      rnorm(n, sd = 0.5) else 1 + 1.5 * (x > 6) + rnorm(n, sd = 0.5)
    d <- data.frame(x = x)
    b <- term_build(if (kind == "jseg") jseg(x, psi = 3) else jump(x, psi = 3), d)
    e <- seg_polish_exact(b, y)
    u <- sort(unique(x))
    m <- (u[-1] + u[-length(u)]) / 2
    lim <- b@blueprint$lim
    m <- m[m > lim[1] & m < lim[2]]
    r <- vapply(m, function(q) rssq(kind, x, y, q), 1)
    expect_equal(rssq(kind, x, y, seg_psi(e)), min(r), tolerance = 1e-10)
  }
})

test_that("seg_polish_exact() agrees with a two-dimensional search under weights", {
  set.seed(2)
  n <- 120
  x <- runif(n, 0, 10)
  y <- 2 * (x > 3) - 1.5 * (x > 7) + rnorm(n, sd = 0.4)
  w <- rexp(n)
  b <- term_build(jump(x, npsi = 2, psi = c(2, 5)), data.frame(x = x))
  e <- seg_polish_exact(b, y, weights = w)
  f <- function(q) {
    Z <- cbind(1, x > q[1], x > q[2])
    sum(w * stats::lm.wfit(Z, y, w)$residuals^2)
  }
  u <- sort(unique(x))
  m <- (u[-1] + u[-length(u)]) / 2
  lim <- b@blueprint$lim
  m <- m[m > lim[1] & m < lim[2]]
  best <- Inf
  for (i in seq_along(m)) for (j in seq_along(m)) if (j > i) {
    best <- min(best, f(c(m[i], m[j])))
  }
  expect_equal(f(seg_psi(e)), best, tolerance = 1e-10)
})

test_that("a held term carries the exact design and holds its position's slot", {
  set.seed(1)
  n <- 100
  x <- sort(runif(n, 0, 10))
  d <- data.frame(x = x)
  b <- term_build(jseg(x, psi = 6), d)
  expect_identical(term_held(b), integer(0))
  expect_false(term_jacobian_block(b))
  h <- seg_hold(b)
  expect_true(term_jacobian_block(h))
  expect_identical(term_held(h), 4L)
  X <- term_matrix(h)
  p <- seg_psi(h)
  expect_equal(unname(X[, 1]), x)
  expect_equal(unname(X[, 2]), pmax(x - p, 0))
  expect_equal(unname(X[, 3]), as.numeric(x > p))
  expect_true(all(X[, 4] == 0))
  # the position does not move with the coefficients, and the contribution is
  # the block times them
  cf <- c(0.2, -0.5, 1.3, 99)
  r <- term_refresh(h, cf)
  expect_equal(seg_psi(r), p)
  expect_equal(term_value(r), as.numeric(term_matrix(r) %*% cf))
  expect_true(term_converged(r))
  expect_false(term_stalled(r))
  # released, the slot carries -delta psi and the read-off returns the position
  rl <- seg_hold(r, hold = FALSE)
  expect_identical(term_held(rl), integer(0))
  expect_equal(rl@blueprint$coef[4], -1.3 * p)
})

test_that("seg_hold() rejects a continuous or smoothed term", {
  d <- data.frame(x = seq(0, 10, length.out = 50))
  expect_error(seg_hold(term_build(seg(x), d)), "sharp jump")
  expect_error(seg_hold(term_build(jump(x, smoothed = numericals7::smooth_probit()), d)),
               "sharp jump")
})

test_that("seg_profile_intervals() is the brute-force profile over every interval", {
  set.seed(4)
  n <- 120
  x <- runif(n, 0, 10)
  y <- 1 + 0.3 * x + 2 * (x > 3) - 1.5 * (x > 7) + rnorm(n, sd = 0.4)
  w <- rexp(n)
  d <- data.frame(x = x)
  u <- sort(unique(x))
  m <- (u[-1] + u[-length(u)]) / 2
  # one break-point of a jseg
  b <- term_build(jseg(x, psi = 4), d)
  pr <- seg_profile_intervals(b, y, weights = w)
  lim <- b@blueprint$lim
  expect_equal(pr$psi, m[m > lim[1] & m < lim[2]])
  ref <- vapply(pr$psi, function(q) {
    Z <- cbind(1, x, pmax(x - q, 0), as.numeric(x > q))
    sum(w * stats::lm.wfit(Z, y, w)$residuals^2)
  }, 1)
  expect_equal(pr$rss, ref, tolerance = 1e-10)
  # the second break-point of a jump, the first held where it is
  b2 <- term_build(jump(x, npsi = 2, psi = c(3.1, 6)), d)
  pr2 <- seg_profile_intervals(b2, y, k = 2)
  q1 <- seg_psi(b2)[1]
  ref2 <- vapply(pr2$psi, function(q) {
    sum(stats::lm.fit(cbind(1, x > q1, x > q), y)$residuals^2)
  }, 1)
  ok <- is.finite(pr2$rss)
  expect_equal(pr2$rss[ok], ref2[ok], tolerance = 1e-10)
  expect_error(seg_profile_intervals(b2, y, k = 3), "'k' must be one of 1 to 2")
})

test_that("a held developed jump is polished to the exhaustive per-group minimum", {
  set.seed(1)
  g <- factor(rep(c("a", "b", "c"), each = 30))
  x <- runif(90, 0, 10)
  gi <- as.integer(g)
  y <- 1 + 1.5 * (x > c(3, 5, 7)[gi]) + rnorm(90, sd = 0.5)
  d <- data.frame(x = x, g = g)
  b <- seg_hold(term_build(jump(x, psi ~ 0 + g), d))
  e <- seg_polish_exact(b, y)
  rss <- function(q) sum(stats::.lm.fit(cbind(1, x > q[gi]), y)$residuals^2)
  cand <- lapply(1:3, function(j) {
    u <- sort(unique(x[gi == j]))
    m <- (u[-1] + u[-length(u)]) / 2
    lim <- stats::quantile(x[gi == j], c(0.05, 0.95))
    m[m > lim[1] & m < lim[2]]
  })
  best <- Inf
  for (a in cand[[1]]) for (bb in cand[[2]]) {
    r <- vapply(cand[[3]], function(v) rss(c(a, bb, v)), 1)
    best <- min(best, r)
  }
  q <- vapply(1:3, function(j) e@blueprint$psi[which(gi == j)[1], 1], 1)
  expect_equal(rss(q), best, tolerance = 1e-10)
  # the held positions are reported one per level
  rd <- term_readable(e, e@blueprint$coef)
  expect_identical(rd$name[grepl("^psi", rd$name)], c("psi1.ga", "psi1.gb", "psi1.gc"))
  expect_equal(rd$value[grepl("^psi", rd$name)], unname(q), tolerance = 1e-12)
})

test_that("a held jump developed on a continuous covariate reaches the grid minimum", {
  set.seed(2)
  n <- 150
  x <- runif(n, 0, 10)
  z <- runif(n, -1, 1)
  y <- 1 + 1.5 * (x > 5 + 1.5 * z) + rnorm(n, sd = 0.4)
  b <- seg_hold(term_build(jump(x, psi ~ z), data.frame(x = x, z = z)))
  e <- seg_polish_exact(b, y)
  rss <- function(p) sum(stats::.lm.fit(cbind(1, x > p[1] + p[2] * z), y)$residuals^2)
  gr <- expand.grid(p0 = seq(3, 7, by = 0.05), p1 = seq(-1, 4, by = 0.05))
  grid_min <- min(apply(gr, 1, rss))
  expect_lte(rss(e@blueprint$pk[[1]]), grid_min + 1e-10)
})

test_that("a developed term is polished only when held", {
  d <- data.frame(x = runif(60, 0, 10), g = factor(rep(1:2, 30)))
  b <- term_build(jump(x, psi ~ 0 + g), d)
  expect_error(seg_polish_exact(b, rnorm(60)), "polished held")
})

test_that("seg_polish_exact() reaches the profile minimum of a change of slope", {
  set.seed(5)
  n <- 150
  x <- runif(n, 0, 10)
  y <- 1 + 0.5 * x - 0.9 * pmax(x - 4, 0) + 0.7 * pmax(x - 7.5, 0) +
    rnorm(n, sd = 0.3)
  d <- data.frame(x = x)
  rss <- function(q) {
    sum(stats::lm.fit(cbind(1, x, sapply(q, function(v) pmax(x - v, 0))),
                      y)$residuals^2)
  }
  # one break-point: the profile over every interval against a fine grid
  b <- term_build(seg(x, psi = 2), d)
  pr <- seg_profile_intervals(b, y)
  expect_equal(pr$rss, vapply(pr$psi, rss, 1), tolerance = 1e-9)
  g <- seq(min(pr$psi), max(pr$psi), length.out = 5000)
  expect_lte(min(pr$rss), min(vapply(g, rss, 1)) + 1e-10)
  # two: from a poor start the polish reaches the two-dimensional minimum
  e <- seg_polish_exact(term_build(seg(x, npsi = 2, psi = c(2, 3)), d), y)
  g2 <- expand.grid(a = seq(3, 5, by = 0.02), b = seq(6.5, 8.5, by = 0.02))
  expect_lte(rss(seg_psi(e)), min(apply(g2, 1, rss)) + 1e-10)
  # weighted
  w <- rexp(n)
  pw <- seg_profile_intervals(b, y, weights = w)
  ref <- vapply(pw$psi, function(q) {
    Z <- cbind(1, x, pmax(x - q, 0))
    sum(w * stats::lm.wfit(Z, y, w)$residuals^2)
  }, 1)
  expect_equal(pw$rss, ref, tolerance = 1e-9)
})

test_that("a break-point term starts from least squares when given a target", {
  set.seed(6)
  x <- 1900 + runif(120, 0, 100)
  y <- 0.01 * (x - 1900) - 0.03 * pmax(x - 1960, 0) + rnorm(120, sd = 0.1)
  b <- term_build(seg(x, psi = 1960), data.frame(x = x))
  expect_identical(unname(term_coef_start(b)[2]), 1)
  cf <- term_coef_start(b, target = y)
  l <- stats::lm.fit(cbind(1, x, pmax(x - 1960, 0)), y)$coefficients
  expect_equal(unname(cf[1:2]), unname(l[2:3]), tolerance = 1e-8)
  expect_equal(unname(cf[3]), 1960)
  j <- term_build(jump(x, psi = 1960), data.frame(x = x))
  cj <- term_coef_start(j, target = y)
  lj <- stats::lm.fit(cbind(1, as.numeric(x > 1960)), y)$coefficients
  expect_equal(unname(cj[1]), unname(lj[2]), tolerance = 1e-8)
  # the slot carries -delta psi, which reads back the starting position
  expect_equal(unname(-cj[2] / cj[1]), 1960)
})
