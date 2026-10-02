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
