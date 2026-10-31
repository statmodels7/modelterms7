# scs(): the contribution against explicit sums, the block and its derivatives
# against central differences of the lower order, the shape beyond the data.

set.seed(21)
scs_n <- 100
scs_d <- data.frame(x = runif(scs_n))
scs_sm <- basis7::bspline_smooth(k = 8)
scs_shapes <- function() {
  list(shape_increasing(), shape_decreasing(min_slope = 0.2), shape_convex(),
       shape_concave(link = linkfunctions7::softplus_link()))
}
scs_build <- function(sh, d = scs_d) term_build(scs(x, sh, scs_sm), d)

# the contribution by explicit loops: beta from b, then the centered curve
scs_ref_value <- function(tm, g, d) {
  bp <- tm@build_info; sh <- tm@shape
  K <- length(bp$xi); xi <- bp$xi; sg <- sh@shape_params$sign
  b <- numeric(K)
  b[-1L] <- .shape_b(bp$link, g, bp$cons, bp$lower, 0L)
  beta <- numeric(K); beta[1L] <- b[1L]
  if (S7::S7_inherits(sh, ShapeMonotone)) {
    for (j in 2:K) beta[j] <- beta[j - 1L] + (xi[j] - xi[j - 1L]) * sg * b[j]
  } else {
    s <- b[2L]; beta[2L] <- beta[1L] + (xi[2L] - xi[1L]) * s
    for (j in 3:K) {
      s <- s + sg * b[j] * (xi[j] - xi[j - 2L]) / 2
      beta[j] <- beta[j - 1L] + (xi[j] - xi[j - 1L]) * s
    }
  }
  f <- drop(basis7::basis_eval(bp$basis, d$x) %*% beta)
  fb <- drop(basis7::basis_eval(bp$basis, scs_d$x) %*% beta)
  f - mean(fb)
}

test_that("coefficient names and count are K - 1", {
  tm <- scs_build(shape_increasing())
  expect_equal(term_coef_names(tm), sprintf("scs(x).g%d", 1:7))
  expect_equal(ncol(term_matrix(tm)), 7L)
  expect_equal(colnames(term_matrix(tm)), term_coef_names(tm))
  expect_equal(length(term_coef_start(tm)), 7L)
  tm12 <- term_build(scs(x, shape_convex(), basis7::bspline_smooth(k = 12)), scs_d)
  expect_equal(length(term_coef_names(tm12)), 11L)
  expect_output(print(tm), "built: 7 coefficients; shape increasing")
  expect_output(print(scs(x)), "specification")
})

test_that("the label can be set and defaults to the call", {
  expect_equal(scs(x)@label, "scs(x)")
  expect_equal(scs(x, label = "mine")@label, "mine")
  expect_error(scs(x, label = ""), "non-empty")
})

test_that("term_value equals the explicit sums and averages to zero", {
  set.seed(22)
  for (sh in scs_shapes()) {
    tm <- scs_build(sh)
    p <- length(term_coef_names(tm))
    for (i in 1:3) {
      g <- rnorm(p, 0, 0.7)
      r <- term_refresh(tm, g)
      ref <- scs_ref_value(r, g, scs_d)
      expect_equal(term_value(r), ref, tolerance = 1e-10)
      expect_equal(term_value(tm, coef = g), ref, tolerance = 1e-10)
      expect_equal(mean(term_value(r)), 0, tolerance = 1e-10)
      expect_equal(term_value(r, newdata = scs_d), ref, tolerance = 1e-10)
    }
  }
})

test_that("the block is the Jacobian of term_value in the coefficients", {
  set.seed(23)
  h <- 1e-6
  for (sh in scs_shapes()) {
    tm <- scs_build(sh); p <- length(term_coef_names(tm))
    g <- rnorm(p, 0, 0.7); r <- term_refresh(tm, g)
    J <- sapply(seq_len(p), function(j) {
      e <- replace(numeric(p), j, h)
      (term_value(tm, coef = g + e) - term_value(tm, coef = g - e)) / (2 * h)
    })
    expect_true(max(abs(term_matrix(r) - J)) / max(abs(J)) < 1e-6)
  }
})

test_that("term_predict at the build rows equals term_matrix", {
  set.seed(24)
  for (sh in scs_shapes()) {
    tm <- scs_build(sh); p <- length(term_coef_names(tm))
    r <- term_refresh(tm, rnorm(p, 0, 0.7))
    P <- term_predict(r, scs_d)
    expect_equal(P, term_matrix(r), tolerance = 1e-12)
    expect_equal(colnames(P), term_coef_names(r))
    expect_equal(unname(term_predict(r, scs_d[5:9, , drop = FALSE])),
                 unname(term_matrix(r)[5:9, ]), tolerance = 1e-12)
  }
})

test_that("term_block_deriv matches a difference of the refreshed block", {
  set.seed(25)
  h <- 1e-6
  for (sh in scs_shapes()) {
    tm <- scs_build(sh); p <- length(term_coef_names(tm))
    g <- rnorm(p, 0, 0.7); r <- term_refresh(tm, g); v <- rnorm(p)
    dX <- (term_matrix(term_refresh(tm, g + h * v)) -
           term_matrix(term_refresh(tm, g - h * v))) / (2 * h)
    D <- term_block_deriv(r, v = v)
    expect_equal(dim(D), c(scs_n, p))
    expect_true(max(abs(dX - D)) / max(abs(dX)) < 1e-6)
    expect_equal(term_block_deriv(r, coef = g, v = v), D)
  }
})

test_that("term_block_contract matches the contraction of the block difference", {
  set.seed(26)
  h <- 1e-6
  for (sh in scs_shapes()) {
    tm <- scs_build(sh); p <- length(term_coef_names(tm))
    g <- rnorm(p, 0, 0.7); r <- term_refresh(tm, g)
    A <- matrix(rnorm(scs_n * p), scs_n, p)
    cn <- sapply(seq_len(p), function(j) {
      e <- replace(numeric(p), j, h)
      sum(A * (term_matrix(term_refresh(tm, g + e)) -
               term_matrix(term_refresh(tm, g - e)))) / (2 * h)
    })
    cc <- term_block_contract(r, A = A)
    expect_length(cc, p)
    expect_true(max(abs(cn - cc)) / max(abs(cn)) < 1e-6)
  }
})

test_that("term_block_deriv2 matches a difference of term_block_deriv", {
  set.seed(27)
  h <- 1e-6
  for (sh in scs_shapes()) {
    tm <- scs_build(sh); p <- length(term_coef_names(tm))
    g <- rnorm(p, 0, 0.7); r <- term_refresh(tm, g)
    v <- rnorm(p); u <- rnorm(p)
    d2 <- (term_block_deriv(tm, coef = g + h * u, v = v) -
           term_block_deriv(tm, coef = g - h * u, v = v)) / (2 * h)
    D2 <- term_block_deriv2(r, v = v, u = u)
    expect_equal(dim(D2), c(scs_n, p))
    expect_true(max(abs(d2 - D2)) / max(abs(d2)) < 1e-6)
    expect_equal(term_block_deriv2(r, v = u, u = v), D2, tolerance = 1e-12)
  }
})

test_that("the curve has the shape on a fine grid, beyond the data", {
  set.seed(28)
  xg <- seq(-0.4, 1.4, length.out = 4001); dx <- diff(xg[1:2])
  dg <- data.frame(x = xg)
  for (sh in scs_shapes()) {
    tm <- scs_build(sh); p <- length(term_coef_names(tm))
    ms <- sh@shape_params$min_slope
    for (i in 1:25) {
      g <- rnorm(p, 0, 2)
      f <- term_value(tm, coef = g, newdata = dg)
      sc <- max(abs(f), 1)
      switch(sh@shape_name,
        increasing = expect_true(min(diff(f)) >= -1e-9 * sc),
        decreasing = expect_true(max(diff(f)) / dx <= -ms + 1e-9 * sc / dx),
        convex = expect_true(min(diff(diff(f))) >= -1e-9 * sc),
        concave = expect_true(max(diff(diff(f))) <= 1e-9 * sc))
    }
  }
})

test_that("beyond the ends the curve continues linearly", {
  set.seed(29)
  tm <- scs_build(shape_convex()); p <- length(term_coef_names(tm))
  g <- rnorm(p)
  lo <- tm@build_info$basis@lower; up <- tm@build_info$basis@upper
  xr <- up + c(0.1, 0.3, 0.5, 0.7)
  f <- term_value(tm, coef = g, newdata = data.frame(x = xr))
  expect_equal(diff(f, differences = 2), rep(0, 2), tolerance = 1e-10)
  xl <- lo - c(0.7, 0.5, 0.3, 0.1)
  f <- term_value(tm, coef = g, newdata = data.frame(x = xl))
  expect_equal(diff(f, differences = 2), rep(0, 2), tolerance = 1e-10)
})

test_that("a minimum slope bounds the decrease from above", {
  set.seed(30)
  tm <- scs_build(shape_decreasing(min_slope = 0.5)); p <- length(term_coef_names(tm))
  xg <- seq(0, 1, length.out = 3001)
  for (i in 1:20) {
    f <- term_value(tm, coef = rnorm(p, 0, 2), newdata = data.frame(x = xg))
    expect_true(max(diff(f) / diff(xg)) <= -0.5 + 1e-6)
  }
})

test_that("term_penalty is the quadratic penalty of differences of the constrained coordinates", {
  tm <- scs_build(shape_increasing())
  pen <- term_penalty(tm)
  expect_true(S7::S7_inherits(pen, penalties7::QuadraticPenalty))
  P <- as.matrix(pen@P)
  D <- matrix(0, 6, 7)
  for (i in 1:6) { D[i, i] <- -1; D[i, i + 1L] <- 1 }
  expect_equal(unname(P), crossprod(D), tolerance = 1e-12)
  # convex: the free slope carries no penalty
  Pc <- as.matrix(term_penalty(scs_build(shape_convex()))@P)
  expect_equal(unname(Pc[1L, ]), rep(0, 7))
  expect_equal(unname(Pc[, 1L]), rep(0, 7))
  Dc <- matrix(0, 5, 6)
  for (i in 1:5) { Dc[i, i] <- -1; Dc[i, i + 1L] <- 1 }
  expect_equal(unname(Pc[-1, -1]), crossprod(Dc), tolerance = 1e-12)
})

test_that("term_coef_start is finite, of the right length and increasing for an increasing target", {
  set.seed(31)
  tm <- scs_build(shape_increasing())
  s0 <- term_coef_start(tm)
  expect_length(s0, 7L); expect_true(all(is.finite(s0)))
  y <- 3 * scs_d$x + rnorm(scs_n, sd = 0.1)
  s1 <- term_coef_start(tm, target = y)
  expect_length(s1, 7L); expect_true(all(is.finite(s1)))
  xg <- seq(0, 1, length.out = 501)
  f <- term_value(tm, coef = s1, newdata = data.frame(x = xg))
  expect_true(all(diff(f) > 0))
  expect_equal(unname(stats::coef(stats::lm(f ~ xg))[2]), 3, tolerance = 0.2)
  for (sh in scs_shapes()) {
    t2 <- scs_build(sh); p2 <- length(term_coef_names(t2))
    for (tg in list(NULL, y, -y, (scs_d$x - 0.5)^2)) {
      s <- term_coef_start(t2, target = tg)
      expect_length(s, p2); expect_true(all(is.finite(s)))
    }
  }
})

test_that("scs() rejects what it cannot fit", {
  expect_error(scs(x, shape = "increasing"), "shape constraint")
  expect_error(scs(x, smoother = basis7::fourier_smooth()), "bspline_smooth")
  expect_error(scs(x, by = g), "'by' is not yet supported")
  expect_error(scs(x, sparse = TRUE), "dense")
  expect_error(shape_increasing(link = linkfunctions7::identity_link()), "onto")
  expect_error(shape_decreasing(min_slope = -1), "nonnegative")
})

test_that("a covariate with a missing value is rejected at the build", {
  d <- scs_d; d$x[3] <- NA
  expect_error(term_build(scs(x), d), "missing or infinite")
  d <- scs_d; d$x[3] <- Inf
  expect_error(term_build(scs(x), d), "missing or infinite")
})

test_that("arguments of the wrong length are rejected", {
  tm <- scs_build(shape_increasing()); p <- 7L
  expect_error(term_refresh(tm, rnorm(p + 1L)), "must have length 7")
  expect_error(term_value(tm, coef = rnorm(p - 1L)), "must have length 7")
  expect_error(term_block_deriv(tm, v = rnorm(p + 1L)), "'v' must have length 7")
  expect_error(term_block_contract(tm, A = matrix(0, scs_n, p + 1L)), "'A' must be")
  expect_error(term_block_contract(tm, A = matrix(0, scs_n - 1L, p)), "'A' must be")
  expect_error(term_block_deriv2(tm, v = rnorm(p), u = rnorm(p - 1L)),
               "must have length 7")
})
