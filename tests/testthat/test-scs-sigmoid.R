# shape_sigmoid() and the valley of shape_unimodal(): the sigmoid map against
# slopes summed from both ends, the shape on a fine grid, the block
# derivatives against central differences, the location from the pilot.

set.seed(81)
sg_d <- data.frame(x = runif(120))
sg_sm <- basis7::bspline_smooth(k = 10)

test_that("the constructors check their arguments", {
  expect_identical(shape_sigmoid()@shape_params$sign, 1)
  expect_identical(shape_sigmoid(decreasing = TRUE)@shape_params$sign, -1)
  expect_identical(shape_unimodal(valley = TRUE)@shape_params$sign, -1)
  expect_error(shape_sigmoid(inflection = "a"), "single finite")
  expect_error(shape_sigmoid(decreasing = NA), "TRUE or FALSE")
  expect_error(shape_unimodal(valley = 1), "TRUE or FALSE")
  expect_true(S7::S7_inherits(shape_sigmoid(), ShapeLocated))
})

test_that("the sigmoid map sums the slopes from both ends toward the inflection", {
  xi <- c(0, 0.1, 0.25, 0.4, 0.55, 0.7, 0.85, 1)
  K <- length(xi)
  w <- c(0, 0, (xi[3:K] - xi[1:(K - 2L)]) / 2)
  set.seed(82)
  for (dec in c(FALSE, TRUE)) for (M in c(2L, 4L, K - 1L, K)) {
    sh <- shape_sigmoid(index = M, decreasing = dec)
    sg <- if (dec) -1 else 1
    m <- shape_map(sh, xi)
    b <- c(0.3, exp(rnorm(K - 1L)))
    s <- numeric(K)
    s[2L] <- sg * b[2L]
    if (M >= 3L) for (j in 3:M) s[j] <- s[j - 1L] + sg * b[j] * w[j]
    if (M < K) {
      s[K] <- sg * b[K]
      if (K - 1L >= M + 1L) for (j in (K - 1L):(M + 1L)) s[j] <- s[j + 1L] + sg * b[j] * w[j + 1L]
    }
    beta <- numeric(K); beta[1L] <- b[1L]
    for (j in 2:K) beta[j] <- beta[j - 1L] + (xi[j] - xi[j - 1L]) * s[j]
    expect_equal(drop(m$A %*% b), beta, tolerance = 1e-12)
    # every slope has the sign of the shape; |slope| grows up to M, shrinks after
    expect_true(all(sg * s[2:K] > 0))
    if (M >= 3L) expect_true(all(diff(abs(s[2:M])) > 0))
    if (M + 1L < K) expect_true(all(diff(abs(s[(M + 1L):K])) < 0))
  }
  expect_error(shape_map(shape_sigmoid(index = 1), xi), "outside")
})

test_that("random coefficients give a monotone curve whose slope has one peak", {
  set.seed(83)
  xg <- seq(-0.2, 1.2, length.out = 3001)
  for (dec in c(FALSE, TRUE)) for (infl in c(0.3, 0.7)) {
    tm <- term_build(scs(x, shape_sigmoid(inflection = infl, decreasing = dec), sg_sm), sg_d)
    sg <- if (dec) -1 else 1
    p <- length(term_coef_names(tm))
    for (i in 1:20) {
      f <- term_value(tm, coef = rnorm(p, 0, 1.5), newdata = data.frame(x = xg))
      sl <- sg * diff(f)
      expect_true(all(sl >= -1e-12 * max(abs(f))))
      d2 <- diff(sl)
      s2 <- sign(d2[abs(d2) > 1e-12 * max(abs(f))])
      expect_lte(sum(diff(s2) != 0), 1L)
    }
  }
})

test_that("the valley has one minimum", {
  set.seed(84)
  tm <- term_build(scs(x, shape_unimodal(mode = 0.6, valley = TRUE), sg_sm), sg_d)
  xg <- seq(0, 1, length.out = 2001)
  for (i in 1:20) {
    f <- term_value(tm, coef = rnorm(length(term_coef_names(tm)), 0, 1.5),
                    newdata = data.frame(x = xg))
    df <- diff(f); s <- sign(df[abs(df) > 1e-12 * max(abs(f))])
    expect_lte(sum(diff(s) != 0), 1L)
    if (length(s)) expect_true(s[1L] <= 0)
  }
})

test_that("the block and its derivatives match central differences", {
  set.seed(85)
  h <- 1e-6
  for (sh in list(shape_sigmoid(inflection = 0.4), shape_sigmoid(inflection = 0.6, decreasing = TRUE),
                  shape_unimodal(mode = 0.5, valley = TRUE))) {
    tm <- term_build(scs(x, sh, sg_sm), sg_d)
    p <- length(term_coef_names(tm))
    g <- rnorm(p, 0, 0.7); v <- rnorm(p); u <- rnorm(p)
    r <- term_refresh(tm, g)
    J <- sapply(seq_len(p), function(j) {
      e <- replace(numeric(p), j, h)
      (term_value(tm, coef = g + e) - term_value(tm, coef = g - e)) / (2 * h)
    })
    expect_lt(max(abs(J - term_matrix(r))) / max(abs(J)), 1e-6)
    dX <- (term_matrix(term_refresh(tm, g + h * v)) -
           term_matrix(term_refresh(tm, g - h * v))) / (2 * h)
    expect_lt(max(abs(dX - term_block_deriv(r, v = v))) / max(abs(dX)), 1e-6)
    d2 <- (term_block_deriv(tm, coef = g + h * u, v = v) -
           term_block_deriv(tm, coef = g - h * u, v = v)) / (2 * h)
    expect_lt(max(abs(d2 - term_block_deriv2(r, v = v, u = u))) / max(abs(d2)), 1e-6)
  }
})

test_that("the pilot places the inflection at the steepest point", {
  tm <- term_build(scs(x, shape_sigmoid(), sg_sm), sg_d)
  expect_true(S7::S7_inherits(term_pilot(tm), SmoothTerm))
  xi <- tm@build_info$xi
  pil <- stats::plogis((sg_d$x - 0.62) / 0.05)
  v <- term_variants(tm, pilot = pil)
  M0 <- (2:length(xi))[which.min(abs(xi[-1L] - 0.62))]
  expect_true(M0 %in% vapply(v, function(t) t@shape@shape_params$index, 1L))
  expect_true(all(vapply(v, function(t) t@shape@shape_params$index, 1L) >= 2L))
  # the valley's pilot is its minimum
  tv <- term_build(scs(x, shape_unimodal(valley = TRUE), sg_sm), sg_d)
  vv <- term_variants(tv, pilot = (sg_d$x - 0.3)^2)
  Mv <- which.min(abs(xi - sg_d$x[which.min((sg_d$x - 0.3)^2)]))
  expect_true(Mv %in% vapply(vv, function(t) t@shape@shape_params$index, 1L))
})
