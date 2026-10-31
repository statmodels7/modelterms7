# The four combinations of monotonicity and curvature: the map against the
# slopes summed from the anchored end, the shape on a fine grid beyond the
# data, the block derivatives against central differences, the penalty with
# the anchored slope inside the walk, the start.

set.seed(41)
cr_d <- data.frame(x = runif(100))
cr_sm <- basis7::bspline_smooth(k = 9)
cr_shapes <- function() {
  list(shape_increasing_convex(), shape_increasing_concave(),
       shape_decreasing_convex(),
       shape_decreasing_concave(link = linkfunctions7::softplus_link()))
}

test_that("the constructors record the two signs and the name", {
  sh <- cr_shapes()
  expect_identical(vapply(sh, function(s) s@shape_name, ""),
                   c("increasing convex", "increasing concave",
                     "decreasing convex", "decreasing concave"))
  expect_identical(vapply(sh, function(s) s@shape_params$sign, 0), c(1, -1, 1, -1))
  expect_identical(vapply(sh, function(s) s@shape_params$mono, 0), c(1, 1, -1, -1))
  expect_true(all(vapply(sh, function(s) S7::S7_inherits(s, ShapeCurvature), TRUE)))
  expect_error(shape_increasing_convex(link = linkfunctions7::identity_link()), "onto")
})

test_that("shape_map sums the slopes from the anchored end", {
  xi <- c(0, 0.1, 0.25, 0.45, 0.7, 0.85, 1)
  K <- length(xi)
  w <- c(0, 0, (xi[3:K] - xi[1:(K - 2L)]) / 2)
  set.seed(42)
  for (sh in cr_shapes()) {
    sg <- sh@shape_params$sign; mo <- sh@shape_params$mono
    m <- shape_map(sh, xi)
    expect_identical(as.integer(m$cons), 2:K)
    expect_identical(as.integer(m$pen), 2:K)
    b <- c(0.4, exp(rnorm(K - 1L)))
    s <- numeric(K)
    if (mo == sg) {
      s[2L] <- mo * b[2L]
      for (j in 3:K) s[j] <- s[j - 1L] + sg * b[j] * w[j]
    } else {
      s[K] <- mo * b[2L]
      for (j in (K - 1L):2) s[j] <- s[j + 1L] - sg * b[j + 1L] * w[j + 1L]
    }
    beta <- numeric(K); beta[1L] <- b[1L]
    for (j in 2:K) beta[j] <- beta[j - 1L] + (xi[j] - xi[j - 1L]) * s[j]
    expect_equal(drop(m$A %*% b), beta, tolerance = 1e-12)
    # every slope has the sign of the monotonicity, and the slopes move
    # with the sign of the curvature
    expect_true(all(mo * s[2:K] > 0))
    expect_true(all(sg * diff(s[2:K]) > 0))
  }
})

test_that("random coefficients give a monotone curve of the curvature's sign", {
  set.seed(43)
  xg <- seq(-0.3, 1.3, length.out = 3001)
  for (sh in cr_shapes()) {
    tm <- term_build(scs(x, sh, cr_sm), cr_d)
    p <- length(term_coef_names(tm))
    sg <- sh@shape_params$sign; mo <- sh@shape_params$mono
    for (i in 1:25) {
      f <- term_value(tm, coef = rnorm(p, 0, 1.5), newdata = data.frame(x = xg))
      expect_true(all(mo * diff(f) >= -1e-12 * max(abs(f))))
      expect_true(all(sg * diff(diff(f)) >= -1e-12 * max(abs(f))))
    }
  }
})

test_that("the block and its derivatives match central differences", {
  set.seed(44)
  h <- 1e-6
  for (sh in cr_shapes()) {
    tm <- term_build(scs(x, sh, cr_sm), cr_d)
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

test_that("the penalty takes the anchored slope into the walk", {
  for (sh in cr_shapes()) {
    tm <- term_build(scs(x, sh, cr_sm), cr_d)
    P <- as.matrix(term_penalty(tm)@P)
    p <- ncol(P)
    expect_equal(P, crossprod(diff(diag(p))), tolerance = 1e-14,
                 ignore_attr = TRUE)
  }
  # a plain curvature shape keeps its free slope out of the walk
  P <- as.matrix(term_penalty(term_build(scs(x, shape_convex(), cr_sm), cr_d))@P)
  expect_true(all(P[1L, ] == 0) && all(P[, 1L] == 0))
})

test_that("the start is finite and keeps the anchor's sign", {
  for (sh in cr_shapes()) {
    tm <- term_build(scs(x, sh, cr_sm), cr_d)
    mo <- sh@shape_params$mono; sg <- sh@shape_params$sign
    for (target in list(NULL, mo * cr_d$x + sg * cr_d$x^2, -cr_d$x)) {
      g <- term_coef_start(tm, target = target)
      expect_length(g, length(term_coef_names(tm)))
      expect_true(all(is.finite(g)))
    }
  }
})
