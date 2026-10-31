# Shape constraints: the maps against cumulative sums written out here, the
# derivatives of b against one central difference per order, and the shape of
# the curves on a fine grid.

shp_xi <- c(0, 0.15, 0.4, 0.5, 0.8, 1.1, 1.5)

# beta from b by explicit loops, sharing no arithmetic with shape_map()
shp_beta <- function(shape, xi, b) {
  K <- length(xi); sg <- shape@shape_params$sign
  beta <- numeric(K); beta[1L] <- b[1L]
  if (S7::S7_inherits(shape, ShapeMonotone)) {
    for (j in 2:K) beta[j] <- beta[j - 1L] + (xi[j] - xi[j - 1L]) * sg * b[j]
  } else {
    s <- b[2L]
    beta[2L] <- beta[1L] + (xi[2L] - xi[1L]) * s
    for (j in 3:K) {
      s <- s + sg * b[j] * (xi[j] - xi[j - 2L]) / 2
      beta[j] <- beta[j - 1L] + (xi[j] - xi[j - 1L]) * s
    }
  }
  beta
}

shp_shapes <- function() {
  list(shape_increasing(), shape_decreasing(min_slope = 0.3), shape_convex(),
       shape_concave(link = linkfunctions7::softplus_link()))
}

test_that("shape constructors record name, sign, link and minimum slope", {
  expect_equal(shape_increasing()@shape_params$sign, 1)
  expect_equal(shape_decreasing()@shape_params$sign, -1)
  expect_equal(shape_decreasing(min_slope = 0.25)@shape_params$min_slope, 0.25)
  expect_equal(shape_convex()@shape_name, "convex")
  expect_equal(shape_concave()@shape_name, "concave")
  expect_true(S7::S7_inherits(shape_increasing(), ShapeMonotone))
  expect_true(S7::S7_inherits(shape_convex(), ShapeCurvature))
  expect_true(S7::S7_inherits(shape_convex(), shape_constraint))
  expect_output(print(shape_decreasing(min_slope = 0.1)), "min_slope = 0.1")
})

test_that("shape constructors reject a bad link or minimum slope", {
  expect_error(shape_increasing(link = linkfunctions7::identity_link()),
               "onto \\(0, Inf\\)")
  expect_error(shape_convex(link = linkfunctions7::logit_link()), "onto \\(0, Inf\\)")
  expect_error(shape_concave(link = "log"), "linkfunctions7 link")
  expect_error(shape_decreasing(min_slope = -0.1), "nonnegative")
  expect_error(shape_increasing(min_slope = c(0, 1)), "single")
  expect_error(shape_increasing(min_slope = NA_real_), "single")
  expect_error(shape_increasing(min_slope = Inf), "single")
})

test_that("shape_map agrees with cumulative sums written out directly", {
  set.seed(11)
  K <- length(shp_xi)
  for (sh in shp_shapes()) {
    m <- shape_map(sh, shp_xi)
    expect_equal(dim(m$A), c(K, K))
    expect_equal(unname(m$A[, 1L]), rep(1, K))
    for (i in 1:5) {
      b <- c(rnorm(1), exp(rnorm(K - 1L)))
      expect_equal(drop(m$A %*% b), shp_beta(sh, shp_xi, b), tolerance = 1e-12)
    }
    if (sh@shape_name %in% c("increasing", "decreasing")) {
      expect_equal(m$cons, 2:K)
      expect_equal(m$lower, c(0, rep(sh@shape_params$min_slope, K - 1L)))
    } else {
      expect_equal(m$cons, 3:K)
      expect_equal(m$lower, rep(0, K))
    }
  }
})

test_that("a constant coordinate is the second divided difference of the coefficients", {
  xi <- c(0, 0.3, 0.5, 0.9, 1.0, 1.6, 2)
  K <- length(xi)
  for (sg in c(1, -1)) {
    sh <- if (sg == 1) shape_convex() else shape_concave()
    beta <- drop(shape_map(sh, xi)$A %*% c(0.5, 1.2, rep(3, K - 2L)))
    s <- diff(beta) / diff(xi)
    mid <- (xi[-1L] + xi[-K]) / 2
    dd <- diff(s) / diff(mid)
    expect_equal(dd, rep(sg * 3, K - 2L), tolerance = 1e-10)
    expect_equal(s[1L], 1.2, tolerance = 1e-12)
  }
})

test_that("shape_map rejects a bad xi", {
  expect_error(shape_map(shape_increasing(), c(0, 1)), "length 3")
  expect_error(shape_map(shape_increasing(), c(0, 1, 1)), "strictly increasing")
  expect_error(shape_map(shape_convex(), c(0, NA, 1)), "strictly increasing")
})

test_that("shape_b derivatives match one central difference per order", {
  set.seed(12)
  K <- 8
  for (link in list(linkfunctions7::log_link(), linkfunctions7::softplus_link())) {
    cons <- 3:K; lower <- c(0, 0, rep(0.4, K - 2L))
    g <- rnorm(K)
    for (k in 1:5) {
      h <- 1e-4
      num <- (.shape_b(link, g + h, cons, lower, k - 1L) -
              .shape_b(link, g - h, cons, lower, k - 1L)) / (2 * h)
      expect_equal(num, .shape_b(link, g, cons, lower, k),
                   tolerance = 1e-6 * max(1, max(abs(num))))
    }
    # free coordinates: identity
    expect_equal(.shape_b(link, g, cons, lower, 0L)[1:2], g[1:2])
    expect_equal(.shape_b(link, g, cons, lower, 1L)[1:2], c(1, 1))
    expect_equal(.shape_b(link, g, cons, lower, 3L)[1:2], c(0, 0))
  }
  expect_error(.shape_b(linkfunctions7::log_link(), c(1, 2, 3), 2:3, rep(0, 3), 6L),
               "order must be 0 to 5")
})

test_that("shape_b of order zero is the lower bound plus the inverse link", {
  g <- c(-1, 0.3, 2)
  expect_equal(.shape_b(linkfunctions7::log_link(), g, 2:3, c(0, 0.5, 0.5), 0L),
               c(-1, 0.5 + exp(0.3), 0.5 + exp(2)), tolerance = 1e-12)
  expect_equal(.shape_b(linkfunctions7::softplus_link(), g, 1:3, rep(0, 3), 0L),
               log1p(exp(g)), tolerance = 1e-12)
})

test_that("shape_gamma inverts shape_b of order zero", {
  set.seed(13)
  K <- 9
  for (link in list(linkfunctions7::log_link(), linkfunctions7::softplus_link())) {
    for (cons in list(2:K, 3:K, integer(0))) {
      lower <- c(0, rep(0.2, K - 1L)); g <- rnorm(K)
      b <- .shape_b(link, g, cons, lower, 0L)
      expect_equal(.shape_gamma(link, b, cons, lower), g, tolerance = 1e-9)
    }
  }
})

test_that("random coefficients give a curve of the shape on a fine grid", {
  set.seed(14)
  x <- runif(200)
  for (K in c(8L, 12L)) {
    b <- basis7::smoother_basis(basis7::bspline_smooth(k = K), x)
    xi <- .scs_greville(b)
    xg <- seq(b@lower, b@upper, length.out = 4001)
    Bx <- basis7::basis_eval(b, xg); dx <- diff(xg[1:2])
    for (sh in shp_shapes()) {
      m <- shape_map(sh, xi)
      ms <- sh@shape_params$min_slope
      for (i in 1:30) {
        g <- rnorm(K, 0, 2)
        f <- drop(Bx %*% m$A %*% .shape_b(sh@link, g, m$cons, m$lower, 0L))
        sc <- max(abs(f), 1)
        switch(sh@shape_name,
          increasing = expect_true(min(diff(f)) / dx >= -1e-9 * sc / dx),
          decreasing = expect_true(max(diff(f)) / dx <= -ms + 1e-9 * sc / dx),
          convex = expect_true(min(diff(diff(f))) >= -1e-9 * sc),
          concave = expect_true(max(diff(diff(f))) <= 1e-9 * sc))
      }
    }
  }
})

test_that("shape_start returns a start of the right length and sign", {
  set.seed(15)
  x <- runif(100); xi <- seq(0, 1, length.out = 7); K <- 7
  yu <- 2 * x + rnorm(100, sd = 0.1)
  yd <- -2 * x + rnorm(100, sd = 0.1)
  yc <- 3 * (x - 0.5)^2 + rnorm(100, sd = 0.05)
  for (sh in shp_shapes()) {
    m <- shape_map(sh, xi)
    for (tg in list(NULL, yu, yd, yc)) {
      s <- shape_start(sh, x, tg, xi)
      expect_length(s, K)
      expect_true(all(is.finite(s)))
      expect_true(all(s[m$cons] > m$lower[m$cons]))
    }
  }
  # an increasing target gives equal positive slopes
  s <- shape_start(shape_increasing(), x, yu, xi)
  expect_true(all(s[-1] > 0))
  expect_equal(s[-1], rep(s[2], K - 1L))
  # a decreasing target under an increasing shape is floored, not negative
  s <- shape_start(shape_increasing(), x, yd, xi)
  expect_true(all(s[-1] > 0))
  # convex start on a convex target reproduces the curvature
  s <- shape_start(shape_convex(), x, yc, xi)
  expect_equal(s[3L], 6, tolerance = 0.3)
  # no target: flat at the floor
  expect_equal(shape_start(shape_increasing(), x, NULL, xi)[1], 0)
  expect_equal(shape_start(shape_convex(), x, NULL, xi)[1:2], c(0, 0))
})
