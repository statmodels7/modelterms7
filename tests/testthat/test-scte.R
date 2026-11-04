# scte(): the block and its derivatives against central differences, the
# shapes and the sign of the interaction on a grid, the reflection of a
# decreasing margin, the coordinates that leave the constant to the
# intercept, and the rejections.

set.seed(51)
st_d <- data.frame(x = runif(200), z = runif(200))
st_sm <- basis7::bspline_smooth(k = 5)
st_grid <- expand.grid(x = seq(0.02, 0.98, length.out = 25), z = seq(0.02, 0.98, length.out = 25))
st_I <- shape_increasing(); st_D <- shape_decreasing()

test_that("the block and its derivatives agree with central differences", {
  set.seed(52)
  for (sh in list(list(st_I, NULL), list(NULL, st_D), list(st_I, st_I), list(st_D, st_I),
                  list(shape_convex(), NULL), list(st_I, shape_concave()),
                  list(shape_convex(), shape_concave()))) {
    tm <- term_build(scte(x, z, shapes = sh, smooths = st_sm), st_d)
    p <- length(term_coef_names(tm)); g <- stats::rnorm(p, sd = 0.7); h <- 1e-5
    e <- function(j) replace(numeric(p), j, h)
    J <- sapply(seq_len(p), function(j) (term_value(tm, coef = g + e(j)) - term_value(tm, coef = g - e(j))) / (2 * h))
    v <- stats::rnorm(p); u <- stats::rnorm(p)
    M <- function(gg) term_matrix(term_refresh(tm, gg))
    rel <- function(a, b) max(abs(a - b)) / max(abs(b))
    r <- term_refresh(tm, g)
    expect_lt(rel(term_matrix(r), J), 1e-8)
    expect_lt(rel(term_block_deriv(tm, coef = g, v = v), (M(g + h * v) - M(g - h * v)) / (2 * h)), 1e-8)
    D2 <- (term_block_deriv(tm, coef = g + h * u, v = v) - term_block_deriv(tm, coef = g - h * u, v = v)) / (2 * h)
    expect_lt(rel(term_block_deriv2(tm, coef = g, v = v, u = u), D2), 1e-8)
    s <- stats::rnorm(p); q <- stats::rnorm(p)
    D3 <- (term_block_deriv2(tm, coef = g + h * s, v = v, u = u) -
           term_block_deriv2(tm, coef = g - h * s, v = v, u = u)) / (2 * h)
    expect_lt(rel(term_block_deriv3(tm, coef = g, v = v, u = u, s = s), D3), 1e-8)
    D4 <- (term_block_deriv3(tm, coef = g + h * q, v = v, u = u, s = s) -
           term_block_deriv3(tm, coef = g - h * q, v = v, u = u, s = s)) / (2 * h)
    expect_lt(rel(term_block_deriv4(tm, coef = g, v = v, u = u, s = s, r = q), D4), 1e-8)
    A <- matrix(stats::rnorm(200 * p), 200)
    ctr <- sapply(seq_len(p), function(k) sum(A * (M(g + e(k)) - M(g - e(k)))) / (2 * h))
    expect_lt(rel(term_block_contract(tm, coef = g, A = A), ctr), 1e-8)
    expect_equal(term_predict(r, st_d), term_matrix(r), tolerance = 1e-14, ignore_attr = TRUE)
    # the constant is left to the intercept: the block beside a column of ones has full rank
    expect_identical(qr(cbind(1, term_matrix(r)))$rank, p + 1L)
  }
})

test_that("every coefficient vector gives the shapes and the sign of the interaction", {
  set.seed(53)
  sgn <- function(sh) if (is.null(sh)) 0 else sh@shape_params$sign
  for (sh in list(list(st_I, NULL), list(NULL, st_I), list(st_D, NULL), list(st_I, st_I),
                  list(st_I, st_D), list(st_D, st_D))) {
    tm <- term_build(scte(x, z, shapes = sh, smooths = st_sm), st_d)
    p <- length(term_coef_names(tm))
    for (r in 1:3) {
      f <- matrix(term_value(tm, coef = stats::rnorm(p, sd = 2), newdata = st_grid), 25, 25)
      if (sgn(sh[[1L]]) != 0) expect_true(all(sgn(sh[[1L]]) * diff(f) >= -1e-10))
      if (sgn(sh[[2L]]) != 0) expect_true(all(sgn(sh[[2L]]) * t(diff(t(f))) >= -1e-10))
      if (sgn(sh[[1L]]) != 0 && sgn(sh[[2L]]) != 0) {
        expect_true(all(sgn(sh[[1L]]) * sgn(sh[[2L]]) * diff(t(diff(t(f)))) >= -1e-10))
      }
    }
  }
})

test_that("every coefficient vector gives the curvature shapes", {
  set.seed(55)
  cx <- shape_convex(); cv <- shape_concave()
  need <- function(sh) {
    if (is.null(sh)) return(c(0, 0))
    p <- sh@shape_params
    if (S7::S7_inherits(sh, ShapeMonotone)) c(p$sign, 0) else c(if (is.null(p$mono)) 0 else p$mono, p$sign)
  }
  for (sh in list(list(cx, NULL), list(NULL, cv), list(shape_increasing_concave(), NULL),
                  list(cx, cx), list(cv, cv), list(st_I, cv), list(st_D, cx),
                  list(shape_decreasing_convex(), shape_increasing_concave()), list(cx, cv))) {
    tm <- term_build(scte(x, z, shapes = sh, smooths = st_sm), st_d)
    p <- length(term_coef_names(tm))
    nx <- need(sh[[1L]]); nz <- need(sh[[2L]])
    for (r in 1:5) {
      f <- matrix(term_value(tm, coef = stats::rnorm(p, -1, 1.5), newdata = st_grid), 25, 25)
      tol <- -1e-10 * max(abs(f))
      if (nx[1L] != 0) expect_true(all(nx[1L] * diff(f) >= tol))
      if (nx[2L] != 0) expect_true(all(nx[2L] * diff(f, differences = 2) >= tol))
      if (nz[1L] != 0) expect_true(all(nz[1L] * t(diff(t(f))) >= tol))
      if (nz[2L] != 0) expect_true(all(nz[2L] * t(diff(t(f), differences = 2)) >= tol))
    }
    # the constant is left to the intercept
    r0 <- term_refresh(tm, stats::rnorm(p, -1, 1))
    expect_identical(qr(cbind(1, term_matrix(r0)))$rank, p + 1L)
  }
  # convex with concave: the 4 x 4 elements constrained in both directions are dropped
  expect_length(term_coef_names(term_build(scte(x, z, shapes = list(cx, cv), smooths = st_sm), st_d)),
                15L)
})

test_that("a minimum slope holds everywhere under random coefficients", {
  set.seed(56)
  xg <- seq(0.02, 0.98, length.out = 49); zg <- seq(0.02, 0.98, length.out = 25)
  for (sh in list(list(shape_increasing(min_slope = 2), NULL),
                  list(shape_decreasing(min_slope = 1.5), shape_increasing()),
                  list(shape_increasing(min_slope = 2), shape_decreasing(min_slope = 1)),
                  list(shape_increasing(min_slope = 2), shape_convex()))) {
    tm <- term_build(scte(x, z, shapes = sh, smooths = st_sm), st_d)
    p <- length(term_coef_names(tm))
    mx <- sh[[1L]]@shape_params$min_slope * sh[[1L]]@shape_params$sign
    mz <- if (S7::S7_inherits(sh[[2L]], ShapeMonotone)) sh[[2L]]@shape_params$min_slope *
      sh[[2L]]@shape_params$sign else 0
    for (r in 1:5) {
      f <- matrix(term_value(tm, coef = stats::rnorm(p, -1, 1.5), newdata = expand.grid(x = xg, z = zg)),
                  length(xg), length(zg))
      expect_true(all(sign(mx) * diff(f) / diff(xg[1:2]) >= abs(mx) * (1 - 1e-8)))
      if (mz != 0) expect_true(all(sign(mz) * t(diff(t(f))) / diff(zg[1:2]) >= abs(mz) * (1 - 1e-8)))
    }
    expect_identical(qr(cbind(1, term_matrix(term_refresh(tm, stats::rnorm(p, -1, 1)))))$rank, p + 1L)
  }
})

test_that("a value along a line or at a point holds for every coefficient vector", {
  set.seed(57)
  sm <- basis7::bspline_smooth(k = 5, lower = 0, upper = 1)
  gv <- seq(0.02, 0.98, length.out = 25)
  need <- function(sh) {
    if (is.null(sh)) return(c(0, 0))
    p <- sh@shape_params
    if (S7::S7_inherits(sh, ShapeMonotone)) c(p$sign, 0) else c(if (is.null(p$mono)) 0 else p$mono, p$sign)
  }
  cases <- list(
    list(sh = list(shape_increasing(value = c(at = 0.3, value = 1)), NULL), line = c(1, 0.3, 1)),
    list(sh = list(NULL, shape_decreasing(value = c(at = 0.7, value = -1))), line = c(2, 0.7, -1)),
    list(sh = list(shape_increasing(value = c(at = 0, value = 0)), shape_increasing()), line = c(1, 0, 0)),
    list(sh = list(shape_convex(value = c(at = 0.5, value = 2)), shape_concave()), line = c(1, 0.5, 2)),
    list(sh = list(shape_increasing(), shape_convex()), point = c(0.2, 0.6, 1)))
  for (cs in cases) {
    tm <- term_build(scte(x, z, shapes = cs$sh, smooths = sm, point = cs$point), st_d)
    p <- length(term_coef_names(tm))
    nx <- need(cs$sh[[1L]]); nz <- need(cs$sh[[2L]])
    for (r in 1:5) {
      g <- stats::rnorm(p, -1, 1.5)
      f <- matrix(term_value(tm, coef = g, newdata = expand.grid(x = gv, z = gv)), 25, 25)
      tol <- -1e-10 * max(abs(f))
      if (nx[1L] != 0) expect_true(all(nx[1L] * diff(f) >= tol))
      if (nx[2L] != 0) expect_true(all(nx[2L] * diff(f, differences = 2) >= tol))
      if (nz[1L] != 0) expect_true(all(nz[1L] * t(diff(t(f))) >= tol))
      if (nz[2L] != 0) expect_true(all(nz[2L] * t(diff(t(f), differences = 2)) >= tol))
      if (!is.null(cs$line)) {
        nd <- if (cs$line[1L] == 1) data.frame(x = cs$line[2L], z = gv) else data.frame(x = gv, z = cs$line[2L])
        expect_equal(term_value(tm, coef = g, newdata = nd), rep(cs$line[3L], 25), tolerance = 1e-12)
      } else {
        expect_equal(term_value(tm, coef = g, newdata = data.frame(x = cs$point[1L], z = cs$point[2L])),
                     cs$point[3L], tolerance = 1e-12)
      }
    }
    expect_identical(qr(cbind(1, term_matrix(term_refresh(tm, stats::rnorm(p, -1, 1)))))$rank, p + 1L)
  }
  # the block and its derivatives under a line keep the contracts of the block
  tm <- term_build(scte(x, z, shapes = list(shape_convex(value = c(at = 0.5, value = 1)), NULL),
                        smooths = sm), st_d)
  p <- length(term_coef_names(tm)); g <- stats::rnorm(p, sd = 0.7); h <- 1e-5
  e <- function(j) replace(numeric(p), j, h)
  J <- sapply(seq_len(p), function(j) (term_value(tm, coef = g + e(j)) - term_value(tm, coef = g - e(j))) / (2 * h))
  expect_lt(max(abs(term_matrix(term_refresh(tm, g)) - J)) / max(abs(J)), 1e-8)
})

test_that("a decreasing margin is the increasing one in the reflected covariate", {
  set.seed(54)
  dd <- transform(st_d, xn = -x)
  td <- term_build(scte(x, z, shapes = list(st_D, NULL), smooths = st_sm), dd)
  ti <- term_build(scte(xn, z, shapes = list(st_I, NULL), smooths = st_sm), dd)
  g <- stats::rnorm(length(term_coef_names(td)))
  expect_equal(term_value(td, coef = g), term_value(ti, coef = g), tolerance = 1e-10)
})

test_that("scte() names two smoothing parameters and prints its shapes", {
  tm <- term_build(scte(x, z, smooths = st_sm), st_d)
  expect_identical(tm@penalty@params, c("lambda1", "lambda2"))
  expect_identical(term_coef_names(tm)[1L], "scte(x, z).g1")
  expect_output(print(tm), "increasing in x, free in z")
  expect_output(print(scte(x, z, shapes = list(st_I, st_D))), "increasing in x, decreasing in z")
  # beyond the basis the surface is held at the boundary
  g <- stats::rnorm(length(term_coef_names(tm)))
  lo <- tm@build_info$mx$basis@lower
  expect_equal(term_value(tm, coef = g, newdata = data.frame(x = c(lo - 1, lo - 0.1), z = 0.5)),
               rep(term_value(tm, coef = g, newdata = data.frame(x = lo, z = 0.5)), 2))
})

test_that("scte() rejects what it cannot fit", {
  expect_error(scte(x, z, shapes = list(NULL, NULL)), "te() is the free surface", fixed = TRUE)
  expect_error(scte(x, z, shapes = list(shape_unimodal(), NULL)), "curvature shape")
  expect_error(scte(x, z, shapes = list(shape_increasing(min_slope = c(a = 1, b = 2)), NULL)),
               "values by level")
  expect_error(scte(x, z, shapes = list(shape_convex(slope = c(at = 0.5, value = 0)), NULL)),
               "slope at a point")
  expect_error(scte(x, z, shapes = list(shape_increasing(value = c(at = 0.5, value = 0)), NULL),
                    point = c(0.5, 0.5, 0)), "one value")
  expect_error(scte(x, z, point = c(a = 0.5, z = 0.5, value = 0)), "named x, z and value")
  expect_error(term_build(scte(x, z, point = c(2, 0.5, 0), smooths = st_sm), st_d),
               "outside the interval")
  expect_error(scte(x, z, shapes = list(st_I)), "two elements")
  expect_error(scte(x, z, smooths = basis7::fourier_smooth()), "bspline_smooth")
  expect_error(term_build(scte(x, z), data.frame(x = c(NA, runif(9)), z = runif(10))), "missing")
})

test_that("scte() takes a factor or a numeric by", {
  set.seed(55)
  dd <- data.frame(x = runif(180), z = runif(180), g = gl(3, 60, labels = c("a", "b", "c")),
                   w = stats::rnorm(180))
  tf <- term_build(scte(x, z, smooths = st_sm, by = g), dd)
  one <- term_build(scte(x, z, smooths = st_sm), dd)
  k <- length(term_coef_names(one))
  expect_length(term_coef_names(tf), 3L * k)
  expect_identical(term_coef_names(tf)[1L], "scte(x, z).a.g1")
  g <- stats::rnorm(3L * k, sd = 0.6)
  sep <- sapply(1:3, function(l) term_value(one, coef = g[(l - 1L) * k + seq_len(k)]))
  expect_equal(term_value(term_refresh(tf, g)), sep[cbind(seq_len(180), as.integer(dd$g))],
               tolerance = 1e-13)
  tn <- term_build(scte(x, z, shapes = list(st_I, st_I), smooths = st_sm, by = w), dd)
  expect_identical(tn@label, "scte(x, z):w")
  expect_identical(term_coef_names(tn)[1:2], c("scte(x, z):w.const", "scte(x, z):w.g1"))
  gn <- stats::rnorm(length(term_coef_names(tn)), sd = 0.6)
  f1 <- term_value(tn, coef = gn, newdata = transform(dd, w = 1))
  expect_equal(term_value(term_refresh(tn, gn)), dd$w * f1, tolerance = 1e-13)
  # the block and its derivative under a by
  h <- 1e-5; p <- length(gn)
  e <- function(j) replace(numeric(p), j, h)
  J <- sapply(seq_len(p), function(j) (term_value(tn, coef = gn + e(j)) - term_value(tn, coef = gn - e(j))) / (2 * h))
  expect_lt(max(abs(term_matrix(term_refresh(tn, gn)) - J)) / max(abs(J)), 1e-8)
  lv <- term_build(scte(x, z, smooths = st_sm, by = g, by_hyper = "level",
                        hyper = c(lambda2.b = 3)), dd)
  ent <- term_penalties(lv)
  expect_identical(vapply(ent, `[[`, "", "name"), c("a", "b", "c"))
  expect_equal(unlist(ent[[2L]]$fixed), c(lambda2 = 3))
  expect_error(scte(x, z, by_hyper = "level"), "has no 'by'")
  expect_error(term_build(scte(x, z, by = w, by_hyper = "level"), dd), "no factor 'by'")
})
