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
  for (sh in list(list(st_I, NULL), list(NULL, st_D), list(st_I, st_I), list(st_D, st_I))) {
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
  expect_error(scte(x, z, shapes = list(shape_convex(), NULL)), "shape_increasing")
  expect_error(scte(x, z, shapes = list(shape_increasing(min_slope = 1), NULL)), "min_slope")
  expect_error(scte(x, z, shapes = list(st_I)), "two elements")
  expect_error(scte(x, z, smooths = basis7::fourier_smooth()), "bspline_smooth")
  expect_error(term_build(scte(x, z), data.frame(x = c(NA, runif(9)), z = runif(10))), "missing")
})
