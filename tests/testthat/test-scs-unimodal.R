# shape_unimodal(): the map against cumulative sums written out, one change of
# sign of the slope on a fine grid, the index from mode, the pilot and the
# versions of term_pilot() and term_variants().

set.seed(61)
un_d <- data.frame(x = runif(120))
un_sm <- basis7::bspline_smooth(k = 10)

test_that("the constructor checks mode and index", {
  expect_null(shape_unimodal()@shape_params$at)
  expect_identical(shape_unimodal(index = 4)@shape_params$index, 4L)
  expect_error(shape_unimodal(mode = c(1, 2)), "single finite")
  expect_error(shape_unimodal(index = 0), "positive whole")
  expect_error(shape_unimodal(index = 2.5), "positive whole")
  expect_error(shape_unimodal(link = linkfunctions7::identity_link()), "onto")
})

test_that("shape_map rises to the index and falls after it", {
  xi <- c(0, 0.1, 0.3, 0.5, 0.6, 0.8, 1)
  K <- length(xi)
  set.seed(62)
  for (M in c(1L, 3L, K)) {
    m <- shape_map(shape_unimodal(index = M), xi)
    b <- c(0.2, exp(rnorm(K - 1L)))
    beta <- numeric(K); beta[1L] <- b[1L]
    for (j in 2:K) beta[j] <- beta[j - 1L] + (xi[j] - xi[j - 1L]) * (if (j <= M) b[j] else -b[j])
    expect_equal(drop(m$A %*% b), beta, tolerance = 1e-12)
    expect_identical(as.integer(m$pen), 2:K)
    expect_identical(which.max(beta), M)
  }
  expect_error(shape_map(shape_unimodal(), xi), "not set")
  expect_error(shape_map(shape_unimodal(index = K + 1L), xi), "outside")
})

test_that("random coefficients give a slope that changes sign at most once", {
  set.seed(63)
  xg <- seq(-0.3, 1.3, length.out = 3001)
  for (mode in c(0.2, 0.5, 0.9)) {
    tm <- term_build(scs(x, shape_unimodal(mode = mode), un_sm), un_d)
    p <- length(term_coef_names(tm))
    for (i in 1:25) {
      f <- term_value(tm, coef = rnorm(p, 0, 1.5), newdata = data.frame(x = xg))
      df <- diff(f)
      s <- sign(df[abs(df) > 1e-12 * max(abs(f))])
      expect_lte(sum(diff(s) != 0), 1L)
    }
  }
})

test_that("mode places the index at the nearest Greville abscissa", {
  tm <- term_build(scs(x, shape_unimodal(mode = 0.37), un_sm), un_d)
  bp <- tm@build_info
  expect_identical(bp$shape@shape_params$index, which.min(abs(bp$xi - 0.37)))
  expect_null(term_pilot(tm))
})

test_that("a unimodal shape without a mode has a pilot and five versions", {
  tm <- term_build(scs(x, shape_unimodal(), un_sm), un_d)
  pl <- term_pilot(tm)
  expect_true(S7::S7_inherits(pl, SmoothTerm))
  expect_identical(pl@label, tm@label)
  xi <- tm@build_info$xi
  v <- term_variants(tm, pilot = -(un_d$x - 0.45)^2)
  M0 <- which.min(abs(xi - un_d$x[which.max(-(un_d$x - 0.45)^2)]))
  expect_identical(names(v), sprintf("M = %d", M0 + (-2):2))
  expect_true(all(vapply(v, function(t) !term_is_built(t), logical(1))))
  expect_identical(unname(vapply(v, function(t) t@shape@shape_params$index, 1L)), M0 + (-2):2)
  # near an end the window is cut to the indices that exist
  v1 <- term_variants(tm, pilot = un_d$x)
  expect_identical(names(v1), sprintf("M = %d", (length(xi) - 2L):length(xi)))
  # the other shapes and s() have no pilot and one version
  expect_null(term_pilot(term_build(scs(x, shape_increasing(), un_sm), un_d)))
  expect_length(term_variants(term_build(s(x), un_d)), 1L)
})
