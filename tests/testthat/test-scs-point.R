# scs() with a value or a slope at a point: the value against the centered
# contribution moved to pass through the point, the slope against an explicit
# curve whose free slope is found by root finding, the block against central
# differences, and the rejections.

set.seed(31)
pt_d <- data.frame(x = runif(120))
pt_sm <- basis7::bspline_smooth(k = 9)
pt_new <- data.frame(x = c(0.2, 0.45, 0.8, 1.3))

# the curve of a plain curvature shape by explicit loops, level 0, from the
# free slope b2 and the second divided differences b3..bK
pt_curve <- function(bp, sg, b2, bcurv, x) {
  K <- length(bp$xi); xi <- bp$xi
  beta <- numeric(K); s <- b2
  beta[2L] <- (xi[2L] - xi[1L]) * s
  for (j in 3:K) {
    s <- s + sg * bcurv[j - 2L] * (xi[j] - xi[j - 2L]) / 2
    beta[j] <- beta[j - 1L] + (xi[j] - xi[j - 1L]) * s
  }
  list(f = drop(basis7::basis_eval(bp$basis, x) %*% beta),
       d = function(a) drop(basis7::basis_deriv(bp$basis, a, order = 1L) %*% beta))
}

test_that("a value at a point is the centered contribution moved through the point", {
  set.seed(32)
  shapes <- list(
    list(shape_increasing(value = c(at = 0.3, value = 1.5)), shape_increasing()),
    list(shape_decreasing(min_slope = 0.2, value = c(value = -1, at = 0.5)),
         shape_decreasing(min_slope = 0.2)),
    list(shape_concave(value = c(0.7, 0)), shape_concave()),
    list(shape_increasing_convex(value = c(at = 0.55, value = 2)),
         shape_increasing_convex()),
    list(shape_unimodal(mode = 0.4, value = c(at = 0.9, value = 3)),
         shape_unimodal(mode = 0.4)),
    list(shape_sigmoid(inflection = 0.5, value = c(at = 0.3, value = 0)),
         shape_sigmoid(inflection = 0.5)))
  for (pair in shapes) {
    tv <- term_build(scs(x, pair[[1L]], pt_sm), pt_d)
    tc <- term_build(scs(x, pair[[2L]], pt_sm), pt_d)
    pv <- pair[[1L]]@shape_params$value
    expect_equal(term_coef_names(tv), term_coef_names(tc))
    for (i in 1:3) {
      g <- rnorm(length(term_coef_names(tv)), 0, 0.7)
      fa <- term_value(tc, coef = g, newdata = data.frame(x = pv$at))
      expect_equal(term_value(tv, coef = g), term_value(tc, coef = g) - fa + pv$value,
                   tolerance = 1e-10)
      expect_equal(term_value(term_refresh(tv, g)), term_value(tv, coef = g))
      expect_equal(term_value(tv, coef = g, newdata = pt_new),
                   term_value(tc, coef = g, newdata = pt_new) - fa + pv$value,
                   tolerance = 1e-10)
      expect_equal(term_value(tv, coef = g, newdata = data.frame(x = pv$at)), pv$value,
                   tolerance = 1e-12)
      expect_equal(term_matrix(term_refresh(tv, g)), term_matrix(term_refresh(tc, g)) -
                     rep(term_predict(term_refresh(tc, g), data.frame(x = pv$at)),
                         each = nrow(pt_d)), tolerance = 1e-10, ignore_attr = TRUE)
    }
  }
})

test_that("a slope at a point removes the free slope and matches an explicit curve", {
  set.seed(33)
  cases <- list(
    list(sh = shape_convex(slope = c(at = 0.3, value = -0.7)), sg = 1),
    list(sh = shape_concave(slope = c(at = 0.9, value = 0.4)), sg = -1),
    list(sh = shape_convex(value = c(at = 0.5, value = 2), slope = c(at = 0.2, value = 1)),
         sg = 1))
  for (cs in cases) {
    tm <- term_build(scs(x, cs$sh, pt_sm), pt_d)
    bp <- tm@build_info
    ps <- cs$sh@shape_params$slope; pv <- cs$sh@shape_params$value
    expect_equal(length(term_coef_names(tm)), length(bp$xi) - 2L)
    for (i in 1:3) {
      g <- rnorm(length(term_coef_names(tm)), 0, 0.7)
      bcurv <- exp(g)
      b2 <- stats::uniroot(function(s) pt_curve(bp, cs$sg, s, bcurv, ps$at)$d(ps$at) -
                             ps$value, c(-1e3, 1e3), tol = 1e-13)$root
      full <- function(x) pt_curve(bp, cs$sg, b2, bcurv, x)$f
      ref <- if (is.null(pv)) full(pt_d$x) - mean(full(pt_d$x))
             else full(pt_d$x) - full(pv$at) + pv$value
      expect_equal(term_value(tm, coef = g), ref, tolerance = 1e-9)
      refn <- if (is.null(pv)) full(pt_new$x[pt_new$x <= bp$basis@upper]) - mean(full(pt_d$x))
              else full(pt_new$x[pt_new$x <= bp$basis@upper]) - full(pv$at) + pv$value
      inside <- pt_new[pt_new$x <= bp$basis@upper, , drop = FALSE]
      expect_equal(term_value(tm, coef = g, newdata = inside), refn, tolerance = 1e-9)
    }
  }
})

test_that("the slope at the point holds and the curvature keeps its sign on a grid", {
  set.seed(34)
  tm <- term_build(scs(x, shape_convex(slope = c(at = 0.6, value = 0.25)), pt_sm), pt_d)
  bp <- tm@build_info
  grid <- data.frame(x = seq(bp$basis@lower, bp$basis@upper, length.out = 2001))
  h <- 1e-5
  for (i in 1:20) {
    g <- rnorm(length(term_coef_names(tm)), 0, 1.5)
    d1 <- (term_value(tm, coef = g, newdata = data.frame(x = 0.6 + h)) -
             term_value(tm, coef = g, newdata = data.frame(x = 0.6 - h))) / (2 * h)
    expect_equal(d1, 0.25, tolerance = 1e-7)
    f <- term_value(tm, coef = g, newdata = grid)
    expect_true(all(diff(f, differences = 2L) >= -1e-10 * max(abs(f))))
  }
})

test_that("the block and term_predict are the Jacobian of term_value with a point constraint", {
  set.seed(35)
  h <- 1e-6
  for (sh in list(shape_increasing(value = c(at = 0.4, value = 1)),
                  shape_convex(value = c(at = 0.5, value = 2), slope = c(at = 0.3, value = -1)),
                  shape_concave(slope = c(at = 0.7, value = 0)))) {
    tm <- term_build(scs(x, sh, pt_sm), pt_d)
    p <- length(term_coef_names(tm))
    g <- rnorm(p, 0, 0.5); r <- term_refresh(tm, g)
    J <- sapply(seq_len(p), function(j) {
      e <- replace(numeric(p), j, h)
      (term_value(tm, coef = g + e) - term_value(tm, coef = g - e)) / (2 * h)
    })
    expect_true(max(abs(term_matrix(r) - J)) / max(abs(J)) < 1e-6)
    Jn <- sapply(seq_len(p), function(j) {
      e <- replace(numeric(p), j, h)
      (term_value(r, coef = g + e, newdata = pt_new) -
         term_value(r, coef = g - e, newdata = pt_new)) / (2 * h)
    })
    expect_true(max(abs(term_predict(r, pt_new) - Jn)) / max(abs(Jn)) < 1e-6)
    expect_equal(length(term_coef_start(tm, target = 2 * pt_d$x)), p)
    expect_true(all(is.finite(term_coef_start(tm, target = 2 * pt_d$x))))
  }
})

test_that("a point constraint is printed and checked", {
  expect_output(print(shape_convex(value = c(0.5, 2), slope = c(at = 0.3, value = -0.7))),
                "f\\(0.5\\) = 2, f'\\(0.3\\) = -0.7")
  expect_equal(shape_increasing(value = c(value = 3, at = 1))@shape_params$value,
               list(at = 1, value = 3))
  expect_error(shape_increasing(slope = c(0, 1)), "unused argument")
  expect_error(shape_increasing_convex(slope = c(0, 1)), "unused argument")
  expect_error(shape_increasing(value = c(at = 0, val = 1)), "two finite numbers")
  expect_error(shape_convex(slope = c(0, NA)), "two finite numbers")
  expect_error(shape_unimodal(value = 1), "two finite numbers")
  expect_error(term_build(scs(x, shape_increasing(value = c(at = 0, value = 0)), pt_sm), pt_d),
               "outside the interval of the basis")
  expect_error(term_build(scs(x, shape_convex(slope = c(at = 2, value = 0)), pt_sm), pt_d),
               "lower' and 'upper'")
  sh <- shape_increasing()
  sh@shape_params$slope <- list(at = 0.5, value = 1)
  expect_error(term_build(scs(x, sh, pt_sm), pt_d), "constrains every slope")
  tm <- term_build(scs(x, shape_increasing(value = c(at = 0, value = 0)),
                       basis7::bspline_smooth(k = 9, lower = 0, upper = 1)), pt_d)
  expect_equal(term_value(tm, newdata = data.frame(x = 0)), 0)
})
