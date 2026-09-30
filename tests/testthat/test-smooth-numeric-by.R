test_that("a numeric by keeps the constant: the block spans z and z x", {
  # The term is f(x) z. Its constant is the coefficient of z itself, so the
  # block must contain z; up to modelterms7 0.85.0 it was centered and did
  # not, and a fit missed the main effect of z unless the formula wrote it.
  set.seed(11)
  d <- data.frame(x = runif(200), z = rnorm(200))
  b <- term_build(s(x, basis7::bspline_smooth(k = 10), by = z), d)
  X <- term_matrix(b)
  expect_identical(ncol(X), 10L)
  expect_identical(term_coef_names(b)[1:2], c("s(x):z.const", "s(x):z.lin"))
  expect_equal(unname(X[, 1]), d$z, tolerance = 1e-12)
  Q <- cbind(d$z, d$z * d$x)
  res <- Q - X %*% qr.solve(X, Q)
  expect_lt(max(abs(res)), 1e-10)
  # the constant and the linear part are not penalized
  S <- as.matrix(penalties7::penalty_hessian(b@penalty, rep(0, 10),
                                             list(lambda = 1)))
  expect_lt(max(abs(S[1:2, ])), 1e-12)
  expect_gt(max(abs(S[3:10, 3:10])), 0.5)
  # and the block is reapplied at new rows, constant included
  expect_equal(term_predict(b, d[1:7, ]), X[1:7, ], tolerance = 1e-12,
               ignore_attr = TRUE)
})

test_that("a factor by stays centered", {
  set.seed(12)
  d <- data.frame(x = runif(120), g = factor(rep(c("a", "b"), 60)))
  b <- term_build(s(x, basis7::bspline_smooth(k = 6), by = g), d)
  expect_false(any(grepl("const", term_coef_names(b))))
  expect_identical(ncol(term_matrix(b)), 10L)
})

test_that("a numeric by names its term after the by variable", {
  set.seed(13)
  d <- data.frame(x = runif(80), z = rnorm(80), w = runif(80))
  expect_identical(term_build(s(x, by = z), d)@label, "s(x):z")
  expect_identical(term_build(s(x, by = z, label = "vc"), d)@label, "vc")
  expect_identical(term_build(s(x), d)@label, "s(x)")
  tb <- term_build(te(x, w, by = z), d)
  expect_identical(tb@label, "te(x,w):z")
  expect_identical(term_coef_names(tb)[1], "te(x,w):z.const")
})
