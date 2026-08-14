# edf and the fitted-point displays.

set.seed(11)
dd <- data.frame(x1 = rnorm(30), x2 = rnorm(30), x3 = rnorm(30))

test_that("an unpenalized term counts its coefficients exactly", {
  built <- term_build(linpar(~ x1 + x2), dd)
  expect_identical(edf(built), 3)
})

test_that("a ridge term's edf is the eigenvalue shrinkage of its block", {
  built <- term_build(ridge(~ x1 + x2 + x3), dd)
  H <- crossprod(term_matrix(built))
  beta <- c(0.4, -1, 2)
  sigma <- 1.7

  got <- edf(built, coef = beta, hessian = H, theta = list(sigma = sigma))
  # independent route: the gaussian penalty's Hessian is I/sigma^2, so the
  # trace is sum d_i / (d_i + 1/sigma^2) over the eigenvalues of H
  d <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  expect_equal(got, sum(d / (d + 1 / sigma^2)), tolerance = 1e-12)

  # the limits: no penalty recovers the count, a hard one removes it
  expect_equal(edf(built, beta, H, list(sigma = 1e6)), 3, tolerance = 1e-6)
  expect_lt(edf(built, beta, H, list(sigma = 1e-4)), 1e-4)
})

test_that("a non-smooth term counts its nonzero coefficients", {
  built <- term_build(lasso(~ x1 + x2 + x3), dd)
  expect_identical(edf(built, coef = c(0.5, 0, -2)), 2)
  expect_identical(edf(built, coef = c(0.5, 1e-12, -2)), 2)
  expect_identical(edf(built, coef = c(0, 0, 0)), 0)
  built_scad <- term_build(scad(~ x1 + x2), dd)
  expect_identical(edf(built_scad, coef = c(1, 0)), 1)
})

test_that("edf counts a partially penalized term parameter by parameter", {
  set.seed(12)
  dx <- data.frame(x = sort(runif(200, 0, 10)))
  dx$y <- 1 + 0.5 * dx$x + 2 * pmax(dx$x - 6, 0) + rnorm(200, sd = 0.3)

  # seg carries a lasso on its one slope change and nothing on its linear
  # effect or its break-point, so those two count exactly
  bl <- term_build(seg(x, gamma ~ 0 + lasso(~1)), dx)
  expect_identical(edf(bl, coef = c(0.5, 2, 6)), 3)
  expect_identical(edf(bl, coef = c(0.5, 0, 6)), 2)
  # nothing is asked of the curvature: the count is read from coef alone
  expect_identical(edf(bl, coef = c(0.5, 2, 6), hessian = NULL), 3)

  # under a ridge the same two count exactly and the change is shrunk, so
  # the answer sits between the two limits and reaches them
  br <- term_build(seg(x, gamma ~ 0 + ridge(~1)), dx)
  H <- crossprod(term_matrix(br))
  expect_equal(edf(br, coef = c(0.5, 2, 6), hessian = H,
                   theta = list(sigma = 1e6)), 3, tolerance = 1e-6)
  expect_equal(edf(br, coef = c(0.5, 2, 6), hessian = H,
                   theta = list(sigma = 1e-8)), 2, tolerance = 1e-6)

  # and the trace is the one an independent assembly gives over the two
  # unpenalized coordinates and the penalized one together
  sig <- 1.3
  S <- matrix(0, 3, 3)
  S[2, 2] <- 1 / sig^2
  expect_equal(edf(br, coef = c(0.5, 2, 6), hessian = H,
                   theta = list(sigma = sig)),
               sum(diag(solve(H + S, H))), tolerance = 1e-12)
})

test_that("a term carrying two penalties keys its hyperparameters by name", {
  set.seed(13)
  dx <- data.frame(x = sort(runif(200, 0, 10)))
  bj <- term_build(jseg(x, gamma ~ 0 + ridge(~1), delta ~ 0 + ridge(~1)),
                   dx)
  H <- crossprod(term_matrix(bj))
  cf <- c(0.5, 2, 1, -6)
  ent <- term_penalties(bj)
  expect_length(ent, 2L)

  th <- stats::setNames(list(list(sigma = 1.1), list(sigma = 0.4)),
                        vapply(ent, function(e) e$name, character(1)))
  got <- edf(bj, coef = cf, hessian = H, theta = th)
  S <- matrix(0, 4, 4)
  S[ent[[1L]]$index, ent[[1L]]$index] <- 1 / 1.1^2
  S[ent[[2L]]$index, ent[[2L]]$index] <- 1 / 0.4^2
  expect_equal(got, sum(diag(solve(H + S, H))), tolerance = 1e-12)

  # a single list cannot say which penalty it belongs to, and is refused
  expect_error(edf(bj, coef = cf, hessian = H, theta = list(sigma = 1.1)),
               "keyed by the penalty names")
})

test_that("edf refuses what it cannot compute and says what is missing", {
  expect_error(edf(linpar(~x1)), "not been built")
  br <- term_build(ridge(~ x1 + x2), dd)
  expect_error(edf(br, coef = c(1, 2)), "'hessian'")
  expect_error(edf(br, coef = 1, hessian = diag(2), theta = list(sigma = 1)),
               "length 2")
  expect_error(edf(br, coef = c(1, 2), hessian = diag(3),
                   theta = list(sigma = 1)), "2 x 2")
  bl <- term_build(lasso(~x1), dd)
  expect_error(edf(bl), "'coef'")
})

test_that("print shows the penalty of a built penalized term", {
  expect_output(print(lasso(~x1)), "specification")
  built <- term_build(lasso(~ x1 + x2), dd)
  expect_output(print(built),
                "2 coefficients; penalty separable \\[fixed laplace2 \\[mu=0\\]\\] \\(lambda\\)")
  bs <- term_build(scad(~x1), dd)
  expect_output(print(bs), "penalty SCAD \\(lambda, a\\)")
})

test_that("plot draws a built term at supplied coefficients", {
  built <- term_build(ridge(~ x1 + x2), dd)
  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  on.exit({
    grDevices::dev.off()
    unlink(tf)
  })
  expect_silent(plot(built, coef = c(0.5, -1)))
  expect_error(plot(built), "'coef' is required")
  expect_error(plot(built, coef = 1), "length 2")
  expect_error(plot(linpar(~x1)), "not been built")
})
