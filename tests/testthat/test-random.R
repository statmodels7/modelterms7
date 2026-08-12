# Grouped random effects: intercepts, slopes, and the effect distributions.

set.seed(21)
dd <- data.frame(y = rnorm(12), x = rnorm(12),
                 g = factor(rep(c("a", "b", "c"), 4)))

test_that("random intercepts build the indicator block with a gaussian default", {
  built <- term_build(random(~ 1 | g), dd)
  ref <- stats::model.matrix(~ 0 + g, dd)
  # the block is SPARSE by construction: a row belongs to one group, so the
  # density is 1/m whatever the data, and the dense form was the whole cost
  # of a random effect. It is compared against the dense reference by value
  Z <- term_matrix(built)
  expect_s4_class(Z, "sparseMatrix")
  expect_equal(unname(as.matrix(Z)), unname(ref), ignore_attr = TRUE)
  expect_identical(term_coef_names(built),
                   c("random.a", "random.b", "random.c"))

  pen <- term_penalty(built)
  expect_identical(pen@params, "sigma")
  expect_true(term_smooth(built))

  res <- check_term(random(~ 1 | g), dd, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
})

test_that("random slopes interact the within-group design with the groups", {
  built <- term_build(random(~ x | g), dd)
  expect_identical(term_npar(built), 6L)
  expect_identical(term_coef_names(built),
                   c("random.a.(Intercept)", "random.a.x",
                     "random.b.(Intercept)", "random.b.x",
                     "random.c.(Intercept)", "random.c.x"))

  # reference construction: indicator times within-group column, group-major
  Z <- as.matrix(term_matrix(built))
  G <- stats::model.matrix(~ 0 + g, dd)
  expect_equal(unname(Z[, 1]), unname(G[, 1]), ignore_attr = TRUE)
  expect_equal(unname(Z[, 2]), unname(G[, 1] * dd$x), ignore_attr = TRUE)
  expect_equal(unname(Z[, 6]), unname(G[, 3] * dd$x), ignore_attr = TRUE)
  # the block reproduces the fixed design when the group effects are summed
  expect_equal(rowSums(Z[, c(1, 3, 5)]), rep(1, 12))
  expect_equal(rowSums(Z[, c(2, 4, 6)]), dd$x)

  res <- check_term(random(~ x | g), dd, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))

  # and the slope alone drops the intercept by the formula convention
  b0 <- term_build(random(~ 0 + x | g), dd)
  expect_identical(term_npar(b0), 3L)
  expect_identical(term_coef_names(b0),
                   c("random.a.x", "random.b.x", "random.c.x"))
})

test_that("the default gaussian is unstructured or diagonal by 'correlated'", {
  bc <- term_build(random(~ x | g), dd)
  ref_c <- penalties7::structured_penalty(
    parameters7::kron_identity(parameters7::log_cholesky(2), 3))
  expect_identical(term_penalty(bc)@params, ref_c@params)

  bu <- term_build(random(~ x | g, correlated = FALSE), dd)
  ref_u <- penalties7::structured_penalty(
    parameters7::kron_identity(parameters7::diagonal_matrix(2), 3))
  expect_identical(term_penalty(bu)@params, ref_u@params)
  expect_identical(length(term_penalty(bc)@params), 3L)
  expect_identical(length(term_penalty(bu)@params), 2L)

  # both are the same penalty penalties7 builds directly
  beta <- rnorm(6)
  th <- stats::setNames(as.list(rep(0.1, 3)), ref_c@params)
  expect_identical(penalties7::penalty_value(term_penalty(bc), beta, th),
                   penalties7::penalty_value(ref_c, beta, th))
  expect_true(term_smooth(bc))

  # edf runs on the kron-structured block
  H <- crossprod(term_matrix(bc))
  e <- edf(bc, coef = beta, hessian = H, theta = th)
  expect_true(is.finite(e) && e > 0 && e < 6)
})

test_that("a per-group precision structure replaces the default", {
  st <- parameters7::compound_symmetry(2)
  built <- term_build(random(~ x | g, precision = st), dd)
  expect_identical(term_penalty(built)@params, st@free_names)

  wrong <- parameters7::log_cholesky(3)
  expect_error(term_build(random(~ x | g, precision = wrong), dd),
               "dimension 3")
  expect_error(term_build(random(~ 1 | g, precision = diag(2)), dd),
               "matrix_parameter")
})

test_that("a distribution on the effects is applied coordinatewise", {
  d0 <- distributions7::fixed(distributions7::laplace2_distrib(), mu = 0)
  built <- term_build(random(~ x | g, distrib = d0, kinks = 0), dd)
  pen <- term_penalty(built)
  expect_identical(pen@params, "lambda")
  expect_identical(pen@n_coef, 6L)
  expect_false(term_smooth(built))
})

test_that("levels are pinned at build and unknown levels are rejected", {
  dc <- data.frame(x = rnorm(6), g = rep(c("u", "v"), 3))
  built <- term_build(random(~ x | g), dc)
  sub <- data.frame(x = 0.5, g = "u")
  expect_identical(dim(term_predict(built, sub)), c(1L, 4L))
  expect_error(term_predict(built, data.frame(x = 1, g = "w")),
               "not present at build time")
})

test_that("degenerate inputs are rejected", {
  expect_error(random(~ g), "grouping bar")
  expect_error(random(y ~ 1 | g), "one-sided")
  expect_error(random(~ 1 | g, correlated = NA), "TRUE or FALSE")
  expect_error(random(~ 1 | g, precision = 1, distrib = 2),
               "mutually exclusive")
  expect_error(term_build(random(~ 1 | g), data.frame(g = rep("a", 4))),
               "at least two levels")
})

test_that("the formula interpreter routes random()", {
  out <- interpret_formula(y ~ x + random(~ x | g), dd)
  expect_named(out$terms, c("linpar", "random(~x | g)"))
  built <- term_build(out$terms[["random(~x | g)"]], dd)
  expect_identical(term_npar(built), 6L)
})
