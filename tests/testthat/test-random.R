# Grouped random intercepts and the three effect distributions.

dd <- data.frame(y = rnorm(9), g = factor(rep(c("a", "b", "c"), 3)),
                 x = rnorm(9))

test_that("random intercepts build the indicator block with a gaussian default", {
  built <- term_build(random(~ 1 | g), dd)
  ref <- stats::model.matrix(~ 0 + g, dd)
  expect_equal(unname(term_matrix(built)), unname(ref), ignore_attr = TRUE)
  expect_identical(term_coef_names(built),
                   c("random.a", "random.b", "random.c"))

  pen <- term_penalty(built)
  ref_pen <- penalties7::ridge_penalty(n_coef = 3)
  expect_identical(pen@params, ref_pen@params)
  beta <- c(0.4, -0.1, 1.2)
  expect_identical(penalties7::penalty_value(pen, beta, list(sigma = 1.5)),
                   penalties7::penalty_value(ref_pen, beta, list(sigma = 1.5)))
  expect_true(term_smooth(built))

  res <- check_term(random(~ 1 | g), dd, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
})

test_that("a character grouping is coerced and levels are pinned at build", {
  dc <- data.frame(g = rep(c("u", "v"), 3))
  built <- term_build(random(~ 1 | g), dc)
  expect_identical(term_coef_names(built), c("random.u", "random.v"))
  # prediction keeps the build-time columns on a subset missing a level
  sub <- data.frame(g = "u")
  expect_identical(dim(term_predict(built, sub)), c(1L, 2L))
  # and rejects a level never seen
  expect_error(term_predict(built, data.frame(g = "w")),
               "not present at build time")
})

test_that("a precision structure gives the structured gaussian prior", {
  skip_if_not_installed("parameters7")
  st <- parameters7::log_cholesky(3)
  built <- term_build(random(~ 1 | g, precision = st), dd)
  pen <- term_penalty(built)
  expect_identical(pen@params, st@free_names)
  expect_true(term_smooth(built))
  # edf runs through the structured penalty's own Hessian
  H <- crossprod(term_matrix(built))
  th <- stats::setNames(as.list(rep(0, length(st@free_names))),
                        st@free_names)
  e <- edf(built, coef = c(0.1, -0.2, 0.3), hessian = H, theta = th)
  expect_true(is.finite(e) && e > 0 && e < 3)

  wrong <- parameters7::log_cholesky(2)
  expect_error(term_build(random(~ 1 | g, precision = wrong), dd),
               "dimension 2")
})

test_that("a distribution on the effects is applied coordinatewise", {
  skip_if_not_installed("distributions7")
  d0 <- distributions7::fixed(distributions7::laplace2_distrib(), mu = 0)
  built <- term_build(random(~ 1 | g, distrib = d0, kinks = 0), dd)
  pen <- term_penalty(built)
  expect_identical(pen@params, "lambda")
  expect_false(term_smooth(built))
})

test_that("the reserved and degenerate cases are rejected", {
  expect_error(random(~ x | g), "block-diagonal composition")
  expect_error(random(~ g), "grouping bar")
  expect_error(random(y ~ 1 | g), "one-sided")
  expect_error(random(~ 1 | g, precision = 1, distrib = 2),
               "mutually exclusive")
  expect_error(term_build(random(~ 1 | g), data.frame(g = rep("a", 4))),
               "at least two levels")
})

test_that("the formula interpreter routes random()", {
  out <- interpret_formula(y ~ x + random(~ 1 | g), dd)
  expect_named(out$terms, c("linpar", "random(~1 | g)"))
  built <- term_build(out$terms[["random(~1 | g)"]], dd)
  expect_identical(term_npar(built), 3L)
})
