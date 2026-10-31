# jfp(): the term, its build and prediction, and the penalty it attaches.

jfp_data <- function() {
  set.seed(3)
  data.frame(x1 = rnorm(20), x2 = rnorm(20),
             g = factor(rep(c("a", "b", "c", "d"), 5)))
}

test_that("jfp() builds a block with a Jeffreys penalty and no hyperparameter", {
  d <- jfp_data()
  b <- term_build(jfp(~ x1 + x2, power = 1), d)
  expect_true(S7::S7_inherits(b, JeffreysTerm))
  expect_true(S7::S7_inherits(b, PenalizedTerm))
  expect_identical(b@coef_names, c("jfp.x1", "jfp.x2"))
  expect_equal(unname(as.matrix(b@X)), unname(cbind(d$x1, d$x2)))
  e <- term_penalties(b)
  expect_length(e, 1L)
  expect_identical(e[[1]]$index, 1:2)
  pen <- e[[1]]$penalty
  expect_true(S7::S7_inherits(pen, JeffreysPenalty))
  expect_identical(pen@params, character(0))
  expect_identical(pen@power, 1)
  expect_null(pen@intercept)
  expect_equal(pen@n_coef, 2)
  expect_length(term_hyper(b), 0L)
  expect_true(term_smooth(b))
})

test_that("a factor is coded with its contrasts, and 0 + g with every level", {
  d <- jfp_data()
  b <- term_build(jfp(~ g), d)
  expect_identical(b@coef_names, c("jfp.gb", "jfp.gc", "jfp.gd"))
  b0 <- term_build(jfp(~ 0 + g), d)
  expect_identical(b0@coef_names, paste0("jfp.g", c("a", "b", "c", "d")))
  # a ridge over the same formula gives every indicator
  expect_identical(term_build(ridge(~ g), d)@coef_names,
                   paste0("ridge.g", c("a", "b", "c", "d")))
  nd <- d[c(4, 1, 2), ]
  expect_equal(unname(as.matrix(term_predict(b, nd))),
               unname(as.matrix(b@X)[c(4, 1, 2), ]))
  bs <- term_build(jfp(~ g, sparse = TRUE), d)
  expect_true(methods::is(bs@X, "sparseMatrix"))
  expect_equal(as.matrix(bs@X), as.matrix(b@X), ignore_attr = TRUE)
})

test_that("a matrix input is kept and predicted from the data", {
  d <- jfp_data()
  M <- cbind(u = d$x1, v = d$x2)
  d$M <- M
  b <- term_build(jfp(M), d)
  expect_identical(b@coef_names, c("jfp.u", "jfp.v"))
  expect_equal(unname(as.matrix(term_predict(b, d[1:3, ]))),
               unname(as.matrix(b@X)[1:3, ]))
})

test_that("jfp() rejects what it does not take", {
  d <- jfp_data()
  expect_error(term_build(jfp(~ x1 + I(2 * x1)), d), "linearly dependent")
  expect_error(jfp(~ x1, power = 0), "positive")
  expect_error(jfp(~ x1, power = c(1, 2)), "positive")
  expect_error(jfp(~ x1, power = NULL), "positive")
  expect_error(jfp(~ x1, intercept = NA), "NULL, TRUE or FALSE")
  expect_error(jfp(~ x1, lambda = 1), "has no argument 'lambda'")
})

test_that("the penalty answers what it can and rejects what needs the model", {
  p <- jeffreys_penalty(3, power = 2, intercept = TRUE)
  expect_identical(penalties7::penalty_kinks(p, list()), numeric(0))
  expect_false(penalties7::is_proper(p))
  expect_false(penalties7::is_quadratic(p))
  expect_error(penalties7::penalty_value(p, c(1, 2, 3), list()),
               "evaluated by the fitting")
  expect_error(penalties7::penalty_gradient(p, c(1, 2, 3), list()),
               "evaluated by the fitting")
  expect_error(penalties7::penalty_hessian(p, c(1, 2, 3), list()),
               "evaluated by the fitting")
  expect_error(penalties7::penalty_draw(p, list()), "evaluated by the fitting")
})
