# Repairs from the documentation audit of 0.95.1 (first tranche). Every block
# fails on 0.95.1.

x <- seq(0, 5, length.out = 40)
dd <- data.frame(x = x, g = factor(rep(1:4, 10)), x1 = sin(x), x2 = cos(3 * x))
dd$y <- 2 * exp(-0.7 * x) + 0.01 * dd$x1
fn <- function(x, th) th$a * exp(-th$r * x)
hess <- function(th, d) {
  e <- exp(-th$r * d$x)
  list(a_a = 0 * d$x, a_r = -d$x * e, r_r = th$a * d$x^2 * e)
}
nl_err <- function(b, order) {
  e <- exp(-0.7 * x)
  z <- 0 * x
  ref <- list(`3` = list(a_a_a = z, a_a_r = z, a_r_r = x^2 * e,
                         r_r_r = -2 * x^3 * e),
              `4` = list(a_a_a_a = z, a_a_a_r = z, a_a_r_r = z,
                         a_r_r_r = -x^3 * e, r_r_r_r = 2 * x^4 * e))[[
                as.character(order)]]
  got <- nl_fderiv(b, coef = c(2, 0.7), order = order)
  max(vapply(names(ref), function(k) max(abs(got[[k]] - ref[[k]])), 0))
}

test_that("an opaque nl() function takes one stencil on its value", {
  b <- term_build(nl(fn, params = c("a", "r"), x = x,
                     start = list(a = 2, r = 0.7)), dd)
  # differencing the numerical gradient gave 3.9e-5 at order four
  expect_lt(nl_err(b, 4), 2e-5)
})

test_that("a supplied Hessian is the base of orders three and four", {
  b <- term_build(nl(fn, params = c("a", "r"), x = x, hessian = hess,
                     start = list(a = 2, r = 0.7)), dd)
  expect_lt(nl_err(b, 3), 1e-9)
  expect_lt(nl_err(b, 4), 1e-7)
})

test_that("nl() rejects start names and params it would ignore", {
  expect_error(term_build(nl(~ a * exp(-r * x), start = list(a = 2, rr = 1)),
                          dd), "'rr'")
  expect_error(nl(~ a * exp(-r * x), params = "zzz"), "params")
})

test_that("a supplied derivative of the wrong length is rejected", {
  g3 <- function(th, d) list(a = d$x[1:4], r = d$x[1:4])
  expect_error(term_build(nl(fn, params = c("a", "r"), x = x, gradient = g3,
                             start = list(a = 2, r = 0.7)), dd), "length")
})

test_that("an nl() term checks the length of its coefficients", {
  b <- term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 0.7)), dd)
  expect_error(term_value(b, coef = 1), "length 2")
  expect_error(term_value(b, coef = c(2, 0.7, 9)), "length 2")
  expect_error(term_block_deriv(b, coef = 1, v = c(1, 0)), "length 2")
  expect_error(nl_fderiv(b, coef = 1, order = 1), "length 2")
})

test_that("check_term() rejects a data frame without rows", {
  expect_error(check_term(linpar(~ x), dd[0, ], verbose = FALSE), "no rows")
})

test_that("a penalized term over a formula without columns is rejected", {
  expect_error(term_build(ridge(~ 0), dd), "no columns")
  expect_error(term_build(lasso(~ -1), dd), "no columns")
})

test_that("a smooth's held hyperparameter must be finite", {
  expect_error(s(x, basis7::bspline_smooth(k = 6), hyper = Inf), "hyper")
})

test_that("a censored response with a missing bound gives the validator's message", {
  expect_error(censored_response(y = 1, lwr = NA_real_, upr = 1,
                                 status = "observed"),
               "lwr and upr must not be missing")
  expect_error(cens(factor(c("1", "5", "3"))), "factor")
})

test_that("interpret_formula() passes on the error of a failing constructor", {
  expect_error(interpret_formula(y ~ ridge(~ g, lamda = 1), dd), "lamda")
})

test_that("every penalty entry carries the same fields", {
  want <- c("name", "index", "penalty", "fixed", "n_values", "values",
            "min_ratio", "search", "ids")
  terms <- list(ridge(~ x1),
                s(x, basis7::bspline_smooth(k = 6), by = g, by_hyper = "level"),
                random(~ 1 + x1 | g),
                seg(x, npsi = 2, gamma ~ 0 + lasso(~ 1)))
  for (tm in terms) {
    for (e in term_penalties(term_build(tm, dd))) expect_named(e, want)
  }
  # a written-out grid inside a break-point term's subformula is kept
  b <- term_build(seg(x, npsi = 2, gamma ~ 0 + lasso(~ 1, lambda = c(0.1, 1))),
                  dd)
  expect_equal(unname(unlist(term_values(b))), c(0.1, 1))
})

test_that("a missing group at new data is an error unless a zero row is asked", {
  b <- term_build(random(~ 1 | g), dd)
  nd <- dd[1:3, ]
  nd$g <- factor(c("1", NA, "3"), levels = levels(dd$g))
  expect_error(term_predict(b, nd), "missing")
  Z <- as.matrix(term_predict(b, nd, unseen = "zero"))
  expect_true(all(Z[2L, ] == 0))
})

test_that("a held hyperparameter has no grid, depth or search", {
  expect_identical(names(term_grid(enet(~ x1 + x2, alpha = 0.5))[[1L]]),
                   "lambda")
  expect_identical(term_grid(lasso(~ x1, lambda = 0.3)), list())
  expect_identical(term_path_min(lasso(~ x1, lambda = 0.3)), list())
  expect_identical(term_search(enet(~ x1 + x2, alpha = 0.5)), list())
})

test_that("a marginal seg() carrying the linear effect absorbs the covariate", {
  d2 <- data.frame(x = rep(seq(0, 1, length.out = 10), 4),
                   id = factor(rep(1:4, each = 10)))
  d2$y <- d2$x
  expect_warning(
    it <- interpret_formula(y ~ x + seg(x, psi ~ random(~ 1 | id),
                                        marginal = TRUE), d2),
    "collinear")
  expect_identical(deparse(it$terms$linpar@formula), "~1")
})

test_that("a path argument for a hyperparameter that is not swept is rejected", {
  expect_error(scad(~ x1, n_a = 6), "'n_a' in 'scad' has no effect")
  expect_error(mcp(~ x1, n_gamma = 7), "'n_gamma' in 'mcp' has no effect")
  expect_error(scad(~ x1, search = "cyclic"), "'search' in 'scad'")
  expect_error(lasso(~ x1, lambda = 0.3, n_lambda = 50), "'n_lambda'")
  expect_error(lasso(~ x1, lambda = 0.3, min_ratio = 1e-6), "'min_ratio'")
  expect_error(enet(~ x1, alpha = 0.5, n_alpha = 9), "'n_alpha'")
  expect_error(lasso(~ x1, lambda = c(0.1, 1, 10), n_lambda = 50),
               "written out")
  # the same arguments are accepted where the hyperparameter is swept
  expect_no_error(scad(~ x1, a = NULL, n_a = 6, search = "cyclic"))
  expect_no_error(lasso(~ x1, n_lambda = 50, min_ratio = 1e-6))
  # a written-out grid reports no grid size and no depth
  expect_identical(term_grid(lasso(~ x1, lambda = c(0.1, 1, 10))), list())
  expect_identical(term_path_min(lasso(~ x1, lambda = c(0.1, 1, 10))), list())
})
