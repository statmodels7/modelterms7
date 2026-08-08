# The penalized quartet. Each term's penalty is pinned against the
# penalties7 constructor called directly: the same object, so no tolerance
# has to be chosen.

set.seed(7)
dd <- data.frame(x1 = rnorm(8), x2 = rnorm(8),
                 g = factor(rep(c("a", "b"), 4)))
dd$R <- matrix(rnorm(24), 8, 3, dimnames = list(NULL, c("r1", "r2", "r3")))

test_that("the four terms attach their penalty and read smoothness from it", {
  cases <- list(
    list(ctor = ridge, ref = penalties7::ridge_penalty, smooth = TRUE),
    list(ctor = lasso, ref = penalties7::lasso_penalty, smooth = FALSE),
    list(ctor = scad, ref = penalties7::scad_penalty, smooth = FALSE),
    list(ctor = mcp, ref = penalties7::mcp_penalty, smooth = FALSE)
  )
  beta <- c(0.6, -1.1, 2.0)
  for (case in cases) {
    built <- term_build(case$ctor(~ x1 + x2 + g), dd)
    pen <- term_penalty(built)
    ref <- case$ref(n_coef = term_npar(built))
    expect_identical(pen@params, ref@params)
    th <- lapply(pen@params_bounds, function(b) {
      if (all(is.finite(b))) mean(b) else if (is.finite(b[1])) b[1] + 1 else 1
    })
    expect_identical(penalties7::penalty_value(pen, beta, th),
                     penalties7::penalty_value(ref, beta, th))
    expect_identical(term_smooth(built), case$smooth)
  }
})

test_that("a formula input removes the intercept and keeps the blueprint discipline", {
  built <- term_build(ridge(~ x1 + x2), dd)
  ref <- stats::model.matrix(~ x1 + x2 - 1, dd)
  expect_equal(unname(term_matrix(built)), unname(ref), ignore_attr = TRUE)
  expect_identical(term_coef_names(built), c("ridge.x1", "ridge.x2"))
  res <- check_term(ridge(~ x1 + x2 + g), dd, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
})

test_that("a matrix column of the data builds, predicts and validates", {
  spec <- interpret_formula(~ lasso(R), dd)$terms[["lasso(R)"]]
  built <- term_build(spec, dd)
  expect_identical(term_coef_names(built),
                   c("lasso.r1", "lasso.r2", "lasso.r3"))
  expect_equal(unname(term_matrix(built)), unname(dd$R), ignore_attr = TRUE)

  # prediction re-evaluates the expression in the new data, so a subset of
  # the data carries the matching rows
  sub <- dd[c(2, 5, 7), , drop = FALSE]
  expect_equal(unname(term_predict(built, sub)), unname(dd$R[c(2, 5, 7), ]),
               ignore_attr = TRUE)

  for (lb in c("ridge(R)", "lasso(R)", "scad(R)", "mcp(R)")) {
    f <- stats::as.formula(paste("~", lb))
    sp <- interpret_formula(f, dd)$terms[[lb]]
    res <- check_term(sp, dd, verbose = FALSE)
    expect_true(all(res$status == "OK"),
                info = paste(lb, paste(res$check[res$status != "OK"],
                                       collapse = ", ")))
  }
})

test_that("a free-standing matrix builds but refuses a mismatched prediction", {
  M <- matrix(rnorm(16), 8, 2)
  built <- term_build(ridge(M), dd)
  expect_identical(term_coef_names(built), c("ridge.1", "ridge.2"))
  # the expression 'M' does not resolve in newdata, and the build-time rows
  # are deliberately not reused
  expect_error(term_predict(built, dd[1:3, , drop = FALSE]),
               "column of the data")
  # and a matrix whose rows do not match the data cannot build at all
  expect_error(term_build(ridge(M[1:5, ]), dd), "5 rows")
})

test_that("by is reserved and inputs are validated", {
  expect_error(ridge(~x1, by = dd$g), "reserved")
  expect_error(lasso(y ~ x1), "one-sided")
  expect_error(scad(matrix("a", 2, 2)), "numeric")
  expect_error(mcp(~x1, label = ""), "non-empty")
})

test_that("the formula interpreter routes the penalized terms", {
  out <- interpret_formula(x1 ~ x2 + ridge(~g, label = "rg") + lasso(R), dd)
  expect_named(out$terms,
               c("linpar", 'ridge(~g, label = "rg")', "lasso(R)"))
  built <- term_build(out$terms[["lasso(R)"]], dd)
  expect_identical(term_npar(built), 3L)
})
