test_that("a linpar() beside the formula's intercept gives its own up", {
  dd <- data.frame(y = rnorm(24), x = rnorm(24),
                   g = factor(rep(c("a", "b", "c"), 8)))
  it <- interpret_formula(y ~ linpar(~ x * g), dd)
  expect_identical(names(it$terms)[1], "linpar")
  ex <- it$terms[[2]]
  expect_true(ex@drop_intercept)
  bt <- term_build(ex, dd)
  ref <- stats::model.matrix(~ x * g, dd)[, -1]
  expect_identical(term_coef_names(bt), colnames(ref))
  expect_equal(unname(as.matrix(term_matrix(bt))), unname(ref))
  # the coding is the one with an intercept, not full indicators
  expect_false(any(grepl("^ga$", term_coef_names(bt))))
  # prediction reapplies the removal
  nd <- dd[c(2, 5, 9), ]
  expect_equal(unname(as.matrix(term_predict(bt, nd))), unname(ref[c(2, 5, 9), ]))
})

test_that("the column survives under 0 +, and a labelled or sparse block drops it too", {
  dd <- data.frame(y = rnorm(24), x = rnorm(24),
                   g = factor(rep(c("a", "b", "c"), 8)))
  keep <- interpret_formula(y ~ 0 + linpar(~ x + g), dd)
  expect_identical(names(keep$terms), "linpar(~x + g)")
  expect_false(keep$terms[[1]]@drop_intercept)
  expect_true("(Intercept)" %in% term_coef_names(term_build(keep$terms[[1]], dd)))
  lab <- interpret_formula(y ~ linpar(~ x + g, label = "a", sparse = TRUE), dd)
  bl <- term_build(lab$terms[[2]], dd)
  expect_identical(term_coef_names(bl), c("a.x", "a.gb", "a.gc"))
  expect_s4_class(term_matrix(bl), "dgCMatrix")
  # a block of nothing but the intercept is left alone
  one <- interpret_formula(y ~ x + linpar(~ 1), dd)
  expect_false(one$terms[[2]]@drop_intercept)
})
