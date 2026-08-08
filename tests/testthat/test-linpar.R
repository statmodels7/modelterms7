# The unpenalized parametric term: model.matrix semantics and the blueprint.

test_that("linpar builds the model.matrix block", {
  dd <- data.frame(x = 1:6, g = factor(rep(c("a", "b", "c"), 2)))
  built <- term_build(linpar(~ x + g), dd)

  ref <- stats::model.matrix(~ x + g, dd)
  expect_equal(unname(term_matrix(built)), unname(ref), ignore_attr = TRUE)
  expect_identical(term_coef_names(built), colnames(ref))
  expect_identical(term_npar(built), ncol(ref))
  expect_true(term_smooth(built))
  expect_null(term_penalty(built))
})

test_that("a label prefixes the coefficient names", {
  dd <- data.frame(x = 1:4)
  built <- term_build(linpar(~x, label = "base"), dd)
  expect_identical(term_coef_names(built), c("base.(Intercept)", "base.x"))
  expect_identical(colnames(term_matrix(built)), term_coef_names(built))
})

test_that("the blueprint reproduces the mapping on new data", {
  dd <- data.frame(x = 1:6, g = factor(rep(c("a", "b", "c"), 2)))
  built <- term_build(linpar(~ x + g), dd)

  # the same data reproduce the block exactly
  expect_identical(unname(term_predict(built, dd)),
                   unname(term_matrix(built)))

  # a subset that drops a whole level keeps the encoding of build time
  sub <- dd[dd$g != "c", , drop = FALSE]
  expect_identical(unname(term_predict(built, sub)),
                   unname(term_matrix(built)[dd$g != "c", , drop = FALSE]))

  # a level unknown to the blueprint is rejected, not re-encoded
  bad <- data.frame(x = 1, g = factor("d"))
  expect_error(term_predict(built, bad), "factor")
})

test_that("missing values propagate instead of dropping rows", {
  dd <- data.frame(x = c(1, NA, 3))
  built <- term_build(linpar(~x), dd)
  expect_identical(nrow(term_matrix(built)), 3L)
  expect_true(is.na(term_matrix(built)[2L, "x"]))
})

test_that("specifications and built terms are told apart", {
  spec <- linpar(~x)
  expect_false(term_is_built(spec))
  expect_error(term_matrix(spec), "not been built")
  expect_error(term_npar(spec), "not been built")
  expect_error(term_predict(spec, data.frame(x = 1)), "not been built")

  built <- term_build(spec, data.frame(x = 1:3))
  expect_true(term_is_built(built))
  # building does not mutate the specification
  expect_false(term_is_built(spec))
})

test_that("constructor and generic inputs are validated", {
  expect_error(linpar(y ~ x), "one-sided")
  expect_error(linpar("~x"), "formula")
  expect_error(linpar(~x, label = c("a", "b")), "single")
  expect_error(term_build(linpar(~x), "not a data frame"), "data frame")
})

test_that("print says what an object is", {
  expect_output(print(linpar(~x)), "specification")
  built <- term_build(linpar(~x), data.frame(x = 1:3))
  expect_output(print(built), "2 coefficients")
})
