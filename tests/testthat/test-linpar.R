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

test_that("a block is built sparse where it is asked for, not compressed after", {
  # a factor of many levels has one non-zero per row, and the DENSE model
  # matrix of it is the memory the choice exists to avoid: building dense and
  # compressing would pay that memory anyway
  set.seed(31)
  n <- 800L
  m <- 120L
  d <- data.frame(g = factor(sample.int(m, n, TRUE)), z = stats::rnorm(n))

  a <- term_build(linpar(~ 0 + g + z), d)
  b <- term_build(linpar(~ 0 + g + z, sparse = TRUE), d)
  expect_true(is.matrix(term_matrix(a)))
  expect_true(methods::is(term_matrix(b), "sparseMatrix"))
  # the same design, only stored differently: same columns, same names,
  # same numbers
  expect_identical(term_coef_names(b), term_coef_names(a))
  expect_equal(unname(as.matrix(term_matrix(b))), unname(term_matrix(a)))
  expect_lt(as.numeric(utils::object.size(term_matrix(b))),
            as.numeric(utils::object.size(term_matrix(a))) / 5)

  # the STORAGE is part of the blueprint, so prediction does not densify
  nd <- d[1:20, , drop = FALSE]
  expect_true(methods::is(term_predict(b, nd), "sparseMatrix"))
  expect_equal(unname(as.matrix(term_predict(b, nd))),
               unname(term_predict(a, nd)))

  expect_error(linpar(~z, sparse = "yes"), "must be TRUE or FALSE")
})

test_that("a penalized term's FORMULA route takes the same storage", {
  # a matrix input needs no argument, being kept in whatever storage it
  # arrives in; what would otherwise be dense whatever the caller passed is
  # the model matrix of a formula
  set.seed(32)
  n <- 400L
  d <- data.frame(g = factor(sample.int(60L, n, TRUE)),
                  y = stats::rnorm(n))
  for (ctor in list(ridge, lasso, enet, scad, mcp)) {
    b <- term_build(ctor(~ 0 + g, sparse = TRUE), d)
    expect_true(methods::is(term_matrix(b), "sparseMatrix"))
    expect_true(methods::is(term_predict(b, d[1:5, , drop = FALSE]),
                            "sparseMatrix"))
  }
  # and standardization does not undo it: it is a diagonal map on the
  # PENALTY and never an operation on the design
  s <- term_build(lasso(~ 0 + g, sparse = TRUE, standardize = TRUE), d)
  expect_true(methods::is(term_matrix(s), "sparseMatrix"))
})

test_that("linpar takes the contrasts for its factors", {
  d <- data.frame(g = factor(rep(letters[1:3], 4)), y = stats::rnorm(12))
  a <- term_build(linpar(~ g), d)
  b <- term_build(linpar(~ g, contrasts = list(g = "contr.sum")), d)
  # the same rank, a different coding, and the blueprint carries it to
  # prediction rather than letting the session's option decide again
  expect_identical(ncol(term_matrix(a)), ncol(term_matrix(b)))
  expect_false(isTRUE(all.equal(unname(term_matrix(a)),
                                unname(term_matrix(b)))))
  expect_equal(unname(term_predict(b, d[1:3, , drop = FALSE])),
               unname(term_matrix(b)[1:3, , drop = FALSE]))
  expect_error(linpar(~g, contrasts = "contr.sum"), "named list")
})
