# check_term(): the validator, with paired injections. A correct term must
# pass every check, and each deliberate defect must be caught by the check
# that exists for it -- otherwise the check could be trivially green.

dd_factor <- data.frame(x = 1:6, g = factor(rep(c("a", "b", "c"), 2)))

test_that("a correct term passes every check", {
  res <- check_term(linpar(~ x + g), dd_factor, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
  expect_identical(res$check,
                   c("build", "names", "npar", "smooth", "reproduce", "subset"))
})

# A term whose predict re-derives the factor encoding from the new data
# instead of reusing the blueprint: the classic mistake the subset check
# exists for.
NoBlueprint <- S7::new_class("NoBlueprint", parent = additive_term,
                             package = NULL,
                             properties = list(formula = S7::class_any))
S7::method(term_build, NoBlueprint) <- function(term, data, ...) {
  X <- stats::model.matrix(term@formula, data)
  term@X <- X
  term@coef_names <- colnames(X)
  term@blueprint <- list(formula = term@formula)
  term
}
S7::method(term_predict, NoBlueprint) <- function(term, newdata, ...) {
  stats::model.matrix(term@blueprint$formula, newdata)
}

test_that("a predict that rebuilds instead of reapplying fails the subset check", {
  bad <- NoBlueprint(label = "", formula = ~ x + g, X = NULL,
                     coef_names = character(0), blueprint = list(),
                     penalty = NULL)
  res <- check_term(bad, dd_factor, verbose = FALSE)
  expect_identical(res$status[res$check == "subset"], "FAILED")
  # and the same defect is invisible on data without factors, which is why
  # the subset check drops a level whenever it can
  bad_num <- NoBlueprint(label = "", formula = ~x, X = NULL,
                         coef_names = character(0), blueprint = list(),
                         penalty = NULL)
  res_num <- check_term(bad_num, data.frame(x = 1:6, g2 = rnorm(6)),
                        verbose = FALSE)
  expect_true(all(res_num$status == "OK"))
})

# A predict wrong by 5 percent must fail the reproduce check.
Skewed <- S7::new_class("Skewed", parent = additive_term, package = NULL,
                        properties = list(formula = S7::class_any))
S7::method(term_build, Skewed) <- function(term, data, ...) {
  X <- stats::model.matrix(term@formula, data)
  term@X <- X
  term@coef_names <- colnames(X)
  term@blueprint <- list(formula = term@formula)
  term
}
S7::method(term_predict, Skewed) <- function(term, newdata, ...) {
  1.05 * stats::model.matrix(term@blueprint$formula, newdata)
}

test_that("a predict wrong by 5 percent fails the reproduce check", {
  bad <- Skewed(label = "", formula = ~x, X = NULL,
                coef_names = character(0), blueprint = list(),
                penalty = NULL)
  res <- check_term(bad, data.frame(x = 1:6), verbose = FALSE)
  expect_identical(res$status[res$check == "reproduce"], "FAILED")
})

test_that("a term that does not build reports the build failure and stops", {
  bad <- linpar(~ not_a_column)
  res <- check_term(bad, dd_factor, verbose = FALSE)
  expect_identical(res$status[res$check == "build"], "FAILED")
  expect_identical(nrow(res), 1L)
})

test_that("verbose printing shows one line per check", {
  expect_output(check_term(linpar(~x), data.frame(x = 1:4)), "\\[OK\\] build")
})
test_that("a structural term is rejected by name", {
  # The second statement of the body reads term_matrix(built), which is
  # registered on additive_term alone, so a structural term built successfully
  # and then died inside an internal generic with no row at all.
  dd <- data.frame(y = rnorm(20), x = rnorm(20))

  expect_error(check_term(gas(p = 1, q = 1), dd), "is structural")
  expect_error(check_term(gas(p = 1, q = 1), dd), "GasTerm")
  expect_error(check_term(regime(k = 2), dd), "is structural")

  # the message names what such a term supplies instead
  expect_error(check_term(gas(p = 1, q = 1), dd), "term_filter")

  # and it is not S7's method-not-found error, which is what it was
  msg <- tryCatch(check_term(gas(p = 1, q = 1), dd),
                  error = function(e) conditionMessage(e))
  expect_false(grepl("Can't find method", msg, fixed = TRUE))

  # the additive branch is untouched: every kind still runs its six checks
  dd2 <- data.frame(y = rnorm(20), x = rnorm(20), g = factor(rep(1:2, 10)))
  for (tm in list(linpar(~ x + g), s(x, k = 5), ridge(~x))) {
    res <- check_term(tm, dd2, verbose = FALSE)
    expect_identical(nrow(res), 6L)
    expect_true(all(res$status == "OK"),
                info = paste(res$check[res$status != "OK"], collapse = ", "))
  }
})
