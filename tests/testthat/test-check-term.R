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
