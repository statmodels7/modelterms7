# The formula interpreter: recognition by evaluation, the collected
# parametric block, and the left-hand side.

# A user-defined term, defined here exactly as a user would define one:
# recognition works because the constructor returns a model_term, not
# because the interpreter knows the name.
MatTerm <- S7::new_class("MatTerm", parent = additive_term, package = NULL,
                         properties = list(mat = S7::class_any))
mat_term <- function(m, label = "mat") {
  MatTerm(label = label, mat = as.matrix(m),
          X = NULL, coef_names = character(0),
          blueprint = list(), penalty = NULL)
}
S7::method(term_build, MatTerm) <- function(term, data, ...) {
  X <- term@mat
  colnames(X) <- paste(term@label, seq_len(ncol(X)), sep = ".")
  term@X <- X
  term@coef_names <- colnames(X)
  term@blueprint <- list(ncol = ncol(X))
  term
}

test_that("calls returning terms are routed as terms, everything else stays a covariate", {
  set.seed(1)
  dd <- data.frame(y = rnorm(6), x1 = 1:6, x2 = runif(6),
                   g = factor(rep(c("a", "b"), 3)))
  R <- matrix(rnorm(12), 6, 2)

  out <- interpret_formula(y ~ x1 + log(x2) + x1:x2 + mat_term(R), dd)

  expect_named(out$terms, c("linpar", "mat_term(R)"))
  expect_s3_class(out$terms$linpar, "modelterms7::LinparTerm")
  expect_true(S7::S7_inherits(out$terms[["mat_term(R)"]], MatTerm))

  # the collected block carries the covariates, the transformation and the
  # interaction, with the usual conventions
  built <- term_build(out$terms$linpar, dd)
  ref <- stats::model.matrix(~ x1 + log(x2) + x1:x2, dd)
  expect_equal(unname(term_matrix(built)), unname(ref),
               ignore_attr = TRUE)

  # the response is the evaluated left-hand side
  expect_identical(out$response, dd$y)
})

test_that("the intercept convention is the formula's", {
  dd <- data.frame(y = 1:4, x = 4:1)
  R <- diag(4)

  # specials only: an intercept-only parametric block survives
  out <- interpret_formula(y ~ mat_term(R), dd)
  expect_named(out$terms, c("linpar", "mat_term(R)"))
  b <- term_build(out$terms$linpar, dd)
  expect_identical(term_coef_names(b), "(Intercept)")

  # and -1 removes it entirely
  out0 <- interpret_formula(y ~ mat_term(R) - 1, dd)
  expect_named(out0$terms, "mat_term(R)")
  expect_false(out0$intercept)
})

test_that("a censored left-hand side evaluates to its response object", {
  dd <- data.frame(y = c(0, 0.4, 1.9), x = 1:3)
  out <- interpret_formula(cens(y, lwr = 0) ~ x, dd)
  expect_true(S7::S7_inherits(out$response, censored_response))
  expect_identical(out$response@status, c("left", "observed", "observed"))
})

test_that("a one-sided formula has no response", {
  dd <- data.frame(x = 1:3)
  out <- interpret_formula(~x, dd)
  expect_null(out$response)
  expect_named(out$terms, "linpar")
})

test_that("inputs are validated", {
  expect_error(interpret_formula("y ~ x", data.frame(x = 1)), "formula")
  expect_error(interpret_formula(~x, "not a data frame"), "data frame")
})

test_that("a structural class without term_build is named, like an additive one", {
  # There is one default and it names the class. The structural default it
  # replaces said structural terms were reserved for a later release, which
  # gas(), regime() and the marginal break-point terms had made false, and
  # which only a class written outside the package could ever have read.
  Gassy <- S7::new_class("Gassy", parent = structural_term, package = NULL)
  g <- Gassy(label = "gas")
  expect_error(term_build(g, data.frame(x = 1)), "Gassy")
  expect_error(term_build(g, data.frame(x = 1)), "does not implement")
})

test_that("a term class without term_build is told which class is missing it", {
  Bare <- S7::new_class("BareTerm", parent = additive_term, package = NULL)
  b <- Bare(label = "", X = NULL, coef_names = character(0),
            blueprint = list(), penalty = NULL)
  expect_error(term_build(b, data.frame(x = 1)), "BareTerm")
})

test_that("a call that is neither a term nor a covariate names itself", {
  # Recognition by evaluation has a third case besides term and
  # covariate, and it is the one a masked name lands in: mgcv exports
  # s() and te() and segmented exports seg(), so a user with either
  # attached writes our formula and gets theirs. Before this the value
  # travelled to model.matrix and failed there, naming neither the call
  # nor the mask.
  dd <- data.frame(y = rnorm(20), x = runif(20))
  foreign <- function(x) list(kind = "not a term", var = substitute(x))
  e <- new.env(parent = environment())
  assign("foreign", foreign, envir = e)
  f <- y ~ foreign(x)
  environment(f) <- e
  expect_error(interpret_formula(f, dd), "neither a model term nor a covariate")
  expect_error(interpret_formula(f, dd), "foreign", fixed = TRUE)

  # and when the masked name is one of ours, it says so
  masking <- function(x) list(1)
  e2 <- new.env(parent = environment())
  assign("s", masking, envir = e2)
  f2 <- y ~ s(x)
  environment(f2) <- e2
  expect_error(interpret_formula(f2, dd), "modelterms7::s()", fixed = TRUE)

  # a covariate that happens to be a call is untouched
  expect_no_error(interpret_formula(y ~ log(x) + poly(x, 2), dd))
})
