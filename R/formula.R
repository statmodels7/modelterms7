#' @include term_classes.R generics.R linpar.R
NULL

# operators whose term labels must never be evaluated directly: a bare
# interaction such as x1:x2 parses as a call to `:`, which on numeric
# vectors is the sequence operator and returns something unrelated to the
# model matrix column.
.formula_operators <- c(":", "*", "^", "%in%", "+", "-", "(", "I")

#' Interpret a Model Formula Into Terms
#'
#' @description
#' Splits a model formula into a response specification and a list of term
#' specifications. Term constructors are recognized by what they return: a
#' call on the right-hand side whose value inherits from
#' \code{\link{model_term}} becomes a term of its own, and everything else
#' (bare covariates, transformations such as \code{log(x)}, interactions)
#' is collected into one \code{\link{linpar}} block with the usual
#' \code{\link[stats]{model.matrix}} conventions.
#'
#' @details
#' Recognition by evaluation is what makes the interpreter extensible: a
#' term class defined outside the package works in a formula the day it is
#' written, with no list of special names to amend. \code{log(x)}
#' evaluates to a numeric vector and stays a covariate; a constructor call
#' evaluates to a term specification and is routed as one. Interaction
#' labels and bare symbols are never evaluated directly.
#'
#' The left-hand side, when present, is evaluated in the data: a plain
#' expression gives a numeric response, and a response constructor such as
#' \code{\link{cens}} gives its response object. The intercept convention
#' is the formula's own, carried into the collected parametric block, so
#' \code{y ~ ridge_like(R)} still produces an intercept-only
#' \code{linpar} block and \code{y ~ ridge_like(R) - 1} produces none.
#'
#' @param formula A model formula.
#' @param data A data frame in which the formula's symbols are evaluated.
#'
#' @return A list with elements \code{response} (the evaluated left-hand
#'   side, or \code{NULL} for a one-sided formula), \code{terms} (a named
#'   list of term specifications, the collected parametric block first
#'   under the name \code{"linpar"}), \code{intercept} (logical) and
#'   \code{formula} (the input).
#'
#' @examples
#' dd <- data.frame(y = rnorm(6), x1 = 1:6, x2 = runif(6))
#' out <- interpret_formula(y ~ x1 + log(x2), dd)
#' names(out$terms)
#'
#' @export
interpret_formula <- function(formula, data) {
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  env <- environment(formula)
  if (is.null(env)) env <- baseenv()

  response <- NULL
  if (length(formula) == 3L) {
    response <- eval(formula[[2L]], data, env)
  }

  tt <- stats::terms(formula, data = data)
  labels <- attr(tt, "term.labels")
  intercept <- attr(tt, "intercept") == 1L

  specials <- list()
  ordinary <- character(0)
  for (lb in labels) {
    ex <- str2lang(lb)
    res <- NULL
    if (is.call(ex) &&
        !as.character(ex[[1L]])[1L] %in% .formula_operators) {
      res <- tryCatch(eval(ex, data, env), error = function(e) NULL)
    }
    if (S7::S7_inherits(res, model_term)) {
      specials[[lb]] <- res
    } else {
      ordinary <- c(ordinary, lb)
    }
  }

  terms_list <- list()
  if (length(ordinary)) {
    f <- stats::reformulate(ordinary, intercept = intercept)
    environment(f) <- env
    terms_list$linpar <- linpar(f)
  } else if (intercept) {
    f <- ~1
    environment(f) <- env
    terms_list$linpar <- linpar(f)
  }
  terms_list <- c(terms_list, specials)

  list(response = response, terms = terms_list,
       intercept = intercept, formula = formula)
}
