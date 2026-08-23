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
#' [model_term()] becomes a term of its own, and everything else
#' (bare covariates, transformations such as `log(x)`, interactions)
#' is collected into one [linpar()] block with the usual
#' [stats::model.matrix()] conventions.
#'
#' @details
#' Recognition by evaluation is what makes the interpreter extensible: a
#' term class defined outside the package works in a formula the day it is
#' written, with no list of special names to amend. `log(x)`
#' evaluates to a numeric vector and stays a covariate; a constructor call
#' evaluates to a term specification and is routed as one. Interaction
#' labels and bare symbols are never evaluated directly.
#'
#' The left-hand side, when present, is evaluated in the data: a plain
#' expression gives a numeric response, and a response constructor such as
#' [cens()] gives its response object. The intercept convention
#' is the formula's own, carried into the collected parametric block, so
#' `y ~ ridge_like(R)` still produces an intercept-only
#' `linpar` block and `y ~ ridge_like(R) - 1` produces none.
#'
#' One covariate is removed rather than collected. A [seg()] or
#' [jseg()] term carrying the linear effect contributes the same
#' column the bare covariate would, so `y ~ x + seg(x)` is rank
#' deficient by one; the term owns that effect, which is what
#' `linear = TRUE` says, so the covariate is dropped from the
#' parametric block and the removal is reported with a warning.
#' `seg(x, linear = FALSE)` keeps the linear effect outside the term
#' instead. Only the bare main effect is removed: an interaction spans no
#' main effect and is left alone, and another term spanning the same
#' direction, as a spline basis does, is reported without being modified.
#'
#' @param formula A model formula.
#' @param data A data frame in which the formula's symbols are evaluated.
#' @param linpar Arguments for the IMPLICIT [linpar()] term, the one
#'   the bare covariates collapse into, as a named list -- `sparse` and
#'   `contrasts`. It is the only place they can be given: that term is
#'   never written by the caller. Empty, the default, leaves the constructor's
#'   own.
#'
#' @return A list with elements `response` (the evaluated left-hand
#'   side, or `NULL` for a one-sided formula), `terms` (a named
#'   list of term specifications, the collected parametric block first
#'   under the name `"linpar"`), `intercept` (logical) and
#'   `formula` (the input).
#'
#' @examples
#' dd <- data.frame(y = rnorm(6), x1 = 1:6, x2 = runif(6))
#' out <- interpret_formula(y ~ x1 + log(x2), dd)
#' names(out$terms)
#'
#' @seealso [cens()], [check_term()]
#' @export
interpret_formula <- function(formula, data, linpar = list()) {
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  if (!is.list(linpar)) {
    stop("'linpar' must be a list of arguments for linpar().", call. = FALSE)
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
    ok <- FALSE
    if (is.call(ex) &&
        !as.character(ex[[1L]])[1L] %in% .formula_operators) {
      res <- tryCatch({ v <- eval(ex, data, env); ok <- TRUE; v },
                      error = function(e) NULL)
    }
    if (S7::S7_inherits(res, model_term)) {
      specials[[lb]] <- res
    } else {
      if (ok) .reject_unusable(res, ex, lb, env)
      ordinary <- c(ordinary, lb)
    }
  }

  ordinary <- .absorb_linear(ordinary, specials)

  # The IMPLICIT linpar is the one a caller never writes -- the bare
  # covariates of the formula collapsed into one term -- so the only place
  # its arguments can come from is here.
  mk <- function(f) do.call(modelterms7::linpar, c(list(f), linpar))
  terms_list <- list()
  if (length(ordinary)) {
    f <- stats::reformulate(ordinary, intercept = intercept)
    environment(f) <- env
    terms_list$linpar <- mk(f)
  } else if (intercept) {
    f <- ~1
    environment(f) <- env
    terms_list$linpar <- mk(f)
  }
  terms_list <- c(terms_list, specials)

  list(response = response, terms = terms_list,
       intercept = intercept, formula = formula)
}

# A break-point term that carries the linear effect carries the column the
# bare covariate would contribute, so the two are EXACTLY collinear and the
# equation is rank deficient by one. The term owns it -- that is what
# `linear = TRUE` says -- so the covariate is removed from the parametric
# part and the removal is reported, `y ~ x + seg(x)` being an ordinary
# thing to write and a silent singularity a poor answer to it. Writing
# `seg(x, linear = FALSE)` is the way to keep the linear effect outside.
#
# Only the bare main effect is removed, matched by the deparsed expression:
# an interaction spans no main effect and is left alone. Another term that
# spans the same direction is reported and NOT modified -- a spline basis
# contains the constant and the line, and its penalty leaves the line
# unpenalized, so `s(x) + seg(x)` is confounded too -- because there is no
# one column to remove there and reshaping another term is not this
# function's business.
.absorb_linear <- function(ordinary, specials) {
  keep <- vapply(specials, function(tm) {
    S7::S7_inherits(tm, SegTerm) && isTRUE(tm@linear)
  }, logical(1))
  if (!any(keep)) return(ordinary)
  segs <- specials[keep]
  owned <- stats::setNames(
    vapply(segs, function(tm) tm@kind, character(1)),
    vapply(segs, function(tm) paste(deparse(tm@var), collapse = ""),
           character(1)))
  for (lb in intersect(ordinary, names(owned))) {
    warning(sprintf(paste("the covariate '%s' is exactly collinear with the",
                          "linear effect that '%s' carries, and has been",
                          "removed from the parametric part. Write %s(%s,",
                          "linear = FALSE) to keep the linear effect outside",
                          "the term instead."),
                    lb, owned[[lb]], owned[[lb]], lb), call. = FALSE)
  }
  # another term spanning the same direction is reported and left alone
  for (nm in names(specials)) {
    tm <- specials[[nm]]
    if (!("vars" %in% S7::prop_names(tm))) next
    v <- intersect(vapply(tm@vars, function(e)
      paste(deparse(e), collapse = ""), character(1)), names(owned))
    if (length(v)) {
      warning(sprintf(paste("'%s' spans the linear effect in '%s' that '%s'",
                            "also carries, so the two are confounded along",
                            "it. Neither is modified; write %s(%s,",
                            "linear = FALSE) if that is what was meant."),
                      nm, v[1L], owned[[v[1L]]], owned[[v[1L]]], v[1L]),
              call. = FALSE)
    }
  }
  setdiff(ordinary, names(owned))
}

# A call in a formula is a term when its value inherits model_term and a
# covariate otherwise, which leaves a third case: a value that is neither.
# It reaches model.matrix and fails there, several frames from the cause
# and without naming it. The commonest way to arrive is a masked name --
# s() and te() are also exported by mgcv, seg() by segmented -- so the
# message says which package supplied the function that was called.
.reject_unusable <- function(res, ex, lb, env) {
  usable <- is.numeric(res) || is.logical(res) || is.character(res) ||
    is.factor(res) || is.matrix(res) || inherits(res, "Date") ||
    inherits(res, "difftime")
  if (usable) return(invisible(NULL))
  fn <- as.character(ex[[1L]])[1L]
  ours <- tryCatch(get(fn, envir = asNamespace("modelterms7"),
                       mode = "function"),
                   error = function(e) NULL)
  theirs <- tryCatch(get(fn, envir = env, mode = "function"),
                     error = function(e) NULL)
  extra <- ""
  if (!is.null(ours) && !is.null(theirs) && !identical(ours, theirs)) {
    where <- environmentName(environment(theirs))
    extra <- sprintf(
      paste0(" The name '%s' is masked here%s, and modelterms7 exports a",
             " term of that name: write modelterms7::%s()."),
      fn, if (nzchar(where)) sprintf(" by '%s'", where) else "", fn)
  }
  stop(sprintf(
    paste0("the term '%s' evaluated to an object of class %s, which is",
           " neither a model term nor a covariate.%s"),
    lb, paste(sQuote(class(res), FALSE), collapse = "/"), extra),
    call. = FALSE)
}

