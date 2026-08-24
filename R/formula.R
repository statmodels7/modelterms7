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
#' Splits a model formula into a response and a named list of term
#' specifications. A call on the right-hand side is evaluated, and if its value
#' inherits from [model_term()] it becomes a term of its own; everything else,
#' meaning bare covariates, transformations such as `log(x)` and interactions,
#' is collected into a single [linpar()] block with the usual
#' [stats::model.matrix()] conventions. The terms come back unbuilt, so the
#' caller decides which data each is built against.
#'
#' @details
#' # Recognition by evaluation
#'
#' A term constructor is identified by what its call returns, so a term class
#' defined outside the package works in a formula the day it is written, with
#' no list of special names to amend. `log(x)` evaluates to a numeric vector
#' and stays a covariate; `s(x, k = 5)` evaluates to a `SmoothTerm` and is
#' routed as a term.
#'
#' Some labels are never evaluated: `:`, `*`, `^`, `%in%`, `+`, `-`, `(` and
#' `I`. A bare interaction `x1:x2` parses as a call to `:`, which on numeric
#' vectors is the sequence operator and would return something unrelated to the
#' column [stats::model.matrix()] builds. Those labels go straight to the
#' parametric block.
#'
#' A call whose value is neither a term nor something a model matrix can hold
#' throws, naming the label and the class it produced. The commonest way to
#' arrive there is a masked name. `mgcv` also exports `s()` and `te()` and
#' `segmented` exports `seg()`, so where the package exports a term of that
#' name and the name resolves elsewhere, the message says where it was found
#' and suggests `modelterms7::`.
#'
#' # The response and the intercept
#'
#' The left-hand side, when there is one, is evaluated in `data`: a plain
#' expression gives a numeric vector, and a response constructor such as
#' [cens()] gives its object. A one-sided formula gives `response = NULL`.
#'
#' The intercept convention is the formula's own, carried into the collected
#' parametric block. `y ~ ridge(~ g)` still produces an intercept-only
#' `linpar` block, and `y ~ ridge(~ g) - 1` produces no `linpar` block at all.
#'
#' # One covariate is removed
#'
#' A [seg()] or [jseg()] term built with `linear = TRUE`, which is the default,
#' contributes the same column the bare covariate would, so `y ~ x + seg(x)` is
#' rank deficient by one. The term owns that effect, so the covariate is
#' dropped from the parametric block and a warning names both. Writing
#' `seg(x, linear = FALSE)` keeps the linear effect outside the term.
#'
#' Only the bare main effect is removed, matched on the deparsed expression: an
#' interaction spans no main effect and is left alone. A term that spans the
#' same direction without being that one column is reported by a warning and
#' left unmodified. A spline basis is the case: it contains the line, and its
#' penalty leaves the line unpenalized, so `s(x) + seg(x)` is confounded too.
#' There is no single column to remove there, and reshaping another term is
#' not this function's business.
#'
#' # What it does not do
#'
#' Nothing is built: every element of `terms` is a specification, and
#' [term_is_built()] is `FALSE` for all of them. `data` is read for
#' [stats::terms()]'s variable classification and for evaluating the term calls
#' and the response, and the term constructors themselves do not touch it.
#'
#' Labels are deduplicated by [stats::terms()] before any of this, so
#' `y ~ ridge(~ x) + ridge(~ x)` gives one term, as `y ~ x + x` gives one
#' column.
#'
#' @param formula A two-sided or one-sided model formula. Anything else throws
#'   `"'formula' must be a formula."`. Its environment is where a term call's
#'   symbols are looked up when `data` does not carry them, and it is carried
#'   onto the formula of the implicit `linpar` block.
#' @param data A data frame in which the formula's symbols are evaluated.
#'   Anything else throws `"'data' must be a data frame."`.
#' @param linpar A named list of arguments for the **implicit** [linpar()]
#'   term, the one the bare covariates collapse into: `sparse` and `contrasts`.
#'   This is the only place they can be given, that term never being written by
#'   the caller. Empty by default, which leaves [linpar()]'s own defaults. A
#'   `linpar()` term the caller writes out takes its own arguments and ignores
#'   this one. Anything that is not a list throws.
#'
#' @return A list of four elements:
#'   \describe{
#'     \item{`response`}{The evaluated left-hand side: a numeric vector, a
#'       [censored_response()], or `NULL` for a one-sided formula.}
#'     \item{`terms`}{A named list of unbuilt term specifications. The
#'       collected parametric block comes first under the name `"linpar"` and
#'       is absent when the formula has no bare covariates and no intercept.
#'       Every other name is the term's label as it appears in the formula,
#'       deparsed, such as `"s(x2, k = 5)"`.}
#'     \item{`intercept`}{`TRUE` unless the formula removes the intercept.}
#'     \item{`formula`}{The input, unchanged.}
#'   }
#'
#' @seealso [linpar()] for the block the covariates collapse into, [cens()] for
#'   the response constructor, [check_term()] for validating one of the
#'   returned specifications, and [term_build()] for building it.
#'
#' @examples
#' dd <- data.frame(y = rnorm(20), x1 = 1:20, x2 = runif(20),
#'                  g = factor(rep(letters[1:4], 5)))
#'
#' # Covariates and transformations collapse into one parametric block.
#' out <- interpret_formula(y ~ x1 + log(x2), dd)
#' names(out$terms)
#' out$terms$linpar@formula
#'
#' # A constructor call becomes a term, keyed by its label in the formula.
#' out2 <- interpret_formula(y ~ x1 + s(x2, k = 5) + ridge(~ g), dd)
#' names(out2$terms)
#' vapply(out2$terms, function(t) class(t)[1], character(1))
#'
#' # Nothing is built yet.
#' vapply(out2$terms, term_is_built, logical(1))
#'
#' # The intercept convention is the formula's own.
#' names(interpret_formula(y ~ ridge(~ g), dd)$terms)
#' names(interpret_formula(y ~ ridge(~ g) - 1, dd)$terms)
#'
#' # An interaction is a covariate: `:` is never evaluated.
#' interpret_formula(y ~ x1:x2 + g, dd)$terms$linpar@formula
#'
#' # seg() carries the linear effect, so the bare covariate is removed.
#' w <- interpret_formula(y ~ x1 + seg(x1), dd)
#' w$terms$linpar@formula          # x1 is gone; the term carries it
#'
#' # Unless the term is told not to own it.
#' interpret_formula(y ~ x1 + seg(x1, linear = FALSE), dd)$terms$linpar@formula
#'
#' # Arguments for the implicit block go through `linpar`.
#' sp <- interpret_formula(y ~ g, dd, linpar = list(sparse = TRUE))
#' class(term_matrix(term_build(sp$terms$linpar, dd)))
#'
#' # A call that returns neither a term nor a covariate is refused by name.
#' e <- new.env(parent = globalenv())
#' assign("s", function(x, ...) structure(list(), class = "gamObject"), envir = e)
#' f <- y ~ s(x1)
#' environment(f) <- e
#' try(interpret_formula(f, dd))
#'
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

