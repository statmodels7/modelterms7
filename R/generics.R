#' @include term_classes.R
NULL

#' @title Build a Term on Data
#'
#' @description
#' Turns a term specification into a built term: the design block is
#' computed from the data, the coefficient names are assigned, and the
#' blueprint that reproduces the mapping on new data is recorded. The
#' returned object is a copy of the specification with those properties
#' filled; the specification itself is unchanged.
#'
#' @param term An object inheriting from class \code{\link{model_term}}.
#' @param data A data frame.
#' @param ... Passed to methods.
#'
#' @return A built term of the same class as \code{term}.
#'
#' @examples
#' built <- term_build(linpar(~x), data.frame(x = 1:4))
#' term_matrix(built)
#'
#' @export
term_build <- S7::new_generic("term_build", "term",
  function(term, data, ...) {
    if (!is.data.frame(data)) {
      stop("'data' must be a data frame.", call. = FALSE)
    }
    S7::S7_dispatch()
  })

S7::method(term_build, model_term) <- function(term, data, ...) {
  stop(sprintf("the term class '%s' does not implement term_build().",
               attr(S7::S7_class(term), "name")), call. = FALSE)
}

S7::method(term_build, structural_term) <- function(term, data, ...) {
  stop("structural terms are reserved for a later release; none is implemented yet.",
       call. = FALSE)
}

#' @title Whether a Term Has Been Built
#'
#' @description
#' \code{TRUE} for a term returned by \code{\link{term_build}} and
#' \code{FALSE} for a bare specification. The accessors
#' \code{\link{term_matrix}}, \code{\link{term_npar}},
#' \code{\link{term_coef_names}} and \code{\link{term_predict}} reject a
#' specification, and this predicate is the test they use.
#'
#' @param term An object inheriting from class \code{\link{model_term}}.
#'
#' @return A logical scalar.
#'
#' @examples
#' term_is_built(linpar(~x))
#' term_is_built(term_build(linpar(~x), data.frame(x = 1:4)))
#'
#' @export
term_is_built <- function(term) {
  if (!S7::S7_inherits(term, model_term)) {
    stop("'term' must inherit from 'model_term'.", call. = FALSE)
  }
  S7::S7_inherits(term, additive_term) && length(term@coef_names) > 0L
}

.assert_built <- function(term) {
  if (!term_is_built(term)) {
    stop("the term has not been built; call term_build(term, data) first.",
         call. = FALSE)
  }
  invisible(term)
}

#' @title Design Block of a Built Term
#'
#' @description
#' The \eqn{n \times k} design block of a built additive term, with the
#' term's coefficient names as column names.
#'
#' @param term A built term (see \code{\link{term_build}}).
#' @param ... Passed to methods.
#'
#' @return A numeric matrix.
#'
#' @examples
#' term_matrix(term_build(linpar(~x), data.frame(x = 1:4)))
#'
#' @export
term_matrix <- S7::new_generic("term_matrix", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_matrix, additive_term) <- function(term, ...) {
  .assert_built(term)
  term@X
}

#' @title Penalty of a Term
#'
#' @description
#' The penalty attached to the term's coefficients, or \code{NULL} for an
#' unpenalized term. The hyperparameters, their bounds and links, and every
#' derivative in the coefficients and the hyperparameters are the penalty
#' object's, not the term's.
#'
#' @param term An object inheriting from class \code{\link{additive_term}}.
#' @param ... Passed to methods.
#'
#' @return A penalty object, or \code{NULL}.
#'
#' @examples
#' term_penalty(linpar(~x))
#'
#' @export
term_penalty <- S7::new_generic("term_penalty", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_penalty, additive_term) <- function(term, ...) {
  term@penalty
}

#' @title Number of Coefficients of a Built Term
#'
#' @description The number of columns of the term's design block.
#'
#' @param term A built term (see \code{\link{term_build}}).
#' @param ... Passed to methods.
#'
#' @return An integer.
#'
#' @examples
#' term_npar(term_build(linpar(~x), data.frame(x = 1:4)))
#'
#' @export
term_npar <- S7::new_generic("term_npar", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_npar, additive_term) <- function(term, ...) {
  .assert_built(term)
  ncol(term@X)
}

#' @title Coefficient Names of a Built Term
#'
#' @description
#' The names of the term's coefficients, prefixed by the term's label when
#' the label is non-empty.
#'
#' @param term A built term (see \code{\link{term_build}}).
#' @param ... Passed to methods.
#'
#' @return A character vector.
#'
#' @examples
#' term_coef_names(term_build(linpar(~x), data.frame(x = 1:4)))
#'
#' @export
term_coef_names <- S7::new_generic("term_coef_names", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_coef_names, additive_term) <- function(term, ...) {
  .assert_built(term)
  term@coef_names
}

#' @title Whether a Term's Penalized Objective Is Smooth
#'
#' @description
#' \code{TRUE} when the term's contribution to the penalized objective is
#' differentiable in the coefficients. The answer is read from the penalty
#' rather than declared by the term: an unpenalized term is smooth, and a
#' penalized one is smooth exactly when its penalty declares no kinks, so a
#' term cannot disagree with its own penalty. The model layer uses this
#' flag to split the coefficient vector into the block the classical
#' optimizers handle and the block that needs non-smooth strategies.
#'
#' @param term An object inheriting from class \code{\link{additive_term}}.
#' @param ... Passed to methods.
#'
#' @return A logical scalar.
#'
#' @examples
#' term_smooth(linpar(~x))
#'
#' @export
term_smooth <- S7::new_generic("term_smooth", "term",
  function(term, ...) S7::S7_dispatch())

# a hyperparameter value inside each domain, at which the kink set is asked
# for; the kinks are structural, so any admissible value answers the question
.penalty_probe_theta <- function(pen) {
  stats::setNames(lapply(pen@params, function(p) {
    b <- pen@params_bounds[[p]]
    if (is.finite(b[1L]) && is.finite(b[2L])) return(mean(b))
    if (is.finite(b[1L])) return(b[1L] + 1)
    if (is.finite(b[2L])) return(b[2L] - 1)
    0
  }), pen@params)
}

S7::method(term_smooth, additive_term) <- function(term, ...) {
  pen <- term@penalty
  if (is.null(pen)) return(TRUE)
  kinks <- penalties7::penalty_kinks(pen, .penalty_probe_theta(pen))
  length(kinks) == 0L
}

#' @title Design Block on New Data
#'
#' @description
#' Applies a built term's mapping to new data, reproducing the block the
#' term would have produced had the new rows been part of the original
#' data: factor levels, contrasts and any constants recorded in the
#' blueprint at build time are reused, never recomputed. New data carrying
#' a factor level unknown to the blueprint is rejected.
#'
#' @param term A built term (see \code{\link{term_build}}).
#' @param newdata A data frame.
#' @param ... Passed to methods.
#'
#' @return A numeric matrix with \code{nrow(newdata)} rows and one column
#'   per coefficient.
#'
#' @examples
#' built <- term_build(linpar(~x), data.frame(x = 1:4))
#' term_predict(built, data.frame(x = c(0.5, 2.5)))
#'
#' @export
term_predict <- S7::new_generic("term_predict", "term",
  function(term, newdata, ...) {
    if (!is.data.frame(newdata)) {
      stop("'newdata' must be a data frame.", call. = FALSE)
    }
    S7::S7_dispatch()
  })

# --- printing ---------------------------------------------------------------

S7::method(print, model_term) <- function(x, ...) {
  cls <- attr(S7::S7_class(x), "name")
  lab <- if (nzchar(x@label)) sprintf(" '%s'", x@label) else ""
  if (term_is_built(x)) {
    cat(sprintf("<%s>%s built: %d coefficient%s\n", cls, lab,
                ncol(x@X), if (ncol(x@X) == 1L) "" else "s"))
  } else {
    cat(sprintf("<%s>%s (specification; call term_build() with data)\n",
                cls, lab))
  }
  invisible(x)
}
