#' @include term_classes.R generics.R
NULL

#' Unpenalized Parametric Term
#'
#' @description
#' Creates the specification of an unpenalized parametric block: the design
#' matrix of a one-sided formula, with the usual \code{\link[stats]{model.matrix}}
#' conventions for factors, contrasts, interactions and the intercept.
#'
#' @details
#' \code{\link{interpret_formula}} collects the bare covariates of a model
#' formula into one term of this kind, so \code{y ~ x1 + x2} and
#' \code{y ~ linpar(~ x1 + x2)} produce the same block; the explicit
#' constructor exists for callers who want several parametric blocks with
#' distinct labels.
#'
#' Building the term records a blueprint: the terms object, the factor
#' levels and the contrasts. \code{\link{term_predict}} reapplies the
#' mapping through that blueprint, so a factor column in new data is
#' encoded against the levels seen at build time, and a level the
#' blueprint does not know is rejected rather than re-encoded. Missing
#' values are propagated (\code{na.pass}), never dropped, so the block
#' stays row-aligned with the response.
#'
#' @param formula A one-sided formula, e.g. \code{~ x1 + x2}.
#' @param label A character string; when non-empty it is prefixed to the
#'   coefficient names as \code{label.name}.
#'
#' @return An object of class \code{\link{LinparTerm}} (a specification;
#'   see \code{\link{term_build}}).
#'
#' @examples
#' dd <- data.frame(x = 1:4, g = factor(c("a", "a", "b", "b")))
#' built <- term_build(linpar(~ x + g), dd)
#' term_matrix(built)
#' term_coef_names(built)
#'
#' @export
linpar <- function(formula, label = "") {
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula.", call. = FALSE)
  }
  if (length(formula) != 2L) {
    stop("'formula' must be one-sided, e.g. ~ x1 + x2.", call. = FALSE)
  }
  if (!is.character(label) || length(label) != 1L || is.na(label)) {
    stop("'label' must be a single character string.", call. = FALSE)
  }
  LinparTerm(label = label, formula = formula,
             X = NULL, coef_names = character(0),
             blueprint = list(), penalty = NULL)
}

S7::method(term_build, LinparTerm) <- function(term, data, ...) {
  mf <- stats::model.frame(term@formula, data,
                           na.action = stats::na.pass,
                           drop.unused.levels = FALSE)
  tt <- attr(mf, "terms")
  X <- stats::model.matrix(tt, mf)
  contr <- attr(X, "contrasts")
  cn <- colnames(X)
  if (nzchar(term@label)) cn <- paste(term@label, cn, sep = ".")
  # the block is a plain numeric matrix: the model.matrix bookkeeping lives
  # in the blueprint, not on the result
  attr(X, "assign") <- NULL
  attr(X, "contrasts") <- NULL
  colnames(X) <- cn
  term@X <- X
  term@coef_names <- cn
  term@blueprint <- list(
    terms = stats::delete.response(tt),
    xlev = stats::.getXlevels(tt, mf),
    contrasts = contr
  )
  term
}

S7::method(term_predict, LinparTerm) <- function(term, newdata, ...) {
  .assert_built(term)
  bp <- term@blueprint
  mf <- stats::model.frame(bp$terms, newdata,
                           na.action = stats::na.pass,
                           xlev = bp$xlev)
  X <- stats::model.matrix(bp$terms, mf, contrasts.arg = bp$contrasts)
  attr(X, "assign") <- NULL
  attr(X, "contrasts") <- NULL
  colnames(X) <- term@coef_names
  X
}
