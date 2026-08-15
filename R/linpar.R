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
#' The block is the model matrix \eqn{X} of the formula and the term's
#' contribution to the predictor is linear in its coefficients,
#'
#' \deqn{\eta = X\beta,}
#'
#' with no penalty attached, so all \eqn{p = \operatorname{ncol}(X)}
#' coefficients are free and \code{\link{edf}} counts every one of them.
#'
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
#' @section Sparse storage:
#' \code{sparse = TRUE} builds the block through
#' \code{\link[Matrix]{sparse.model.matrix}}, which BUILDS it sparse rather
#' than building a dense matrix and compressing it -- the second would cost
#' the memory the choice exists to avoid. Measured at 20000 rows and a factor
#' of 1000 levels, 0.002 s and 1.8 MB against 0.100 s and 161.5 MB, the
#' numbers identical; and a design that would be 32 GB dense builds in 0.02 s
#' and 19 MB, which is what says there is no dense intermediate.
#'
#' It pays where the formula carries a FACTOR OF MANY LEVELS, whose indicator
#' columns hold one non-zero per row. On numeric covariates the block is dense
#' whatever is asked for, and the sparse storage then costs more than it
#' saves; measured on a grouping indicator, the crossover is around a hundred
#' levels. The storage is part of the blueprint, so
#' \code{\link{term_predict}} builds new data the same way.
#'
#' @param formula A one-sided formula, e.g. \code{~ x1 + x2}.
#' @param label A character string; when non-empty it is prefixed to the
#'   coefficient names as \code{label.name}.
#' @param sparse Whether to build the block as a \code{dgCMatrix}. See the
#'   section below.
#' @param contrasts The contrasts for the formula's factors, as a named list
#'   of the kind \code{\link[stats]{model.matrix}}'s \code{contrasts.arg}
#'   takes. \code{NULL}, the default, leaves them to the session's
#'   \code{options("contrasts")}.
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
#' @seealso \code{\link{ridge}}, \code{\link{lasso}}, \code{\link{scad}}, \code{\link{mcp}}, \code{\link{enet}}
#' @export
linpar <- function(formula, label = "", sparse = FALSE, contrasts = NULL) {
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula.", call. = FALSE)
  }
  if (length(formula) != 2L) {
    stop("'formula' must be one-sided, e.g. ~ x1 + x2.", call. = FALSE)
  }
  if (!is.character(label) || length(label) != 1L || is.na(label)) {
    stop("'label' must be a single character string.", call. = FALSE)
  }
  o <- .check_design_opts(sparse, contrasts, "linpar")
  LinparTerm(label = label, formula = formula, sparse = o$sparse,
             contrasts = o$contrasts,
             X = NULL, coef_names = character(0),
             blueprint = list(), penalty = NULL)
}

#' The Model Matrix of a Formula, in Either Storage
#'
#' @description
#' \code{\link[stats]{model.matrix}} or
#' \code{\link[Matrix]{sparse.model.matrix}} on the same terms object, with
#' the bookkeeping stripped either way.
#'
#' @details
#' The sparse route BUILDS the matrix sparse; it does not build a dense one
#' and compress it, which would cost the memory the choice exists to avoid.
#' Measured at 20000 rows and a factor of 1000 levels, 0.002 s and 1.8 MB
#' against \code{stats::model.matrix}'s 0.100 s and 161.5 MB, the numbers
#' identical; and a design that would be 32 GB dense builds in 0.02 s and
#' 19 MB, which is what says there is no dense intermediate.
#'
#' It is worth it where the formula carries a factor of MANY LEVELS, whose
#' indicator columns are one non-zero per row. On a formula of numeric
#' covariates the block is dense whatever is asked for, and the sparse
#' storage then costs more than it saves.
#'
#' @param tt A terms object.
#' @param mf The model frame.
#' @param contrasts The contrasts, or \code{NULL} for the session's.
#' @param sparse Whether to build a \code{dgCMatrix}.
#'
#' @return A numeric matrix or a \code{dgCMatrix}.
#'
#' @seealso \code{\link{linpar}}
#'
#' @keywords internal
.design_matrix <- function(tt, mf, contrasts = NULL, sparse = FALSE) {
  X <- if (isTRUE(sparse)) {
    Matrix::sparse.model.matrix(tt, mf, contrasts.arg = contrasts)
  } else {
    stats::model.matrix(tt, mf, contrasts.arg = contrasts)
  }
  ctr <- attr(X, "contrasts")
  attr(X, "assign") <- NULL
  attr(X, "contrasts") <- NULL
  list(X = X, contrasts = ctr)
}


#' Check a Term's Storage and Contrasts
#'
#' @param sparse What the constructor was given.
#' @param contrasts What the constructor was given.
#' @param what The term's label, for the message.
#'
#' @return A list of the two, validated.
#'
#' @keywords internal
.check_design_opts <- function(sparse, contrasts, what = "this term") {
  if (!is.logical(sparse) || length(sparse) != 1L || is.na(sparse)) {
    stop(sprintf("'sparse' in '%s' must be TRUE or FALSE.", what),
         call. = FALSE)
  }
  if (!is.null(contrasts) && !is.list(contrasts)) {
    stop(sprintf(paste0("'contrasts' in '%s' must be a named list, one entry",
                        " per factor, or NULL."), what), call. = FALSE)
  }
  list(sparse = sparse,
       contrasts = if (is.null(contrasts)) list() else contrasts)
}


S7::method(term_build, LinparTerm) <- function(term, data, ...) {
  mf <- stats::model.frame(term@formula, data,
                           na.action = stats::na.pass,
                           drop.unused.levels = FALSE)
  tt <- attr(mf, "terms")
  # the block is a plain numeric matrix, or a dgCMatrix where the caller
  # asked for one: the model.matrix bookkeeping lives in the blueprint, not
  # on the result
  b <- .design_matrix(tt, mf,
                      if (length(term@contrasts)) term@contrasts else NULL,
                      term@sparse)
  X <- b$X
  cn <- colnames(X)
  if (nzchar(term@label)) cn <- paste(term@label, cn, sep = ".")
  colnames(X) <- cn
  term@X <- X
  term@coef_names <- cn
  term@blueprint <- list(
    terms = stats::delete.response(tt),
    xlev = stats::.getXlevels(tt, mf),
    contrasts = b$contrasts,
    sparse = term@sparse
  )
  term
}

S7::method(term_predict, LinparTerm) <- function(term, newdata, ...) {
  .assert_built(term)
  bp <- term@blueprint
  mf <- stats::model.frame(bp$terms, newdata,
                           na.action = stats::na.pass,
                           xlev = bp$xlev)
  # the STORAGE is part of the blueprint: a prediction that densified would
  # spend at new data what the build was careful not to
  X <- .design_matrix(bp$terms, mf, bp$contrasts, isTRUE(bp$sparse))$X
  colnames(X) <- term@coef_names
  X
}
