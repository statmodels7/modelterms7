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
#' saves. \code{sparse = NULL}, the default, settles it at build from the
#' design: the dense indicator part holds \code{n} times its column count in
#' cells against one non-zero per row, and the two routes cross at about
#' \eqn{10^5} of those cells, which is the rule \code{\link{.resolve_sparse}}
#' applies. \code{TRUE} and \code{FALSE} override it. The storage that was
#' settled is part of the blueprint, so \code{\link{term_predict}} builds new
#' data the same way.
#'
#' @param formula A one-sided formula, e.g. \code{~ x1 + x2}.
#' @param label A character string; when non-empty it is prefixed to the
#'   coefficient names as \code{label.name}.
#' @param sparse Whether to build the block as a \code{dgCMatrix}.
#'   \code{NULL}, the default, settles it from the design. See the section
#'   below.
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
linpar <- function(formula, label = "", sparse = NULL, contrasts = NULL) {
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
#' @param sparse What the constructor was given: \code{TRUE}, \code{FALSE}, or
#'   \code{NULL} for the storage to be settled at build from the design.
#' @param contrasts What the constructor was given.
#' @param what The term's label, for the message.
#'
#' @return A list of the two, validated.
#'
#' @keywords internal
.check_design_opts <- function(sparse, contrasts, what = "this term") {
  if (!is.null(sparse) &&
      (!is.logical(sparse) || length(sparse) != 1L || is.na(sparse))) {
    stop(sprintf(paste0("'sparse' in '%s' must be TRUE, FALSE, or NULL to",
                        " settle it from the design."), what), call. = FALSE)
  }
  if (!is.null(contrasts) && !is.list(contrasts)) {
    stop(sprintf(paste0("'contrasts' in '%s' must be a named list, one entry",
                        " per factor, or NULL."), what), call. = FALSE)
  }
  list(sparse = sparse,
       contrasts = if (is.null(contrasts)) list() else contrasts)
}


#' How Many Columns of a Model Matrix Come From Factors
#'
#' @description
#' The number of columns a terms object contributes through indicator
#' variables, counted from the model frame without building the matrix.
#'
#' @details
#' A term contributes the product of its variables' level counts, a numeric
#' variable counting one; a term carrying no factor contributes columns that
#' are dense whatever the storage, and is not counted. The count is an upper
#' bound, contrasts dropping one level per factor, which is the right side to
#' err on: it is read against a threshold below which the sparse route loses
#' little and above which it wins by orders.
#'
#' @param tt A terms object.
#' @param mf The model frame it was built from.
#'
#' @return A single number, zero when no term carries a factor.
#'
#' @seealso \code{\link{.resolve_sparse}}
#'
#' @keywords internal
.indicator_cols <- function(tt, mf) {
  fac <- attr(tt, "factors")
  if (is.null(fac) || !length(fac) || !ncol(fac)) return(0)
  nlev <- vapply(rownames(fac), function(v) {
    x <- mf[[v]]
    if (is.null(x)) 1
    else if (is.factor(x)) as.numeric(nlevels(x))
    else if (is.character(x) || is.logical(x)) as.numeric(length(unique(x)))
    else 1
  }, numeric(1))
  tot <- 0
  for (j in seq_len(ncol(fac))) {
    inn <- fac[, j] > 0
    if (!any(inn) || !any(nlev[inn] > 1)) next
    tot <- tot + prod(nlev[inn])
  }
  tot
}


#' Settle Whether a Block Is Built Sparse
#'
#' @description
#' Passes an explicit \code{TRUE} or \code{FALSE} through, and where the
#' caller left \code{NULL} decides from the size of the block.
#'
#' @details
#' The dense indicator part holds \code{n * ncol_ind} cells where the sparse
#' one holds one non-zero per row, so that product is what the two routes are
#' separated by, and the threshold is read off it rather than off a count of
#' levels. Measured end to end on \code{y ~ 0 + g + s(x)} over eighteen
#' combinations of sample size and level count, the routes cross at about
#' \eqn{10^5} cells: at \eqn{n = 1000} the sparse route loses at every level
#' count up to sixty (\eqn{6 \times 10^4} cells, 0.93 times the dense route),
#' at \eqn{n = 5000} it crosses between fifteen and twenty-five levels, and at
#' \eqn{n = 20000} between six and ten. The same threshold accounts for the
#' large cases: four hundred levels at \eqn{n = 20000} are \eqn{8 \times 10^6}
#' cells and run 43.75 times faster, with the log-likelihood identical.
#'
#' A design carrying no factor has no indicator part, so the product is zero
#' and the block is built dense, which is what the measurements ask for there
#' (0.66 to 0.90 times the dense route on purely continuous covariates).
#'
#' @param sparse \code{TRUE}, \code{FALSE}, or \code{NULL} to decide here.
#' @param n The number of rows.
#' @param ncol_ind The columns coming from indicators, from
#'   \code{\link{.indicator_cols}}.
#'
#' @return A single logical.
#'
#' @seealso \code{\link{.indicator_cols}}, \code{\link{linpar}}
#'
#' @keywords internal
.resolve_sparse <- function(sparse, n, ncol_ind) {
  if (!is.null(sparse)) return(isTRUE(sparse))
  isTRUE(is.finite(n) && is.finite(ncol_ind) && n * ncol_ind > 1e5)
}


S7::method(term_build, LinparTerm) <- function(term, data, ...) {
  mf <- stats::model.frame(term@formula, data,
                           na.action = stats::na.pass,
                           drop.unused.levels = FALSE)
  tt <- attr(mf, "terms")
  # the block is a plain numeric matrix, or a dgCMatrix where the caller
  # asked for one or where the design is large enough that one pays: the
  # model.matrix bookkeeping lives in the blueprint, not on the result
  sp <- .resolve_sparse(term@sparse, nrow(mf), .indicator_cols(tt, mf))
  b <- .design_matrix(tt, mf,
                      if (length(term@contrasts)) term@contrasts else NULL,
                      sp)
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
    # the SETTLED storage, not what the constructor was given: new data may
    # carry any number of rows, and a prediction that decided again could
    # build a block of a different kind from the one that was fitted
    sparse = sp
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
