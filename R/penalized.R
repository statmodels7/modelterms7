#' @include term_classes.R generics.R
NULL

#' @title S7 Class for Penalized Parametric Terms
#' @name PenalizedTerm
#'
#' @description
#' A subclass of \code{\link{additive_term}} for a parametric block whose
#' coefficients carry a \pkg{penalties7} penalty. Constructed by
#' \code{\link{ridge}}, \code{\link{lasso}}, \code{\link{scad}} and
#' \code{\link{mcp}}; the four differ only in the penalty their factory
#' attaches at build time, and every derivative, hyperparameter, bound,
#' link and kink is the penalty object's.
#'
#' @inheritParams additive_term
#' @param input The block as given: a one-sided formula or a numeric
#'   matrix.
#' @param input_expr The expression that produced a matrix input, kept so
#'   \code{\link{term_predict}} can re-evaluate it in new data.
#' @param factory The function mapping a coefficient count to the penalty
#'   object.
#'
#' @return An object of class \code{PenalizedTerm}.
#'
#' @seealso \code{\link{ridge}}
#' @examples
#' S7::S7_inherits(ridge(~x), PenalizedTerm)
#' @export
PenalizedTerm <- S7::new_class(
  name = "PenalizedTerm",
  parent = additive_term,
  properties = list(
    input = S7::class_any,
    input_expr = S7::class_any,
    factory = S7::class_function
  )
)

.penalized_spec <- function(x, expr, label, by, factory) {
  if (!is.null(by)) {
    stop("'by' is reserved for a later release and is not implemented yet.",
         call. = FALSE)
  }
  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !nzchar(label)) {
    stop("'label' must be a single non-empty character string.",
         call. = FALSE)
  }
  is_formula <- inherits(x, "formula")
  if (is_formula) {
    if (length(x) != 2L) {
      stop("a formula input must be one-sided, e.g. ~ x1 + x2.",
           call. = FALSE)
    }
  } else {
    x <- as.matrix(x)
    if (!is.numeric(x)) {
      stop("a matrix input must be numeric.", call. = FALSE)
    }
  }
  PenalizedTerm(label = label, input = x, input_expr = expr,
                factory = factory,
                X = NULL, coef_names = character(0),
                blueprint = list(), penalty = NULL)
}

#' Penalized Parametric Terms
#'
#' @description
#' The four classical penalized blocks as model terms: ridge, lasso, SCAD
#' and MCP. Each takes its block as a one-sided formula or as a numeric
#' matrix, and attaches the corresponding \pkg{penalties7} object to the
#' block's coefficients at build time -- \code{\link[penalties7]{ridge_penalty}},
#' \code{\link[penalties7]{lasso_penalty}},
#' \code{\link[penalties7]{elasticnet_penalty}},
#' \code{\link[penalties7]{scad_penalty}},
#' \code{\link[penalties7]{mcp_penalty}} -- so the hyperparameters, their
#' bounds and links, the derivatives and the kink set are the penalty's,
#' never restated by the term.
#'
#' @details
#' A formula input goes through the \code{\link[stats]{model.matrix}}
#' machinery with the intercept removed (a penalized block does not
#' penalize an intercept; the model's intercept lives in the parametric
#' block), and its blueprint records the terms, the factor levels and the
#' contrasts, exactly as \code{\link{linpar}} does. A matrix input is used
#' as given, and its columns are named after the matrix's own column
#' names, or numbered when it has none.
#'
#' Prediction for a matrix input re-evaluates the expression that produced
#' the matrix in the new data, and ONLY there, so the intended use is a
#' matrix column of the model data frame (\code{dd$R <- R}; then
#' \code{ridge(R)} in the formula): a subset of the data then carries the
#' matching rows. A free-standing matrix from the calling environment
#' builds, since its value was captured, but prediction is rejected --
#' resolving it outside the new data would silently reuse the build-time
#' rows.
#'
#' \code{\link{term_smooth}} is \code{TRUE} for \code{ridge} and
#' \code{FALSE} for \code{lasso}, \code{enet}, \code{scad} and
#' \code{mcp}, read from each penalty's kink set.
#'
#' @section The penalties:
#' Writing \eqn{\beta} for the block's coefficients and \eqn{p} for their
#' number, the five attach
#'
#' \deqn{\rho_{\mathrm{ridge}}(\beta) =
#'   \frac{\lVert\beta\rVert_2^2}{2\sigma^2}
#'   + p\log\!\left(\sigma\sqrt{2\pi}\right),}
#'
#' \deqn{\rho_{\mathrm{lasso}}(\beta) = \lambda\lVert\beta\rVert_1
#'   - p\log\!\left(\frac{\lambda}{2}\right),}
#'
#' \deqn{\rho_{\mathrm{enet}}(\beta) = \lambda\left\{
#'   \alpha\lVert\beta\rVert_1
#'   + \frac{1-\alpha}{2}\lVert\beta\rVert_2^2\right\}
#'   + p\log Z(\lambda, \alpha),}
#'
#' and the two non-convex ones, which are defined by their derivative
#' rather than by their value,
#'
#' \deqn{\rho'_{\mathrm{scad}}(t) = \lambda\min\!\left\{1,
#'   \frac{(a\lambda - t)_+}{(a-1)\lambda}\right\},
#'   \qquad
#'   \rho'_{\mathrm{mcp}}(t) = \left(\lambda - \frac{t}{\gamma}\right)_+ ,}
#'
#' for \eqn{t = \lvert\beta_j\rvert \ge 0}, summed over the coefficients.
#' The first three are negative log-priors and keep their normalizing
#' constants, which is what makes their hyperparameters estimable by a
#' marginal criterion; the last two are improper by construction and have
#' none. All five, and the arithmetic behind them, belong to
#' \pkg{penalties7}: the term attaches the object and restates nothing.
#'
#' @param x A one-sided formula or a numeric matrix.
#' @param label A single non-empty string prefixed to the coefficient
#'   names.
#' @param by Reserved for a later release; must be \code{NULL}.
#'
#' @return An object of class \code{\link{PenalizedTerm}} (a
#'   specification; see \code{\link{term_build}}).
#'
#' @examples
#' dd <- data.frame(x1 = rnorm(8), x2 = rnorm(8))
#' built <- term_build(lasso(~ x1 + x2), dd)
#' term_coef_names(built)
#' term_penalty(built)@params
#' term_smooth(built)
#'
#' @seealso \code{\link{linpar}}, \code{\link{s}}, \code{\link{random}}, \code{\link{term_penalty}}, \code{\link{edf}}
#' @export
ridge <- function(x, label = "ridge", by = NULL) {
  .penalized_spec(x, substitute(x), label, by,
                  function(k) penalties7::ridge_penalty(n_coef = k))
}

#' @rdname ridge
#' @export
lasso <- function(x, label = "lasso", by = NULL) {
  .penalized_spec(x, substitute(x), label, by,
                  function(k) penalties7::lasso_penalty(n_coef = k))
}

#' @rdname ridge
#' @export
enet <- function(x, label = "enet", by = NULL) {
  .penalized_spec(x, substitute(x), label, by,
                  function(k) penalties7::elasticnet_penalty(n_coef = k))
}

#' @rdname ridge
#' @export
scad <- function(x, label = "scad", by = NULL) {
  .penalized_spec(x, substitute(x), label, by,
                  function(k) penalties7::scad_penalty(n_coef = k))
}

#' @rdname ridge
#' @export
mcp <- function(x, label = "mcp", by = NULL) {
  .penalized_spec(x, substitute(x), label, by,
                  function(k) penalties7::mcp_penalty(n_coef = k))
}

S7::method(term_build, PenalizedTerm) <- function(term, data, ...) {
  if (inherits(term@input, "formula")) {
    f <- stats::update(term@input, ~ . - 1)
    environment(f) <- environment(term@input)
    mf <- stats::model.frame(f, data,
                             na.action = stats::na.pass,
                             drop.unused.levels = FALSE)
    tt <- attr(mf, "terms")
    X <- stats::model.matrix(tt, mf)
    contr <- attr(X, "contrasts")
    base_names <- colnames(X)
    attr(X, "assign") <- NULL
    attr(X, "contrasts") <- NULL
    term@blueprint <- list(
      kind = "formula",
      terms = stats::delete.response(tt),
      xlev = stats::.getXlevels(tt, mf),
      contrasts = contr
    )
  } else {
    X <- term@input
    if (nrow(X) != nrow(data)) {
      stop(sprintf("the matrix input has %d rows and 'data' has %d.",
                   nrow(X), nrow(data)), call. = FALSE)
    }
    base_names <- colnames(X)
    if (is.null(base_names)) base_names <- as.character(seq_len(ncol(X)))
    term@blueprint <- list(
      kind = "matrix",
      expr = term@input_expr,
      base_names = base_names
    )
  }
  cn <- paste(term@label, base_names, sep = ".")
  colnames(X) <- cn
  rownames(X) <- NULL
  term@X <- X
  term@coef_names <- cn
  term@penalty <- term@factory(ncol(X))
  term
}

S7::method(term_predict, PenalizedTerm) <- function(term, newdata, ...) {
  .assert_built(term)
  bp <- term@blueprint
  if (identical(bp$kind, "formula")) {
    mf <- stats::model.frame(bp$terms, newdata,
                             na.action = stats::na.pass,
                             xlev = bp$xlev)
    X <- stats::model.matrix(bp$terms, mf, contrasts.arg = bp$contrasts)
    attr(X, "assign") <- NULL
    attr(X, "contrasts") <- NULL
  } else {
    X <- tryCatch(as.matrix(eval(bp$expr, newdata, baseenv())),
                  error = function(e) {
                    stop(sprintf(paste("the matrix expression `%s` could not be",
                                       "evaluated in 'newdata' (%s); supply the",
                                       "matrix as a column of the data."),
                                 deparse(bp$expr), conditionMessage(e)),
                         call. = FALSE)
                  })
    if (nrow(X) != nrow(newdata)) {
      stop(sprintf("the matrix expression yields %d rows against %d in 'newdata'.",
                   nrow(X), nrow(newdata)), call. = FALSE)
    }
    if (ncol(X) != length(term@coef_names)) {
      stop(sprintf("the matrix expression yields %d columns against %d at build time.",
                   ncol(X), length(term@coef_names)), call. = FALSE)
    }
  }
  colnames(X) <- term@coef_names
  rownames(X) <- NULL
  X
}
