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
#'   object. It is called with the diagonal map as a second argument where
#'   \code{standardize} asks for one, so a factory that will never be
#'   standardized may take the count alone.
#' @param standardize Whether the block's columns are put on a common scale
#'   by the penalty's diagonal map.
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
    factory = S7::class_function,
    standardize = S7::new_property(S7::class_logical, default = FALSE)
  )
)

.penalized_spec <- function(x, expr, label, by, standardize, factory,
                            hyper = list(), extra = list(), grid = list(),
                            min_ratio = NULL) {
  # An argument named after ANOTHER penalty's hyperparameter is the mistake
  # this catches: `mcp(x, a = 3)` writes SCAD's shape on an MCP, whose own
  # is gamma, and reaches nothing. R would report it as an unused argument,
  # which says that it was not read and not what to write instead.
  if (length(extra)) {
    pen <- tryCatch(factory(1L), error = function(e) NULL)
    stop(sprintf(paste0("'%s' has no argument '%s'. Its hyperparameters are:",
                        " %s,
  and each is held by naming it and estimated",
                        " when left NULL."),
                 label, names(extra)[[1L]],
                 if (is.null(pen)) "none" else paste(pen@params,
                                                     collapse = ", ")),
         call. = FALSE)
  }
  if (!is.logical(standardize) || length(standardize) != 1L ||
      is.na(standardize)) {
    stop("'standardize' must be TRUE or FALSE.", call. = FALSE)
  }
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
  } else if (isS4(x) && methods::is(x, "Matrix")) {
    # A Matrix is KEPT as it is. Densifying a caller's sparse design here
    # would undo exactly the saving it was passed for, and the factor is
    # 1/density: measured on a 4000 x 60 indicator design at density 0.017,
    # 0.050 MB becomes 1.920 MB. A logical Matrix is carried to double
    # rather than rejected, an indicator being the commonest sparse input.
    if (!methods::is(x, "dMatrix")) {
      x <- tryCatch(methods::as(x, "dMatrix"), error = function(e) {
        stop("a Matrix input must be numeric or logical.", call. = FALSE)
      })
    }
  } else {
    x <- as.matrix(x)
    if (!is.numeric(x)) {
      stop("a matrix input must be numeric.", call. = FALSE)
    }
  }
  PenalizedTerm(label = label, input = x, input_expr = expr,
                factory = factory, standardize = standardize,
                hyper = check_hyper(hyper, factory, label),
                grid = check_grid(grid, factory, label),
                min_ratio = check_min_ratio(min_ratio, label),
                X = NULL, coef_names = character(0),
                blueprint = list(), penalty = NULL)
}

# The spread of each column, read without densifying: for a sparse block
# the two-pass form would have to subtract the mean from every entry, so
# the sum of squares is used there and the two-pass elsewhere, where the
# cancellation it avoids is real for a column with a large mean.
.block_sd <- function(X) {
  if (isS4(X) && methods::is(X, "Matrix")) {
    n <- nrow(X)
    m <- Matrix::colMeans(X)
    s <- sqrt(pmax(Matrix::colSums(X^2) - n * m^2, 0) / max(n - 1L, 1L))
  } else {
    s <- apply(X, 2L, stats::sd)
  }
  s <- as.numeric(s)
  # a constant column has nothing to standardize by, and its coefficient is
  # penalized on its own scale rather than divided by zero
  s[!is.finite(s) | s <= 0] <- 1
  s
}

#' What the Penalized Terms Share
#'
#' @description
#' The input handling, the standardization and the prediction of the five
#' penalized terms, documented once. Each of them -- \code{\link{ridge}},
#' \code{\link{lasso}}, \code{\link{enet}}, \code{\link{scad}}, \code{\link{mcp}} -- has a page
#' of its own carrying its formula, its hyperparameters and where those may
#' lie. Each takes its block as a one-sided formula or as a numeric
#' matrix and attaches the corresponding \pkg{penalties7} object to the
#' block's coefficients at build time, so the hyperparameters, their
#' bounds and links, the derivatives and the kink set are the penalty's,
#' never restated by the term.
#'
#' @details
#' A formula input goes through the \code{\link[stats]{model.matrix}}
#' machinery with the intercept removed (a penalized block does not
#' penalize an intercept; the model's intercept lives in the parametric
#' block), and its blueprint records the terms, the factor levels and the
#' contrasts, exactly as \code{\link{linpar}} does. The exception is a
#' formula whose intercept is all it has: \code{ridge(~1)} is a block of
#' that one column under the penalty, since removing it would leave no
#' block. That is the form a subformula on another term's parameter uses,
#' \code{gamma ~ lasso(~1)} saying that the parameter itself carries a
#' lasso, there being no parametric block of its own to hold an
#' unpenalized intercept. A matrix input is used
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
#' @section Standardization:
#' A hyperparameter is comparable across coordinates only where the
#' coordinates share a scale: without \code{standardize} a lasso penalizes a
#' column measured in metres more than the same column measured in
#' kilometres, and a reader of \eqn{\lambda} has no way to know.
#'
#' \code{standardize = TRUE} divides each coefficient by the standard
#' deviation of its own column, and it does so through the penalty's
#' diagonal map rather than by touching the design. With \eqn{z_j = x_j/s_j}
#' the coefficient satisfies \eqn{\beta_{z,j} = s_j\beta_{x,j}}, so
#'
#' \deqn{\lambda\sum_j \lvert\beta_{z,j}\rvert
#'   = \lambda\sum_j s_j\lvert\beta_{x,j}\rvert
#'   = \rho(S\beta_x), \qquad S = \mathrm{diag}(s),}
#'
#' which is the standardized penalty read on the original scale. Three
#' things follow. The design is never rescaled, so a sparse block stays
#' sparse; \eqn{\lambda} stays one number and the coefficients are already
#' on the scale the data came in, with nothing to map back; and centring,
#' which is what would destroy sparsity, is not needed, the fit being
#' invariant to a translation of a penalized column wherever an intercept
#' is free.
#'
#' The spread is computed from the built block and frozen in the blueprint,
#' so the same term standardizes identically in every equation of a
#' distributional model and does not move with the working weights of a
#' fit. A constant column takes \eqn{s_j = 1}. \code{\link{print}} shows the
#' values, a number that changes the meaning of \eqn{\lambda} having to be
#' legible.
#'
#' For SCAD and MCP the diagonal map is not a rescaling of \eqn{\lambda}
#' alone: substituting \eqn{s_j\beta_j} gives \eqn{\lambda_j = \lambda s_j}
#' AND \eqn{a_j = a/s_j} (or \eqn{\gamma_j = \gamma/s_j}), a composition of
#' both hyperparameters per coordinate, which the map expresses exactly.
#'
#' \code{\link{random}} does not standardize and takes no such argument. Its
#' columns are grouping indicators and its penalty is a variance component
#' with a meaning of its own; weighting it by the size of the groups would
#' change the model rather than its parametrization.
#'
#' @param x A one-sided formula or a numeric matrix.
#' @param label A single non-empty string prefixed to the coefficient
#'   names.
#' @param by Reserved for a later release; must be \code{NULL}.
#' @param standardize A single logical: whether to penalize each
#'   coefficient on the scale of its own column. See the section below.
#' @param ... Not used, and reported: an argument named after another
#'   penalty's hyperparameter is the mistake this catches.
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
#' # the same block penalized on a common scale
#' dd$x3 <- 1000 * dd$x2
#' term_penalty(term_build(lasso(~ x1 + x3, standardize = TRUE), dd))@map
#'
#' @seealso \code{\link{linpar}}, \code{\link{s}}, \code{\link{random}}, \code{\link{term_penalty}}, \code{\link{edf}}
#' @name penalized_terms
#' @keywords internal
NULL

S7::method(term_build, PenalizedTerm) <- function(term, data, ...) {
  if (inherits(term@input, "formula")) {
    # The intercept is dropped because the model's intercept lives in the
    # parametric block and is not penalized -- unless it is all the formula
    # has, where dropping it leaves no block at all. `ridge(~1)` then means
    # the single coefficient under the penalty, which is what a subformula
    # on a term's own parameter needs: `gamma ~ lasso(~1)` says the change
    # itself carries a lasso, there being no parametric block of its own to
    # hold an unpenalized intercept.
    f <- if (.intercept_only(term@input)) term@input else {
      stats::update(term@input, ~ . - 1)
    }
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
  # standardization is a weight per coordinate on the penalty and never an
  # operation on X: with z_j = x_j/s_j the coefficient is beta_z,j = s_j
  # beta_x,j, so rho(s beta) is the standardized penalty read on the
  # original scale. The design keeps its scale and its sparsity, lambda
  # stays one number, and the coefficients need no mapping back.
  map <- NULL
  if (isTRUE(term@standardize)) {
    s <- stats::setNames(.block_sd(X), cn)
    term@blueprint$standardize <- s
    map <- Matrix::Diagonal(x = as.numeric(s))
  }
  # the factory's second argument is the map, and it is passed only when
  # there is one: a subclass written before standardization existed carries a
  # one-argument factory, and it is never standardized, so the old contract
  # goes on working unchanged
  term@penalty <- if (is.null(map)) term@factory(ncol(X)) else
    term@factory(ncol(X), map)
  term
}

# Whether a formula carries an intercept and nothing else, which is the one
# case where removing the intercept would leave no block. `~1` and `~+1`
# reach here; `~0` and `~-1` do not, having neither an intercept nor a term
# and being an empty block whichever rule is applied.
.intercept_only <- function(f) {
  tt <- stats::terms(f)
  attr(tt, "intercept") == 1L && length(attr(tt, "term.labels")) == 0L
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
    # kept sparse for the reason the constructor keeps it: a prediction that
    # densified would spend at new data what the build was careful not to
    X <- tryCatch(.as_block(eval(bp$expr, newdata, baseenv())),
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
