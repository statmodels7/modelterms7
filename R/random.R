#' @include term_classes.R generics.R
NULL

#' @title S7 Class for Grouped Random-Effect Terms
#' @name RandomTerm
#'
#' @description
#' A subclass of \code{\link{additive_term}} for grouped coefficients with
#' a distribution on the effects: the within-group design interacted with
#' the grouping indicators, one coefficient per group and per within-group
#' column, with the penalty carrying the effects' distribution.
#' Constructed by \code{\link{random}}.
#'
#' @inheritParams additive_term
#' @param formula The bar formula, e.g. \code{~ 1 | g} or \code{~ x | g}.
#' @param correlated Logical; whether the default Gaussian lets the
#'   within-group effects correlate.
#' @param precision A \pkg{parameters7} matrix parameter for the
#'   per-group precision, or \code{NULL}.
#' @param distrib A univariate \pkg{distributions7} object for the
#'   effects, or \code{NULL}.
#' @param kinks The declared kink set of \code{distrib}'s log-density.
#'
#' @return An object of class \code{RandomTerm}.
#'
#' @seealso \code{\link{random}}
#' @examples
#' S7::S7_inherits(random(~ 1 | g), RandomTerm)
#' @export
RandomTerm <- S7::new_class(
  name = "RandomTerm",
  parent = additive_term,
  properties = list(
    formula = S7::class_any,
    correlated = S7::class_logical,
    precision = S7::class_any,
    distrib = S7::class_any,
    kinks = S7::class_numeric
  )
)

#' Grouped Random-Effect Term
#'
#' @description
#' Random intercepts and slopes for a grouping factor:
#' \code{random(~ 1 | g)} builds one coefficient per level of \code{g},
#' and \code{random(~ x | g)} one intercept and one slope per level, with
#' the distribution of the effects attached as the penalty on those
#' coefficients -- which is what a random effect is under penalized
#' likelihood.
#'
#' @details
#' The left side of the bar is an ordinary one-sided formula for the
#' within-group design, with the usual intercept convention:
#' \code{~ x | g} carries an intercept and a slope per group and
#' \code{~ 0 + x | g} the slope alone. The block interacts that design
#' with the group indicators, ordered group by group, so the coefficients
#' of one group are adjacent.
#'
#' Three choices of effect distribution are available. By default the
#' effects are Gaussian with one covariance shared by every group:
#' unstructured over the within-group coefficients when
#' \code{correlated = TRUE} (a \code{\link[parameters7]{log_cholesky}}
#' chart, so intercepts and slopes may correlate), diagonal when
#' \code{correlated = FALSE} (independent effects, one variance per
#' within-group column), and a single scale when the bar's left side is
#' the intercept alone. With \code{precision}, a \pkg{parameters7} matrix
#' parameter of the within-group dimension replaces that default (an
#' AR(1) over ordered within-group columns, a compound symmetry); it
#' enters as the per-group precision, replicated across groups by
#' \code{\link[parameters7]{kron_identity}} into
#' \code{\link[penalties7]{structured_penalty}}, and its hyperparameters
#' are the structure's free values. With \code{distrib}, a univariate
#' \pkg{distributions7} object -- holding its own location, typically
#' through \code{\link[distributions7]{fixed}} at zero -- is applied
#' coordinatewise to every effect through
#' \code{\link[penalties7]{distrib_penalty}}. A distribution used as a
#' penalty gives joint-mode (penalized likelihood) estimation of the
#' effects; the marginal likelihood, which integrates them out, is not
#' provided here.
#'
#' Prediction maps new data onto the levels seen at build time; a level
#' the term has not seen is rejected.
#'
#' A random effect is not standardized, and there is no \code{standardize}
#' argument to ask for it with; passing one is an error. The term is a
#' ridge where the effects are Gaussian and independent (see below), but
#' its columns are grouping indicators rather than measured covariates, and
#' its hyperparameter is a variance component with a meaning of its own.
#' Dividing each coefficient by the spread of its indicator would weight
#' the effects by the sizes of the groups, which changes the model rather
#' than the scale its hyperparameter is read on.
#'
#' @section The block and its penalty:
#' With \eqn{m} levels and a within-group design \eqn{Z_i} of \eqn{d}
#' columns, the block is the interaction of that design with the group
#' indicators, ordered group by group,
#'
#' \deqn{Z = \operatorname{diag}(Z_1, \dots, Z_m),
#'   \qquad b = (b_1', \dots, b_m')',}
#'
#' so the \eqn{d} coefficients of one group occupy adjacent positions.
#' Under the default Gaussian the penalty is the negative log-density of
#' \eqn{b} at a precision replicated across groups,
#'
#' \deqn{\rho(b; \eta) = \tfrac{1}{2} b'
#'   \left(I_m \otimes \Omega_d(\eta)\right) b
#'   - \tfrac{m}{2}\log\lvert \Omega_d(\eta)\rvert
#'   + \tfrac{md}{2}\log(2\pi),}
#'
#' whose hyperparameters are the free values \eqn{\eta} of the
#' within-group parameter. Minimizing the penalized least squares in
#' \eqn{(\beta, b)} is the mixed-model equation at the variance ratio the
#' precision encodes, so the minimizer is the best linear unbiased
#' predictor; \eqn{\Omega_d = I_d/\sigma_b^2} recovers the ordinary random
#' intercept, where the penalty is exactly the ridge.
#'
#' @param formula A bar formula, \code{~ lhs | g}, with \code{g}
#'   evaluating to the grouping variable in the data.
#' @param correlated Logical; whether the default Gaussian lets the
#'   within-group effects correlate. Ignored when \code{precision} or
#'   \code{distrib} is given.
#' @param precision A \pkg{parameters7} matrix parameter of the
#'   within-group dimension, or \code{NULL}.
#' @param distrib A univariate \pkg{distributions7} object for the
#'   effects, or \code{NULL}; exclusive with \code{precision}.
#' @param kinks The kink set of \code{distrib}'s log-density in its
#'   argument, passed to \code{\link[penalties7]{distrib_penalty}}.
#' @param label A single non-empty string prefixed to the coefficient
#'   names.
#'
#' @return An object of class \code{\link{RandomTerm}} (a specification;
#'   see \code{\link{term_build}}).
#'
#' @examples
#' dd <- data.frame(y = rnorm(9), x = rnorm(9),
#'                  g = factor(rep(c("a", "b", "c"), 3)))
#' built <- term_build(random(~ x | g), dd)
#' term_coef_names(built)
#' term_penalty(built)@params
#'
#' @references
#' Laird, N. M. and Ware, J. H. (1982). Random-effects models for
#' longitudinal data. \emph{Biometrics} 38, 963-974.
#'
#' @param hyper The hyperparameters of the effects' distribution to HOLD, as
#'   a named vector or list; those not named are estimated. Which names there
#'   are depends on what the term was given: \code{c(lambda = 4)} for the
#'   default, the free names of the structure for a \code{precision}, and the
#'   distribution's own parameters for a \code{distrib}. A name the penalty
#'   does not carry is reported when the term is built, which is where the
#'   penalty first exists.
#'
#' @seealso \code{\link{s}}, \code{\link{te}}, \code{\link{nl}}
#' @export
random <- function(formula, correlated = TRUE, precision = NULL,
                   distrib = NULL, kinks = numeric(0), label = "random",
                   hyper = NULL) {
  if (!inherits(formula, "formula") || length(formula) != 2L) {
    stop("'formula' must be a one-sided bar formula, e.g. ~ 1 | g.",
         call. = FALSE)
  }
  e <- formula[[2L]]
  if (!is.call(e) || !identical(e[[1L]], as.name("|"))) {
    stop("'formula' must contain a grouping bar, e.g. ~ 1 | g.",
         call. = FALSE)
  }
  if (!is.logical(correlated) || length(correlated) != 1L ||
      is.na(correlated)) {
    stop("'correlated' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(precision) && !is.null(distrib)) {
    stop("'precision' and 'distrib' are mutually exclusive.", call. = FALSE)
  }
  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !nzchar(label)) {
    stop("'label' must be a single non-empty character string.",
         call. = FALSE)
  }
  RandomTerm(label = label, formula = formula, correlated = correlated,
             precision = precision, distrib = distrib, kinks = kinks,
             hyper = as_hyper(hyper, label),
             X = NULL, coef_names = character(0),
             blueprint = list(), penalty = NULL)
}

.random_group <- function(expr, data, levels = NULL) {
  v <- eval(expr, data, baseenv())
  if (is.null(levels)) return(factor(v))
  f <- factor(v, levels = levels)
  if (any(is.na(f) & !is.na(v))) {
    bad <- unique(as.character(v)[is.na(f) & !is.na(v)])
    stop(sprintf("grouping level '%s' was not present at build time.",
                 bad[1L]), call. = FALSE)
  }
  f
}

# The within-group design and the indicators, interacted group by group so
# the coefficients of one group are adjacent -- the order I_m %x% S assumes.
#
# The block is SPARSE by construction and is built as such rather than
# assembled dense and converted: a row belongs to one group, so exactly d of
# its m*d entries are non-zero and the density is 1/m whatever the data. The
# dense form was quadratic in the wrong place -- at n = 20000 and m = 1000 it
# is 152.6 MB against 0.23 MB, built in 1.76 s against 0.0011 s, and the
# crossprod every iteration takes is 12.77 s against 0.0006 s. The
# intermediate `outer(g, levels(g), ==)` was itself a dense n x m.
.random_block <- function(g, W) {
  m <- nlevels(g)
  d <- ncol(W)
  n <- nrow(W)
  gi <- as.integer(g)
  # one entry per (row, within-column) pair: row i lands in the d columns of
  # its own group, and in no other
  Matrix::sparseMatrix(
    i = rep(seq_len(n), each = d),
    j = as.vector(t(outer((gi - 1L) * d, seq_len(d), `+`))),
    x = as.vector(t(W)),
    dims = c(n, m * d))
}

.random_names <- function(label, glevels, wnames) {
  if (length(wnames) == 1L && wnames == "(Intercept)") {
    return(paste(label, glevels, sep = "."))
  }
  as.character(t(outer(glevels, wnames,
                       function(a, b) paste(label, a, b, sep = "."))))
}

S7::method(term_build, RandomTerm) <- function(term, data, ...) {
  e <- term@formula[[2L]]
  g <- .random_group(e[[3L]], data)
  m <- nlevels(g)
  if (m < 2L) {
    stop("the grouping variable must have at least two levels.",
         call. = FALSE)
  }

  wf <- stats::as.formula(call("~", e[[2L]]),
                          env = environment(term@formula))
  mf <- stats::model.frame(wf, data, na.action = stats::na.pass,
                           drop.unused.levels = FALSE)
  tt <- attr(mf, "terms")
  W <- stats::model.matrix(tt, mf)
  contr <- attr(W, "contrasts")
  d <- ncol(W)
  if (d < 1L) {
    stop("the left side of the bar yields no columns.", call. = FALSE)
  }

  Z <- .random_block(g, W)
  cn <- .random_names(term@label, levels(g), colnames(W))
  colnames(Z) <- cn

  pen <- if (!is.null(term@precision)) {
    st <- term@precision
    if (!S7::S7_inherits(st, parameters7::matrix_parameter)) {
      stop("'precision' must be a parameters7 matrix_parameter.",
           call. = FALSE)
    }
    if (st@dimension != d) {
      stop(sprintf("'precision' has dimension %d and the bar's left side %d columns.",
                   st@dimension, d), call. = FALSE)
    }
    penalties7::structured_penalty(parameters7::kron_identity(st, m))
  } else if (!is.null(term@distrib)) {
    penalties7::distrib_penalty(term@distrib, n_coef = m * d,
                                kinks = term@kinks)
  } else if (d == 1L) {
    penalties7::ridge_penalty(n_coef = m)
  } else {
    st <- if (term@correlated) parameters7::log_cholesky(d)
          else parameters7::diagonal_matrix(d)
    penalties7::structured_penalty(parameters7::kron_identity(st, m))
  }

  term@X <- Z
  term@coef_names <- cn
  term@blueprint <- list(
    gexpr = e[[3L]], glevels = levels(g),
    terms = stats::delete.response(tt),
    xlev = stats::.getXlevels(tt, mf),
    contrasts = contr
  )
  # the names are checked HERE, the first point at which the penalty this
  # term builds exists: which hyperparameters it has depends on whether
  # the effects carry a structure, a distribution or neither
  term@hyper <- check_hyper(term@hyper, pen, term@label)
  term@penalty <- pen
  term
}

S7::method(term_predict, RandomTerm) <- function(term, newdata, ...) {
  .assert_built(term)
  bp <- term@blueprint
  g <- .random_group(bp$gexpr, newdata, levels = bp$glevels)
  mf <- stats::model.frame(bp$terms, newdata, na.action = stats::na.pass,
                           xlev = bp$xlev)
  W <- stats::model.matrix(bp$terms, mf, contrasts.arg = bp$contrasts)
  Z <- .random_block(g, W)
  colnames(Z) <- term@coef_names
  Z
}
