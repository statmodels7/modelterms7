#' @include term_classes.R generics.R penalized.R
NULL

#' Effective Degrees of Freedom of a Term
#'
#' @description
#' The effective degrees of freedom of a built term, computed from the
#' pieces a fitted model supplies. The counting rule follows the term's
#' penalty. An unpenalized term counts its coefficients exactly. A term
#' with a smooth penalty uses the trace of
#' \eqn{(H + S)^{-1} H} over its block, where \eqn{H} is the term's
#' unpenalized curvature (the weighted crossproduct of its design block at
#' the fit) and \eqn{S} is the penalty's Hessian in the coefficients at
#' the estimated hyperparameters; without a penalty the trace reduces to
#' the coefficient count, and as the penalty grows it falls toward the
#' penalty's null space. A term with a non-smooth penalty counts its
#' nonzero coefficients, which for the lasso is the unbiased estimator of
#' its degrees of freedom [Zou, Hastie & Tibshirani (2007)].
#'
#' @param term A built term (see \code{\link{term_build}}).
#' @param coef The fitted coefficients of the term's block.
#' @param hessian The unpenalized curvature of the fit restricted to the
#'   term's block, a \eqn{k \times k} matrix; required only for a smooth
#'   penalized term.
#' @param theta The estimated hyperparameters of the term's penalty, as a
#'   named list; required only for a smooth penalized term.
#' @param tol The threshold below which a coefficient counts as zero for
#'   a non-smooth penalty.
#' @param ... Passed to methods.
#'
#' @return A single number.
#'
#' @references
#' Zou, H., Hastie, T. and Tibshirani, R. (2007). On the "degrees of
#' freedom" of the lasso. \emph{The Annals of Statistics}, 35(5),
#' 2173--2192.
#'
#' @examples
#' dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
#' built <- term_build(ridge(~ x1 + x2), dd)
#' H <- crossprod(term_matrix(built))
#' edf(built, coef = c(0.5, -0.2), hessian = H, theta = list(sigma = 2))
#'
#' @export
edf <- S7::new_generic("edf", "term",
  function(term, coef = NULL, hessian = NULL, theta = NULL,
           tol = 1e-8, ...) S7::S7_dispatch())

S7::method(edf, additive_term) <- function(term, coef = NULL, hessian = NULL,
                                           theta = NULL, tol = 1e-8, ...) {
  .assert_built(term)
  k <- ncol(term@X)
  pen <- term@penalty
  if (is.null(pen)) {
    return(as.numeric(k))
  }
  if (!is.null(coef) && length(coef) != k) {
    stop(sprintf("'coef' must have length %d.", k), call. = FALSE)
  }
  if (term_smooth(term)) {
    if (is.null(coef) || is.null(hessian) || is.null(theta)) {
      stop("a smooth penalized term needs 'coef', 'hessian' and 'theta'.",
           call. = FALSE)
    }
    hessian <- as.matrix(hessian)
    if (!all(dim(hessian) == k)) {
      stop(sprintf("'hessian' must be a %d x %d matrix.", k, k),
           call. = FALSE)
    }
    S <- penalties7::penalty_hessian(pen, as.numeric(coef), theta)
    return(sum(diag(solve(hessian + S, hessian))))
  }
  if (is.null(coef)) {
    stop("a non-smooth penalized term counts its nonzero coefficients and needs 'coef'.",
         call. = FALSE)
  }
  as.numeric(sum(abs(coef) > tol))
}

# --- printing and plotting --------------------------------------------------

S7::method(print, PenalizedTerm) <- function(x, ...) {
  lab <- if (nzchar(x@label)) sprintf(" '%s'", x@label) else ""
  if (term_is_built(x)) {
    pen <- x@penalty
    cat(sprintf("<PenalizedTerm>%s built: %d coefficient%s; penalty %s (%s)\n",
                lab, ncol(x@X), if (ncol(x@X) == 1L) "" else "s",
                pen@penalty_name, paste(pen@params, collapse = ", ")))
  } else {
    cat(sprintf("<PenalizedTerm>%s (specification; call term_build() with data)\n",
                lab))
  }
  invisible(x)
}

S7::method(plot, additive_term) <- function(x, coef = NULL, ...) {
  .assert_built(x)
  if (is.null(coef)) {
    stop("'coef' is required: a term is displayed at fitted coefficients.",
         call. = FALSE)
  }
  k <- ncol(x@X)
  if (length(coef) != k) {
    stop(sprintf("'coef' must have length %d.", k), call. = FALSE)
  }
  cn <- term_coef_names(x)
  op <- graphics::par(mar = c(7, 4, 3, 1))
  on.exit(graphics::par(op))
  main <- if (nzchar(x@label)) x@label else attr(S7::S7_class(x), "name")
  graphics::plot(seq_len(k), coef, pch = 19, xaxt = "n", xlab = "",
                 ylab = "coefficient", main = main, ...)
  graphics::abline(h = 0, lty = 3)
  graphics::segments(seq_len(k), 0, seq_len(k), coef)
  graphics::axis(1, at = seq_len(k), labels = cn, las = 2, cex.axis = 0.8)
  invisible(x)
}
