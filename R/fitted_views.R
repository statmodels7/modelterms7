#' @include term_classes.R generics.R penalized.R
NULL

#' Effective Degrees of Freedom of a Term
#'
#' @description
#' The effective degrees of freedom of a built term, computed from the
#' pieces a fitted model supplies. The counting rule follows the penalties
#' the term declares through [term_penalties()], and applies to
#' each of them over the parameters it covers.
#'
#' @details
#' A parameter no penalty reaches counts one, exactly. A parameter under a
#' **non-smooth** penalty counts one when it is away from zero and
#' nothing when it is at it, which for the lasso is the unbiased estimator
#' of its degrees of freedom (Zou, Hastie & Tibshirani, 2007). The
#' remaining parameters -- those unpenalized and those under a
#' **smooth** penalty -- are counted together by the trace of
#' \eqn{(H + S)^{-1} H} over the sub-block they occupy, where \eqn{H} is
#' the term's unpenalized curvature there (the weighted crossproduct of its
#' design block at the fit) and \eqn{S} carries each smooth penalty's
#' Hessian in the coefficients at the estimated hyperparameters, placed at
#' the parameters that penalty covers and zero elsewhere. An unpenalized
#' parameter contributes a zero row and column to \eqn{S}, so the trace
#' returns its one; as a penalty grows the trace falls toward the dimension
#' of its null space.
#'
#' The rules compose because they partition the term's parameters, and each
#' reduces to what the term reported before when one penalty covers the
#' whole block: the trace over every column for a smooth penalty, the
#' nonzero count for a kinked one, the coefficient count for none.
#'
#' `hessian` is asked for over the whole block, and is used at the rows
#' and columns the trace runs over. It is not needed at all when every
#' penalty is kinked, since the count is then read from `coef` alone.
#'
#' @param term A built term (see [term_build()]).
#' @param coef The fitted coefficients of the term's block.
#' @param hessian The unpenalized curvature of the fit restricted to the
#'   term's block, a \eqn{k \times k} matrix; required whenever some
#'   parameter is not under a kinked penalty.
#' @param theta The estimated hyperparameters. For a term carrying one
#'   penalty, that penalty's hyperparameters as a named list; for a term
#'   carrying several, a list of such lists keyed by the penalty names
#'   [term_penalties()] gives.
#' @param tol The threshold below which a coefficient counts as zero for
#'   a non-smooth penalty.
#' @param ... Passed to methods.
#'
#' @return A single number.
#'
#' @references
#' Zou, H., Hastie, T. and Tibshirani, R. (2007). On the "degrees of
#' freedom" of the lasso. *The Annals of Statistics*, 35(5),
#' 2173--2192.
#'
#' @examples
#' dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
#' built <- term_build(ridge(~ x1 + x2), dd)
#' H <- crossprod(term_matrix(built))
#' edf(built, coef = c(0.5, -0.2), hessian = H, theta = list(lambda = 0.25))
#'
#' @seealso [term_penalties()], [term_penalty()], [term_smooth()]
#' @export
edf <- S7::new_generic("edf", "term",
  function(term, coef = NULL, hessian = NULL, theta = NULL,
           tol = 1e-8, ...) S7::S7_dispatch())

# The hyperparameters of one entry. A term carrying a single penalty is
# handed that penalty's own list, as it always was; a term carrying several
# is handed one list per penalty name, and the two spellings are told apart
# by the value being a list, a hyperparameter never being one.
.edf_theta <- function(entries, e, theta) {
  if (nzchar(e$name) && is.list(theta) && is.list(theta[[e$name]])) {
    return(theta[[e$name]])
  }
  if (length(entries) > 1L) {
    stop(sprintf(paste("'theta' must be keyed by the penalty names (%s):",
                       "this term carries %d penalties."),
                 paste(vapply(entries, function(z) z$name, character(1)),
                       collapse = ", "), length(entries)),
         call. = FALSE)
  }
  theta
}

S7::method(edf, model_term) <- function(term, coef = NULL, hessian = NULL,
                                        theta = NULL, tol = 1e-8, ...) {
  k <- term_npar(term)
  ent <- term_penalties(term)
  if (!length(ent)) {
    return(as.numeric(k))
  }
  if (!is.null(coef) && length(coef) != k) {
    stop(sprintf("'coef' must have length %d.", k), call. = FALSE)
  }
  kinked <- vapply(ent, function(e) {
    length(penalties7::penalty_kinks(e$penalty,
                                     .penalty_probe_theta(e$penalty))) > 0L
  }, logical(1))

  hard <- sort(unique(unlist(lapply(ent[kinked], function(e) e$index))))
  count <- 0
  if (length(hard)) {
    if (is.null(coef)) {
      stop("a non-smooth penalty counts nonzero coefficients and needs 'coef'.",
           call. = FALSE)
    }
    count <- sum(abs(coef[hard]) > tol)
  }
  rest <- setdiff(seq_len(k), hard)
  if (!length(rest)) {
    return(as.numeric(count))
  }
  soft <- ent[!kinked]
  if (!length(soft)) {
    return(as.numeric(count + length(rest)))
  }
  if (is.null(coef) || is.null(hessian) || is.null(theta)) {
    stop("a smooth penalty needs 'coef', 'hessian' and 'theta'.",
         call. = FALSE)
  }
  hessian <- as.matrix(hessian)
  if (!all(dim(hessian) == k)) {
    stop(sprintf("'hessian' must be a %d x %d matrix.", k, k), call. = FALSE)
  }
  S <- matrix(0, k, k)
  for (e in soft) {
    i <- e$index
    S[i, i] <- S[i, i] +
      penalties7::penalty_hessian(e$penalty, as.numeric(coef[i]),
                                  .edf_theta(ent, e, theta))
  }
  H <- hessian[rest, rest, drop = FALSE]
  count + sum(diag(solve(H + S[rest, rest, drop = FALSE], H)))
}

# --- printing and plotting --------------------------------------------------

S7::method(print, PenalizedTerm) <- function(x, ...) {
  lab <- if (nzchar(x@label)) sprintf(" '%s'", x@label) else ""
  if (term_is_built(x)) {
    pen <- x@penalty
    cat(sprintf("<PenalizedTerm>%s built: %d coefficient%s; penalty %s (%s)\n",
                lab, ncol(x@X), if (ncol(x@X) == 1L) "" else "s",
                pen@penalty_name, paste(pen@params, collapse = ", ")))
    # the scale a hyperparameter is read against belongs on the page
    s <- x@blueprint$standardize
    if (!is.null(s)) {
      cat("  standardized by: ",
          paste(sprintf("%s = %.4g", names(s), s), collapse = ", "),
          "\n", sep = "")
    }
  } else {
    cat(sprintf("<PenalizedTerm>%s%s (specification; call term_build() with data)\n",
                lab, if (isTRUE(x@standardize)) ", standardized" else ""))
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
