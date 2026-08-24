#' @include term_classes.R generics.R penalized.R
NULL

#' Effective Degrees of Freedom of a Term
#'
#' @description
#' Counts what a built term spends, given the coefficients a fit reached, the
#' curvature at that fit and the hyperparameters that were estimated. An
#' unpenalized parameter costs one; a parameter under a penalty costs less, and
#' how much less depends on the kind of penalty. The three rules partition the
#' term's parameters, so a term carrying several penalties is counted piece by
#' piece and the pieces add.
#'
#' The one method is registered on [model_term()] and reads the penalties the
#' term declares through [term_penalties()], so a term class written outside
#' the package is counted by the same rule with nothing to register.
#'
#' @details
#' # The three rules
#'
#' **A parameter no penalty reaches counts one**, exactly. A term with no
#' penalties at all therefore returns `term_npar(term)` and reads none of the
#' other arguments.
#'
#' **A parameter under a non-smooth penalty counts one when it is away from
#' zero and nothing when it is at it.** For the lasso that count is an unbiased
#' estimator of the degrees of freedom (Zou, Hastie and Tibshirani, 2007). A
#' penalty is treated as non-smooth when [penalties7::penalty_kinks()] reports
#' any point at a probe value of its hyperparameters, so lasso, SCAD, MCP and
#' the elastic net go here and ridge does not.
#'
#' **Everything else is counted together by one trace.** Over the parameters
#' left, meaning those unpenalized and those under a smooth penalty,
#'
#' \deqn{\mathrm{edf} = \mathrm{tr}\{(H + S)^{-1} H\},}
#'
#' where \eqn{H} is `hessian` restricted to those rows and columns and \eqn{S}
#' carries each smooth penalty's Hessian in the coefficients, evaluated at
#' `coef` and at the estimated hyperparameters and placed at the parameters
#' that penalty covers. An unpenalized parameter contributes a zero row and
#' column to \eqn{S}, so the trace returns its one. As a smoothing parameter
#' grows the trace falls toward the dimension of the penalty's null space: for
#' [s()] that limit is 1, the straight line the Demmler-Reinsch penalty leaves
#' free.
#'
#' Each rule reduces to what a term reported before it carried more than one
#' penalty: the trace over every column for a single smooth penalty, the
#' nonzero count for a single kinked one, the coefficient count for none.
#'
#' # Which arguments are needed when
#'
#' `hessian` is asked for over the whole block and read at the rows and columns
#' the trace runs over. It is never needed when every parameter sits under a
#' kinked penalty, the count being read from `coef` alone, and `coef` is never
#' needed when the term carries no penalty. Everything the rule in force does
#' not reach may be left `NULL`; leaving out something it does reach throws
#' with the missing arguments named.
#'
#' @param term A built term (see [term_build()]). An unbuilt one throws.
#' @param coef The fitted coefficients of the term's block, a numeric vector of
#'   length `term_npar(term)`. Any other length throws with the required length
#'   named. `NULL` is allowed for an unpenalized term.
#' @param hessian The unpenalized curvature of the fit restricted to the term's
#'   own block, a \eqn{k \times k} matrix with \eqn{k} = `term_npar(term)`,
#'   typically the weighted crossproduct of the block at the fitted weights.
#'   Any other shape throws. `NULL` is allowed unless the trace runs.
#' @param theta The estimated hyperparameters. A term carrying one penalty
#'   takes that penalty's own named list, `list(lambda = 0.25)`. A term
#'   carrying several takes a list of such lists, keyed by the names
#'   [term_penalties()] gives; an unkeyed list throws, naming them. The two
#'   spellings are told apart by the value being a list, a hyperparameter never
#'   being one.
#' @param tol The magnitude below which a coefficient counts as zero under a
#'   non-smooth penalty, `1e-8` by default. A proximal step returns exact
#'   zeros, so the threshold matters only for coefficients a different route
#'   left small.
#' @param ... Passed to methods.
#'
#' @return A single number, between 0 and `term_npar(term)`. Not an integer in
#'   general: only the unpenalized and the non-smooth rules give whole numbers.
#'
#' @references
#' Zou, H., Hastie, T. and Tibshirani, R. (2007). On the "degrees of freedom"
#' of the lasso. *The Annals of Statistics*, 35(5), 2173--2192.
#'
#' @seealso [term_penalties()] for the entries the count runs over,
#'   [term_npar()] for the ceiling, [term_smooth()] for whether a term's
#'   penalized objective is twice differentiable, and
#'   [penalties7::penalty_hessian()] for the \eqn{S} block.
#'
#' @examples
#' set.seed(2)
#' dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20), x3 = rnorm(20))
#'
#' # No penalty: the coefficient count, and nothing else is read.
#' edf(term_build(linpar(~ x1 + x2), dd))
#'
#' # A ridge spends between 2 and 0 as lambda grows.
#' b <- term_build(ridge(~ x1 + x2), dd)
#' H <- crossprod(term_matrix(b))
#' vapply(c(1e-6, 0.25, 10, 1e6),
#'        function(l) edf(b, coef = c(0.5, -0.2), hessian = H,
#'                        theta = list(lambda = l)), numeric(1))
#'
#' # And that is the trace, computed apart.
#' all.equal(edf(b, coef = c(0.5, -0.2), hessian = H, theta = list(lambda = 0.25)),
#'           sum(diag(solve(H + 0.25 * diag(2), H))))
#'
#' # A lasso counts survivors, and needs no hessian.
#' bl <- term_build(lasso(~ x1 + x2 + x3), dd)
#' edf(bl, coef = c(1, 0, -2), theta = list(lambda = 1))
#'
#' # A smooth falls from k toward the dimension of its null space, which
#' # for s() is the one straight line the penalty leaves free.
#' d2 <- data.frame(x = seq(0, 1, length.out = 60))
#' bs <- term_build(s(x, k = 8), d2)
#' Hs <- crossprod(term_matrix(bs))
#' cf <- rnorm(term_npar(bs))
#' c(k = term_npar(bs),
#'   vapply(c(1e-8, 1, 1e10),
#'          function(l) edf(bs, coef = cf, hessian = Hs,
#'                          theta = list(lambda = l)), numeric(1)))
#'
#' # Two penalties on one term: theta is keyed by the entry names.
#' d3 <- data.frame(x = runif(40, 0, 5), g = factor(rep(1:4, 10)))
#' nb <- term_build(nl(~ a * exp(-r * x), a ~ 0 + ridge(~ g),
#'                     r ~ 0 + lasso(~ g)), d3)
#' vapply(term_penalties(nb), function(e) e$name, character(1))
#' Hn <- crossprod(as.matrix(term_matrix(nb)))
#' edf(nb, coef = c(1, 2, 0, 3, 0.5, 0, 0.2, 0), hessian = Hn,
#'     theta = list(`a::ridge(~g)` = list(lambda = 2),
#'                  `r::lasso(~g)` = list(lambda = 1)))
#'
#' @export
#' @aliases edf.model_term
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

#' @title Print a Penalized Term
#' @name print.PenalizedTerm
#'
#' @description
#' Prints one line describing a [ridge()], [lasso()], [enet()], [scad()] or
#' [mcp()] term. An unbuilt specification reports its label and whether it will
#' standardize; a built one reports the number of coefficients, the penalty's
#' name and the hyperparameters that penalty carries, and adds a second line
#' giving the spread each column was standardized by when there is one.
#'
#' @details
#' The built form reads `ncol(x@X)`, `x@penalty@penalty_name` and
#' `x@penalty@params`, so the hyperparameters shown are the penalty's own
#' names, `lambda` for a ridge or a lasso, `lambda, alpha` for the elastic net,
#' `lambda, a` for SCAD and `lambda, gamma` for MCP.
#'
#' The standardization line matters for reading a fitted hyperparameter: with
#' `standardize = TRUE` the penalty acts on the coefficients of columns divided
#' by those spreads, and the spreads are frozen in the blueprint at build time,
#' so prediction at new rows uses the same numbers.
#'
#' @param x A [PenalizedTerm()], built or not.
#' @param ... Unused, and accepted so that the signature matches [print()]'s.
#'
#' @return `x`, invisibly. Called for the line it writes.
#'
#' @seealso [penalized_terms()] for the five constructors and what they share.
#'
#' @examples
#' set.seed(5)
#' d <- data.frame(x1 = rnorm(30), x2 = rnorm(30) * 20)
#'
#' # A specification says only what it is.
#' ridge(~ x1 + x2)
#' lasso(~ x1 + x2, standardize = TRUE)
#'
#' # A built term names its penalty and that penalty's hyperparameters.
#' term_build(ridge(~ x1 + x2), d)
#' term_build(scad(~ x1), d)
#'
#' # Standardizing adds the spreads: x2 was simulated twenty times wider.
#' term_build(enet(~ x1 + x2, standardize = TRUE), d)
#'
#' @keywords internal
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

#' @title Plot an Additive Term's Coefficients
#' @name plot.additive_term
#'
#' @description
#' Draws the coefficients of a built additive term as a stem plot: one point
#' per coefficient against its position in the block, a stem down to zero, a
#' dotted line at zero, and the coefficient names rotated along the horizontal
#' axis. It shows which coefficients a penalty has driven to zero and how large
#' the survivors are, so it reads a lasso or a smooth at a glance.
#'
#' The coefficients are supplied by the caller. A term holds a design block and
#' no fit, so there is nothing to display without them, and `coef = NULL`
#' throws.
#'
#' @details
#' The panel's title is the term's `label` when it has one, so `"ridge"` or
#' `"s(x)"`, and the class name otherwise. The bottom margin is widened to
#' seven lines so that the coefficient names fit, and [graphics::par()] is
#' restored on exit.
#'
#' The horizontal axis is the coefficient's position in the block, so the
#' picture is of the block and not of the covariate. For a smooth, whose block
#' is a Demmler-Reinsch reparametrization, the columns are ordered from least
#' to most wiggly, and a penalized fit shows a decaying profile.
#'
#' @param x A built additive term. An unbuilt one throws
#'   `"the term has not been built; call term_build(term, data) first."`.
#' @param coef The coefficients to draw, a numeric vector of length
#'   `term_npar(x)`. Required: `NULL` throws
#'   `"'coef' is required: a term is displayed at fitted coefficients."`, and
#'   any other length throws with the required length named.
#' @param ... Passed to [graphics::plot()], so `col`, `ylim`, `cex` and the
#'   rest of the usual graphical arguments work. `xaxt`, `xlab`, `ylab`,
#'   `main` and `pch` are set here.
#'
#' @return `x`, invisibly. Called for the plot.
#'
#' @seealso [term_coef_names()] for the axis labels, [edf()] for what the same
#'   coefficients cost in degrees of freedom.
#'
#' @examples
#' set.seed(5)
#' d <- data.frame(x1 = rnorm(30), x2 = rnorm(30), x3 = rnorm(30))
#' b <- term_build(lasso(~ x1 + x2 + x3), d)
#'
#' # Two survivors and one coefficient at exactly zero.
#' plot(b, coef = c(0.8, 0, -0.35))
#'
#' # A smooth's block runs from least to most wiggly.
#' d2 <- data.frame(x = seq(0, 1, length.out = 60))
#' bs <- term_build(s(x, k = 8), d2)
#' plot(bs, coef = c(1.2, 0.9, -0.4, 0.2, -0.1, 0.05, -0.02))
#'
#' # Coefficients are required, and must fit the block.
#' try(plot(b))
#' try(plot(b, coef = 1))
#'
#' @keywords internal
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
