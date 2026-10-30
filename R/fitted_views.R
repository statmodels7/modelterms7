#' @include term_classes.R generics.R penalized.R ggplot.R
NULL

#' Effective Degrees of Freedom of a Term
#'
#' @description
#' Computes the effective degrees of freedom of a built term from the
#' coefficients of a fit, the curvature at that fit and the estimated
#' hyperparameters. An unpenalized parameter counts one. A parameter under a
#' non-smooth penalty counts one when it is nonzero and zero otherwise. The
#' parameters unpenalized or under a smooth penalty are counted together by a
#' single trace, which is less than their number when a penalty acts on them.
#' The three rules partition the term's parameters, and the counts add.
#'
#' The method is registered on [model_term()] and reads the penalties that the
#' term declares through [term_penalties()], so a term class defined outside
#' the package is counted by the same rules without registering a method.
#'
#' @details
#' # The three rules
#'
#' A parameter outside every penalty counts one, exactly. A term with no
#' penalties at all therefore returns `term_npar(term)` and does not read any
#' of the other arguments.
#'
#' A parameter under a non-smooth penalty counts one when it is away from
#' zero and counts zero when it is at zero. For the lasso that count is an
#' unbiased estimator of the degrees of freedom (Zou, Hastie and Tibshirani,
#' 2007). A penalty is treated as non-smooth when [penalties7::penalty_kinks()]
#' reports any point at a probe value of its hyperparameters. The lasso, SCAD,
#' MCP and the elastic net are therefore counted by this rule, and the ridge is
#' not.
#'
#' Everything else is counted together by one trace. Over the parameters
#' that remain, meaning those that are unpenalized and those under a smooth
#' penalty, the count is
#'
#' \deqn{\mathrm{edf} = \mathrm{tr}\{(H + S)^{-1} H\},}
#'
#' where \eqn{H} is `hessian` restricted to those rows and columns and \eqn{S}
#' carries each smooth penalty's Hessian in the coefficients, evaluated at
#' `coef` and at the estimated hyperparameters and placed at the parameters
#' that the penalty covers. An unpenalized parameter contributes a zero row and
#' column to \eqn{S}, so the trace returns its one. As a smoothing parameter
#' grows, the trace falls toward the dimension of the penalty's null space.
#' For [s()] with the B-spline, P-spline or Legendre smoother that limit is 1,
#' the straight line that the Demmler-Reinsch penalty leaves free; for the
#' first-order B-spline, the cyclic and the Fourier smoothers it is 0.
#'
#' When one penalty covers the whole block, the result is the trace over every
#' column if the penalty is smooth, the number of nonzero coefficients if it is
#' non-smooth, and the number of coefficients if there is no penalty.
#'
#' # Which arguments are needed when
#'
#' `hessian` is supplied for the whole block and read at the rows and columns
#' over which the trace runs. It is not needed when every penalty of the term
#' is non-smooth, since the count is then read from `coef` alone, and `coef` is
#' not needed when the term carries no penalty. An argument that the rule in
#' force does not read may be left `NULL`. If an argument that the rule reads
#' is missing, `edf()` signals an error, and the message names the arguments
#' that the rule needs.
#'
#' @param term A built term (see [term_build()]). An unbuilt term signals an
#'   error.
#' @param coef The fitted coefficients of the term's block, a numeric vector of
#'   length `term_npar(term)`. A vector of any other length signals an error
#'   that gives the required length. `NULL` is allowed for an unpenalized term.
#' @param hessian The unpenalized curvature of the fit restricted to the term's
#'   own block, a \eqn{k \times k} matrix with \eqn{k} = `term_npar(term)`,
#'   typically the weighted crossproduct of the block at the fitted weights.
#'   When the trace runs, a matrix of any other shape signals an error. `NULL`
#'   is allowed unless the trace runs.
#' @param theta The estimated hyperparameters. A term carrying one penalty
#'   takes that penalty's own named list, `list(lambda = 0.25)`. A term
#'   carrying several takes a list of such lists, keyed by the names that
#'   [term_penalties()] gives; an unkeyed list signals an error that lists the
#'   names. The two spellings are distinguished by whether the value is a list,
#'   since a hyperparameter is never a list.
#' @param tol The magnitude below which a coefficient counts as zero under a
#'   non-smooth penalty, `1e-8` by default. A proximal step returns exact
#'   zeros, so the threshold matters only for coefficients that a different
#'   route left small.
#' @param ... Passed to methods.
#'
#' @return A single number, between 0 and `term_npar(term)`. The number is not
#'   an integer in general: only the unpenalized and the non-smooth rules give
#'   whole numbers.
#'
#' @references
#' Zou, H., Hastie, T. and Tibshirani, R. (2007). On the "degrees of freedom"
#' of the lasso. *The Annals of Statistics*, 35(5), 2173--2192.
#'
#' @seealso [term_penalties()] for the entries over which the count runs,
#'   [term_npar()] for the ceiling, [term_smooth()] for whether a term's
#'   penalized objective is twice differentiable, and
#'   [penalties7::penalty_hessian()] for the \eqn{S} block.
#'
#' @examples
#' set.seed(2)
#' dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20), x3 = rnorm(20))
#'
#' # Without a penalty, edf() returns the coefficient count and ignores the
#' # other arguments.
#' edf(term_build(linpar(~ x1 + x2), dd))
#'
#' # For a ridge the count falls from 2 toward 0 as lambda grows.
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
#' # For a smooth the count falls from k toward the dimension of the penalty's
#' # null space, which for this smoother is the straight line.
#' d2 <- data.frame(x = seq(0, 1, length.out = 60))
#' bs <- term_build(s(x, basis7::bspline_smooth(k = 8)), d2)
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
# handed that penalty's own list; a term carrying several is handed one list
# per penalty name, and the two spellings are distinguished by the value being
# a list, a hyperparameter never being one.
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
#' name and the hyperparameters that the penalty carries, and, when the term is
#' standardized, adds a second line giving the spread by which each column is
#' standardized.
#'
#' @details
#' The built form reads `ncol(x@X)`, `x@penalty@penalty_name` and
#' `x@penalty@params`, so the hyperparameters shown are the penalty's own
#' names: `lambda` for a ridge or a lasso, `lambda, alpha` for the elastic net,
#' `lambda, a` for SCAD and `lambda, gamma` for MCP.
#'
#' The standardization line is needed to read a fitted hyperparameter. With
#' `standardize = TRUE` the penalty is evaluated at the coefficients multiplied
#' by those spreads; these products are the coefficients of the columns after
#' each is divided by its spread. The spreads are computed once at build time
#' and recorded in the `build_info`. They enter the penalty only: the design block
#' is not rescaled, so prediction at new rows does not use them.
#'
#' @param x A [PenalizedTerm()], built or not.
#' @param ... Unused, and accepted so that the signature matches [print()]'s.
#'
#' @return `x`, invisibly. The method is called for the line that it writes.
#'
#' @seealso [penalized_terms()] for the five constructors and what they share.
#'
#' @examples
#' set.seed(5)
#' d <- data.frame(x1 = rnorm(30), x2 = rnorm(30) * 20)
#'
#' # A specification prints what it is, without any quantity from data.
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
    s <- x@build_info$standardize
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
#' `autoplot()` returns a graph of the coefficients of a built additive term
#' as a stem plot, and `plot()` prints that graph: one point per coefficient
#' against its position in the block, a stem down to zero, a grey line at
#' zero, and the coefficient names along the horizontal axis. It shows which
#' coefficients a penalty has driven to zero and how large the survivors are,
#' so the effect of a lasso or a smooth is visible at a glance.
#'
#' The coefficients are supplied by the caller. A term holds a design block and
#' no fit, so there is nothing to display without them, and `coef = NULL` is
#' rejected.
#'
#' @details
#' The title of the graph is the term's `label` when it has one, so `"ridge"`
#' or `"s(x)"`, and the class name otherwise.
#'
#' The horizontal axis is the coefficient's position in the block, so the
#' graph shows the block and not the covariate. For [s()], whose block is a
#' Demmler-Reinsch reparametrization, the columns are ordered from the smoothest
#' to the most wiggly, and a penalized fit shows a decaying profile.
#'
#' The graph is drawn with \pkg{ggplot2}, which is suggested; the methods
#' signal an error that names it when it is not installed. The graph has one
#' element, the coefficients, which the vertical axis names, and so it has no
#' legend.
#'
#' @param x,object A built additive term. An unbuilt term signals the error
#'   `"the term has not been built; call term_build(term, data) first."`.
#' @param coef The coefficients to draw, a numeric vector of length
#'   `term_npar(x)`. Required: `NULL` signals the error
#'   `"'coef' is required: a term is displayed at fitted coefficients."`, and a
#'   vector of any other length signals an error that gives the required
#'   length.
#' @param ... Unused; an argument passed here signals an error.
#'
#' @return `autoplot()` returns a ggplot2 object; `plot()` prints the same
#'   object and returns it invisibly.
#'
#' @seealso [term_coef_names()] for the axis labels, [edf()] for what the same
#'   coefficients cost in degrees of freedom.
#'
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' set.seed(5)
#' d <- data.frame(x1 = rnorm(30), x2 = rnorm(30), x3 = rnorm(30))
#' b <- term_build(lasso(~ x1 + x2 + x3), d)
#'
#' # The graph shows two survivors and one coefficient at exactly zero.
#' plot(b, coef = c(0.8, 0, -0.35))
#'
#' # A smooth's block runs from least to most wiggly.
#' d2 <- data.frame(x = seq(0, 1, length.out = 60))
#' bs <- term_build(s(x, basis7::bspline_smooth(k = 8)), d2)
#' plot(bs, coef = c(1.2, 0.9, -0.4, 0.2, -0.1, 0.05, -0.02))
#'
#' # Coefficients are required, and must fit the block.
#' try(plot(b))
#' try(plot(b, coef = 1))
#'
#' @keywords internal
S7::method(plot, additive_term) <- function(x, coef = NULL, ...) {
  p <- autoplot_additive_term(x, coef = coef, ...)
  print(p)
  invisible(p)
}

#' @rdname plot.additive_term
autoplot_additive_term <- function(object, coef = NULL, ...) {
  if (...length()) {
    stop("The graph of a term takes no further arguments; modify the ",
         "returned object instead.", call. = FALSE)
  }
  x <- object
  .assert_built(x)
  if (is.null(coef)) {
    stop("'coef' is required: a term is displayed at fitted coefficients.",
         call. = FALSE)
  }
  k <- ncol(x@X)
  if (length(coef) != k) {
    stop(sprintf("'coef' must have length %d.", k), call. = FALSE)
  }
  need_ggplot()
  cn <- term_coef_names(x)
  d <- data.frame(pos = seq_len(k), coef = as.numeric(coef))
  main <- if (nzchar(x@label)) x@label else attr(S7::S7_class(x), "name")
  .data <- ggplot2::.data
  ggplot2::ggplot(d, ggplot2::aes(x = .data$pos, y = .data$coef)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$pos, yend = 0),
                          colour = "#1f5fa8", linewidth = 0.7) +
    ggplot2::geom_point(colour = "#1f5fa8", size = 2) +
    ggplot2::scale_x_continuous(breaks = seq_len(k), labels = cn,
                                minor_breaks = NULL) +
    ggplot2::labs(title = main, x = NULL, y = "coefficient") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                       hjust = 1,
                                                       vjust = 0.5))
}

S7::method(autoplot, additive_term) <- autoplot_additive_term
