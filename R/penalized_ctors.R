#' @include penalized.R
NULL

# One page per penalty. The five share their input handling, their
# standardization and their prediction, which are documented once on
# penalized_terms; what differs is the function each one attaches, the
# hyperparameters it carries and where those may lie, and that is what a
# reader comes to one of these pages for.

#' Ridge Penalty on a Block of Coefficients
#'
#' @description
#' A block of coefficients under a Gaussian prior at zero: every one of them
#' shrunk towards zero, none of them set to it.
#'
#' @details
#' Writing \eqn{\beta} for the block's coefficients and \eqn{p} for their
#' number,
#' \deqn{\rho(\beta) = \frac{\lambda\lVert\beta\rVert_2^2}{2}
#'   - \frac{p}{2}\log\!\left(\frac{\lambda}{2\pi}\right).}
#' The constant is kept, which is what makes \eqn{\lambda} estimable by a
#' marginal criterion: it is minus the log density of
#' \eqn{N(0, \lambda^{-1}I)}, so \eqn{\lambda} is the PRECISION of that
#' prior and a larger value shrinks harder.
#'
#' The penalty is twice differentiable everywhere, so the block is fitted in
#' the same system as the unpenalized terms and \eqn{\lambda} is estimated by
#' \code{\link[statmodels7]{reml}()} rather than swept along a path.
#'
#' \strong{Hyperparameter.} \code{lambda}, admissible on \eqn{(0, \infty)}.
#'
#' @inheritParams penalized_terms
#' @param lambda The precision of the prior. One number holds it and
#'   \code{NULL}, the default, has it ESTIMATED. A ridge has no kink and no
#'   path, so several numbers are not a grid it could visit. Must lie in
#'   \eqn{(0, \infty)}.
#' @param n_lambda Unused by this term: a ridge has no kink and its
#'   hyperparameter is estimated by a criterion rather than swept over a
#'   grid. Accepted so that the five constructors read alike.
#'
#' @param min_ratio Unused by this term, which has no path. Accepted so
#'   that the five constructors read alike.
#' @return An object of class \code{\link{PenalizedTerm}} (a specification;
#'   see \code{\link{term_build}}).
#'
#' @examples
#' dd <- data.frame(x1 = rnorm(8), x2 = rnorm(8))
#' term_penalty(term_build(ridge(~ x1 + x2), dd))@params
#' term_hyper(ridge(~ x1 + x2, lambda = 2))
#'
#' @references
#' Hoerl, A. E. and Kennard, R. W. (1970). Ridge regression: biased
#' estimation for nonorthogonal problems. \emph{Technometrics} 12, 55--67.
#'
#' @seealso \code{\link{penalized_terms}} for what the five share,
#'   \code{\link{lasso}}, \code{\link{enet}}, \code{\link{scad}},
#'   \code{\link{mcp}}, \code{\link[penalties7]{ridge_penalty}}
#' @export
ridge <- function(x, label = "ridge", by = NULL, standardize = FALSE,
                  lambda = NULL, n_lambda = NULL, min_ratio = NULL, ...) {
  .penalized_spec(x, substitute(x), label, by, standardize,
                  function(k, map = NULL) penalties7::ridge_penalty(map = map,
                                                                    n_coef = k),
                  list(lambda = lambda),
                  list(...), list(lambda = n_lambda), min_ratio)
}


#' Lasso Penalty on a Block of Coefficients
#'
#' @description
#' A block of coefficients under a Laplace prior at zero: the penalty has a
#' kink there, so coefficients are set EXACTLY to zero and the term selects.
#'
#' @details
#' \deqn{\rho(\beta) = \lambda\lVert\beta\rVert_1
#'   - p\log\!\left(\frac{\lambda}{2}\right).}
#' The constant is kept, so this is minus the log density of a Laplace at
#' zero with rate \eqn{\lambda}, and a larger \eqn{\lambda} shrinks harder
#' and keeps fewer coefficients.
#'
#' The kink is at zero, so the block is fitted by a proximal method or by a
#' coordinate descent with the other terms held, and \eqn{\lambda} is chosen
#' by a PATH over its own values -- \code{\link[statmodels7]{bic}()} by
#' default, or \code{\link[statmodels7]{aic}()} or
#' \code{\link[statmodels7]{cv}()} -- because a marginal criterion is a
#' Laplace expansion at a mode that sits on the kink.
#'
#' \strong{Hyperparameter.} \code{lambda}, admissible on \eqn{(0, \infty)},
#' swept over \code{n_lambda} values from the one that empties the block
#' down to \code{min_ratio} of it.
#'
#' @inheritParams penalized_terms
#' @param lambda The rate of the prior. One number holds it, several are the
#'   grid the path visits as they stand, and \code{NULL}, the default, has the
#'   path build one. Must lie in \eqn{(0, \infty)}.
#' @param n_lambda How many values the path visits, at least 2. \code{NULL},
#'   the default, leaves it to the criterion.
#'
#' @param min_ratio How far down the path reaches, as a fraction of the
#'   kink that empties the block: smaller reaches a denser fit, larger
#'   stops sooner. Must lie in (0, 1). NULL, the default, leaves it to the
#'   criterion. Only the sweep by kink size uses it.
#' @return An object of class \code{\link{PenalizedTerm}} (a specification;
#'   see \code{\link{term_build}}).
#'
#' @examples
#' dd <- data.frame(x1 = rnorm(8), x2 = rnorm(8))
#' built <- term_build(lasso(~ x1 + x2), dd)
#' term_penalty(built)@params
#' term_smooth(built)
#'
#' # a finer path for a wide block
#' term_grid(lasso(~ x1 + x2, n_lambda = 60))
#'
#' @references
#' Tibshirani, R. (1996). Regression shrinkage and selection via the lasso.
#' \emph{Journal of the Royal Statistical Society, Series B} 58, 267--288.
#'
#' @seealso \code{\link{penalized_terms}} for what the five share,
#'   \code{\link{ridge}}, \code{\link{enet}}, \code{\link{scad}},
#'   \code{\link{mcp}}, \code{\link[penalties7]{lasso_penalty}}
#' @export
lasso <- function(x, label = "lasso", by = NULL, standardize = FALSE,
                  lambda = NULL, n_lambda = NULL, min_ratio = NULL, ...) {
  .penalized_spec(x, substitute(x), label, by, standardize,
                  function(k, map = NULL) penalties7::lasso_penalty(map = map,
                                                                    n_coef = k),
                  list(lambda = lambda),
                  list(...), list(lambda = n_lambda), min_ratio)
}


#' Elastic Net Penalty on a Block of Coefficients
#'
#' @description
#' The lasso and the ridge mixed: a kink at zero, so the term still selects,
#' and a quadratic part that keeps correlated coefficients together instead
#' of choosing arbitrarily among them.
#'
#' @details
#' \deqn{\rho(\beta) = \lambda\left\{\alpha\lVert\beta\rVert_1
#'   + \frac{1-\alpha}{2}\lVert\beta\rVert_2^2\right\}
#'   + p\log Z(\lambda, \alpha),}
#' the normalizing constant being that of the product of a Laplace and a
#' Gaussian at zero
#' (\code{\link[distributions7]{enet_distrib}}). It depends on BOTH
#' hyperparameters, which is what makes them estimable rather than merely
#' settable, and what a penalty written as a formula would not have.
#'
#' \eqn{\alpha} is the mixing weight: at \eqn{\alpha \to 1} the penalty is
#' the lasso and at \eqn{\alpha \to 0} the ridge, and the kink at zero has
#' half-width \eqn{\lambda\alpha}, so both hyperparameters scale it.
#'
#' \strong{Hyperparameters.} \code{lambda} on \eqn{(0, \infty)}, swept over
#' \code{n_lambda} values by kink size; \code{alpha} on \eqn{(0, 1)}, swept
#' over \code{n_alpha} values across that interval, the ends excluded
#' because the penalty there is one of the other two. The two are swept
#' cyclically, one at a time, which keeps the cost linear in their number
#' where a product grid would be exponential in it.
#'
#' @inheritParams penalized_terms
#' @param lambda The overall rate. One number holds it, several are the grid
#'   the path visits as they stand, and \code{NULL}, the default, has the path
#'   build one. Must lie in \eqn{(0, \infty)}.
#' @param alpha The mixing weight, in the same three states and settled
#'   independently of \code{lambda}. Must lie in \eqn{(0, 1)}.
#' @param n_lambda,n_alpha How many values the path visits for each, at
#'   least 2. \code{NULL}, the default, leaves it to the criterion.
#'
#' @param min_ratio How far down the path reaches, as a fraction of the
#'   kink that empties the block: smaller reaches a denser fit, larger
#'   stops sooner. Must lie in (0, 1). NULL, the default, leaves it to the
#'   criterion. Only the sweep by kink size uses it.
#' @return An object of class \code{\link{PenalizedTerm}} (a specification;
#'   see \code{\link{term_build}}).
#'
#' @examples
#' dd <- data.frame(x1 = rnorm(8), x2 = rnorm(8))
#' term_penalty(term_build(enet(~ x1 + x2), dd))@params
#'
#' # alpha held at the halfway mixture, lambda still estimated
#' term_hyper(enet(~ x1 + x2, alpha = 0.5))
#'
#' @references
#' Zou, H. and Hastie, T. (2005). Regularization and variable selection via
#' the elastic net. \emph{Journal of the Royal Statistical Society, Series B}
#' 67, 301--320.
#'
#' @seealso \code{\link{penalized_terms}} for what the five share,
#'   \code{\link{ridge}}, \code{\link{lasso}}, \code{\link{scad}},
#'   \code{\link{mcp}}, \code{\link[penalties7]{elasticnet_penalty}}
#' @export
enet <- function(x, label = "enet", by = NULL, standardize = FALSE,
                 lambda = NULL, alpha = NULL,
                 n_lambda = NULL, n_alpha = NULL, min_ratio = NULL, ...) {
  .penalized_spec(x, substitute(x), label, by, standardize,
                  function(k, map = NULL)
                    penalties7::elasticnet_penalty(map = map, n_coef = k),
                  list(lambda = lambda, alpha = alpha),
                  list(...), list(lambda = n_lambda, alpha = n_alpha),
                  min_ratio)
}


#' SCAD Penalty on a Block of Coefficients
#'
#' @description
#' Smoothly clipped absolute deviation: the lasso's kink at zero, so the
#' term selects, and a penalty that FLATTENS beyond a threshold, so a large
#' coefficient is not shrunk at all.
#'
#' @details
#' It is defined by its derivative rather than by its value, for
#' \eqn{t = \lvert\beta_j\rvert \ge 0},
#' \deqn{\rho'(t) = \lambda\min\!\left\{1,
#'   \frac{(a\lambda - t)_+}{(a-1)\lambda}\right\},}
#' summed over the coefficients. It rises like the lasso near zero, bends
#' from \eqn{t = \lambda}, and is flat past \eqn{t = a\lambda}: the bias the
#' lasso puts on a large coefficient is what this removes. Being improper it
#' carries no normalizing constant, and is therefore not a log prior and not
#' reachable by a marginal criterion.
#'
#' \strong{Hyperparameters.} \code{lambda} on \eqn{(0, \infty)}, swept over
#' \code{n_lambda} values by kink size; \code{a} on \eqn{(2, \infty)} --
#' below 2 the penalty is not what its definition intends -- swept over
#' \code{n_a} values on a geometric grid above that bound, since the shape
#' leaves the kink at zero unchanged and no kink-size path can reach it. The
#' literature's value is \eqn{a = 3.7}, and holding it there is what
#' \pkg{ncvreg} does.
#'
#' @inheritParams penalized_terms
#' @param lambda The scale of the penalty. One number holds it, several are
#'   the grid the path visits as they stand, and \code{NULL}, the default, has
#'   the path build one. Must lie in \eqn{(0, \infty)}.
#' @param a The shape, in the same three states and settled independently of
#'   \code{lambda}. Must lie in \eqn{(2, \infty)}.
#' @param n_lambda,n_a How many values the path visits for each, at least 2.
#'   \code{NULL}, the default, leaves it to the criterion.
#'
#' @param min_ratio How far down the path reaches, as a fraction of the
#'   kink that empties the block: smaller reaches a denser fit, larger
#'   stops sooner. Must lie in (0, 1). NULL, the default, leaves it to the
#'   criterion. Only the sweep by kink size uses it.
#' @return An object of class \code{\link{PenalizedTerm}} (a specification;
#'   see \code{\link{term_build}}).
#'
#' @examples
#' dd <- data.frame(x1 = rnorm(8), x2 = rnorm(8))
#' term_penalty(term_build(scad(~ x1 + x2), dd))@params
#' term_hyper(scad(~ x1 + x2, a = 3.7))
#'
#' @references
#' Fan, J. and Li, R. (2001). Variable selection via nonconcave penalized
#' likelihood and its oracle properties. \emph{Journal of the American
#' Statistical Association} 96, 1348--1360.
#'
#' @seealso \code{\link{penalized_terms}} for what the five share,
#'   \code{\link{ridge}}, \code{\link{lasso}}, \code{\link{enet}},
#'   \code{\link{mcp}}, \code{\link[penalties7]{scad_penalty}}
#' @export
scad <- function(x, label = "scad", by = NULL, standardize = FALSE,
                 lambda = NULL, a = NULL,
                 n_lambda = NULL, n_a = NULL, min_ratio = NULL, ...) {
  .penalized_spec(x, substitute(x), label, by, standardize,
                  function(k, map = NULL) penalties7::scad_penalty(map = map,
                                                                   n_coef = k),
                  list(lambda = lambda, a = a),
                  list(...), list(lambda = n_lambda, a = n_a), min_ratio)
}


#' MCP Penalty on a Block of Coefficients
#'
#' @description
#' The minimax concave penalty: like SCAD it selects and then flattens, and
#' it begins to flatten IMMEDIATELY rather than after a first threshold.
#'
#' @details
#' Defined by its derivative, for \eqn{t = \lvert\beta_j\rvert \ge 0},
#' \deqn{\rho'(t) = \left(\lambda - \frac{t}{\gamma}\right)_+ ,}
#' summed over the coefficients: it starts at \eqn{\lambda}, falls linearly,
#' and is flat past \eqn{t = \gamma\lambda}. Improper by construction, so it
#' carries no normalizing constant and is not reachable by a marginal
#' criterion.
#'
#' \strong{Hyperparameters.} \code{lambda} on \eqn{(0, \infty)}, swept over
#' \code{n_lambda} values by kink size; \code{gamma} on \eqn{(1, \infty)} --
#' at \eqn{\gamma \le 1} the penalized objective need not be convex even for
#' an orthogonal design -- swept over \code{n_gamma} values on a geometric
#' grid above that bound. The literature's value is \eqn{\gamma = 3}.
#'
#' @inheritParams penalized_terms
#' @param lambda The scale of the penalty. One number holds it, several are
#'   the grid the path visits as they stand, and \code{NULL}, the default, has
#'   the path build one. Must lie in \eqn{(0, \infty)}.
#' @param gamma The shape, in the same three states and settled independently
#'   of \code{lambda}. Must lie in \eqn{(1, \infty)}.
#' @param n_lambda,n_gamma How many values the path visits for each, at
#'   least 2. \code{NULL}, the default, leaves it to the criterion.
#'
#' @param min_ratio How far down the path reaches, as a fraction of the
#'   kink that empties the block: smaller reaches a denser fit, larger
#'   stops sooner. Must lie in (0, 1). NULL, the default, leaves it to the
#'   criterion. Only the sweep by kink size uses it.
#' @return An object of class \code{\link{PenalizedTerm}} (a specification;
#'   see \code{\link{term_build}}).
#'
#' @examples
#' dd <- data.frame(x1 = rnorm(8), x2 = rnorm(8))
#' term_penalty(term_build(mcp(~ x1 + x2), dd))@params
#' term_hyper(mcp(~ x1 + x2, gamma = 3))
#'
#' @references
#' Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax
#' concave penalty. \emph{The Annals of Statistics} 38, 894--942.
#'
#' @seealso \code{\link{penalized_terms}} for what the five share,
#'   \code{\link{ridge}}, \code{\link{lasso}}, \code{\link{enet}},
#'   \code{\link{scad}}, \code{\link[penalties7]{mcp_penalty}}
#' @export
mcp <- function(x, label = "mcp", by = NULL, standardize = FALSE,
                lambda = NULL, gamma = NULL,
                n_lambda = NULL, n_gamma = NULL, min_ratio = NULL, ...) {
  .penalized_spec(x, substitute(x), label, by, standardize,
                  function(k, map = NULL) penalties7::mcp_penalty(map = map,
                                                                  n_coef = k),
                  list(lambda = lambda, gamma = gamma),
                  list(...), list(lambda = n_lambda, gamma = n_gamma),
                  min_ratio)
}
