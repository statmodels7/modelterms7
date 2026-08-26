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
#' The constant is kept, and that is what a marginal criterion needs to
#' estimate \eqn{\lambda}: the value is minus the log density of
#' \eqn{N(0, \lambda^{-1}I)}, so \eqn{\lambda} is the precision of that
#' prior and a larger value shrinks harder.
#'
#' The penalty is twice differentiable everywhere, so the block is fitted in
#' the same system as the unpenalized terms and \eqn{\lambda} is estimated by
#' [statmodels7::reml()] rather than swept along a path.
#'
#' **Hyperparameter.** `lambda`, admissible on \eqn{(0, \infty)}.
#'
#' @inheritParams penalized_terms
#' @param label A single non-empty character string prefixed to the
#'   coefficient names as `label.name`, `"ridge"` by default, so a block
#'   over `x1` reads `ridge.x1`. Two penalized terms in one formula stay
#'   apart by their labels.
#' @param lambda The precision of the prior. One number holds it and
#'   `NULL`, the default, has it estimated. A ridge has no kink and no
#'   path, so several numbers are not a grid it could visit. Must lie in
#'   \eqn{(0, \infty)}.
#' @return An unbuilt [PenalizedTerm()]: a specification, with `X`,
#'   `coef_names`, `blueprint` and `penalty` empty until [term_build()]
#'   fills them, and the penalty attached there over as many coefficients
#'   as the block turns out to have.
#'
#' @examples
#' set.seed(3)
#' dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
#' b <- term_build(ridge(~ x1 + x2), dd)
#'
#' # One hyperparameter, positive, on the log scale for an optimizer.
#' p <- term_penalty(b)
#' p@params
#' p@params_bounds
#'
#' # The value is exactly minus a Gaussian log-density at precision lambda.
#' beta <- c(0.4, -1.1)
#' all.equal(penalties7::penalty_value(p, beta, list(lambda = 2.5)),
#'           -sum(dnorm(beta, 0, 1 / sqrt(2.5), log = TRUE)))
#'
#' # No kink, so the block is smooth and lambda is estimated at the mode.
#' c(smooth = term_smooth(b),
#'   kinks = length(penalties7::penalty_kinks(p, list(lambda = 1))))
#'
#' # Holding it, and the refusal of a grid it has no path to visit.
#' term_hyper(ridge(~ x1 + x2, lambda = 2))
#' try(ridge(~ x1 + x2, lambda = c(1, 2)))
#'
#'
#' # Fitted. The data are simulated from a known truth, so the
#' # estimates below can be read against it.
#' if (requireNamespace("statmodels7", quietly = TRUE)) {
#'   set.seed(11)
#'   XX <- matrix(rnorm(150 * 8), 150, 8)
#'   fd <- data.frame(y = as.numeric(XX %*% c(2, -1.5, 1, rep(0, 5))) +
#'                      rnorm(150, sd = 0.4))
#'   fd$X <- XX
#'   cf <- coef(statmodels7::statmod(y ~ ridge(X),
#'                                   distributions7::gaussian1_distrib(), fd))$mu
#'   # truth: the first three columns carry 2, -1.5 and 1, the other five
#'   # nothing. A ridge shrinks and keeps every column.
#'   round(cf[2:5], 2)
#' }
#' @references
#' Hoerl, A. E. and Kennard, R. W. (1970). Ridge regression: biased
#' estimation for nonorthogonal problems. *Technometrics* 12, 55--67.
#'
#' @seealso [penalized_terms()] for what the five share,
#'   [lasso()], [enet()], [scad()],
#'   [mcp()], [penalties7::ridge_penalty()]
#' @export
ridge <- function(x, label = "ridge", standardize = FALSE,
                  lambda = NULL, sparse = NULL, ...) {
  .penalized_spec(x, substitute(x), label, standardize,
                  function(k, map = NULL) penalties7::ridge_penalty(map = map,
                                                                    n_coef = k),
                  list(lambda = lambda), list(...), sparse = sparse)
}


#' Lasso Penalty on a Block of Coefficients
#'
#' @description
#' A block of coefficients under a Laplace prior at zero: the penalty has a
#' kink there, so coefficients are set exactly to zero and the term selects.
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
#' by a path over its own values, scored by [statmodels7::bic()] by default
#' or by [statmodels7::aic()] or [statmodels7::cv()]. A marginal criterion
#' cannot be used: it is a Laplace expansion at a mode that sits on the kink.
#'
#' **Hyperparameter.** `lambda`, admissible on \eqn{(0, \infty)},
#' swept over `n_lambda` values from the one that empties the block
#' down to `min_ratio` of it.
#'
#' @inheritParams penalized_terms
#' @param label A single non-empty character string prefixed to the
#'   coefficient names as `label.name`, `"lasso"` by default, so a block
#'   over `x1` reads `lasso.x1`. Two penalized terms in one formula stay
#'   apart by their labels.
#' @param lambda The rate of the prior. One number holds it, several are the
#'   grid the path visits as they stand, and `NULL`, the default, has the
#'   path build one. Must lie in \eqn{(0, \infty)}.
#' @param n_lambda How many values the path visits, a whole number of at
#'   least 2, `25` by default. The axis descends four decades of kink
#'   size, and that many points are what covers it.
#' @param min_ratio How far down the path reaches, as a fraction of the
#'   kink that empties the block: smaller reaches a denser fit, larger
#'   stops sooner. A single number in \eqn{(0, 1)}, `1e-4` by default.
#'   Only the sweep by kink size reads it.
#' @return An unbuilt [PenalizedTerm()]: a specification, with `X`,
#'   `coef_names`, `blueprint` and `penalty` empty until [term_build()]
#'   fills them, and the penalty attached there over as many coefficients
#'   as the block turns out to have.
#'
#' @examples
#' set.seed(3)
#' dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
#' b <- term_build(lasso(~ x1 + x2), dd)
#' p <- term_penalty(b)
#'
#' # The kink is at zero, so the block is not smooth and needs a path.
#' c(smooth = term_smooth(b))
#' penalties7::penalty_kinks(p, list(lambda = 1))
#'
#' # The value is minus a Laplace log-density at rate lambda.
#' beta <- c(0.4, -1.1)
#' all.equal(penalties7::penalty_value(p, beta, list(lambda = 2)),
#'           2 * sum(abs(beta)) - 2 * log(2 / 2))
#'
#' # The path's length and depth, at the defaults and set.
#' term_grid(lasso(~ x1 + x2))
#' term_grid(lasso(~ x1 + x2, n_lambda = 60))
#' term_path_min(lasso(~ x1 + x2, min_ratio = 1e-3))
#'
#'
#' # Fitted. The data are simulated from a known truth, so the
#' # estimates below can be read against it.
#' if (requireNamespace("statmodels7", quietly = TRUE)) {
#'   set.seed(11)
#'   XX <- matrix(rnorm(150 * 8), 150, 8)
#'   fd <- data.frame(y = as.numeric(XX %*% c(2, -1.5, 1, rep(0, 5))) +
#'                      rnorm(150, sd = 0.4))
#'   fd$X <- XX
#'   cf <- coef(statmodels7::statmod(y ~ lasso(X),
#'                                   distributions7::gaussian1_distrib(), fd))$mu
#'   # truth: the first three columns carry 2, -1.5 and 1, the other five
#'   # nothing. The hyperparameter is chosen by BIC.
#'   c(round(cf[2:5], 2), kept = sum(cf[-1] != 0))
#' }
#' @references
#' Tibshirani, R. (1996). Regression shrinkage and selection via the lasso.
#' *Journal of the Royal Statistical Society, Series B* 58, 267--288.
#'
#' @seealso [penalized_terms()] for what the five share,
#'   [ridge()], [enet()], [scad()],
#'   [mcp()], [penalties7::lasso_penalty()]
#' @export
lasso <- function(x, label = "lasso", standardize = FALSE,
                  lambda = NULL, n_lambda = 25, min_ratio = 1e-4,
                  sparse = NULL, ...) {
  .penalized_spec(x, substitute(x), label, standardize,
                  function(k, map = NULL) penalties7::lasso_penalty(map = map,
                                                                    n_coef = k),
                  list(lambda = lambda),
                  list(...), list(lambda = n_lambda), min_ratio,
                  sparse = sparse)
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
#' ([distributions7::enet_distrib()]). It depends on both
#' hyperparameters, so both are estimable, where a merely settable one would
#' be all a dropped constant leaves. A penalty
#' written as a formula, with the constant dropped, would not have that.
#'
#' \eqn{\alpha} is the mixing weight: at \eqn{\alpha \to 1} the penalty is
#' the lasso and at \eqn{\alpha \to 0} the ridge, and the kink at zero has
#' half-width \eqn{\lambda\alpha}, so both hyperparameters scale it.
#'
#' **Hyperparameters.** `lambda` on \eqn{(0, \infty)}, swept over
#' `n_lambda` values by kink size; `alpha` on \eqn{(0, 1)}, swept
#' over `n_alpha` values across that interval, the ends excluded
#' because the penalty there is one of the other two. Every combination of
#' the two is visited, `n_lambda * n_alpha` fits, unless
#' `search = "cyclic"` asks for one at a time instead.
#'
#' @inheritParams penalized_terms
#' @param label A single non-empty character string prefixed to the
#'   coefficient names as `label.name`, `"enet"` by default, so a block
#'   over `x1` reads `enet.x1`. Two penalized terms in one formula stay
#'   apart by their labels.
#' @param lambda The overall rate. One number holds it, several are the grid
#'   the path visits as they stand, and `NULL`, the default, has the path
#'   build one. Must lie in \eqn{(0, \infty)}.
#' @param alpha The mixing weight, in the same three states and settled
#'   independently of `lambda`. Must lie in \eqn{(0, 1)}.
#' @param n_lambda,n_alpha How many values the path visits for each, at
#'   least 2. They differ because the axes do: \eqn{\lambda} descends the
#'   size of the kink over four decades and wants that many points, while
#'   \eqn{\alpha} spans one bounded interval and does not.
#' @param search `"grid"` to visit every combination of \eqn{\lambda}
#'   and \eqn{\alpha}, `"cyclic"` to sweep one at a time with the other
#'   held. See [term_search()].
#'
#' @param min_ratio How far down the path reaches, as a fraction of the
#'   kink that empties the block: smaller reaches a denser fit, larger
#'   stops sooner. A single number in \eqn{(0, 1)}, `1e-4` by default.
#'   Only the sweep by kink size reads it.
#' @return An unbuilt [PenalizedTerm()]: a specification, with `X`,
#'   `coef_names`, `blueprint` and `penalty` empty until [term_build()]
#'   fills them, and the penalty attached there over as many coefficients
#'   as the block turns out to have.
#'
#' @examples
#' set.seed(3)
#' dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
#' pe <- term_penalty(term_build(enet(~ x1 + x2), dd))
#' pe@params
#' pe@params_bounds
#'
#' # The two ends really are the other two penalties.
#' pl <- term_penalty(term_build(lasso(~ x1 + x2), dd))
#' pr <- term_penalty(term_build(ridge(~ x1 + x2), dd))
#' beta <- c(0.4, -1.1)
#' all.equal(penalties7::penalty_value(pe, beta, list(lambda = 2, alpha = 1 - 1e-9)),
#'           penalties7::penalty_value(pl, beta, list(lambda = 2)))
#' all.equal(penalties7::penalty_value(pe, beta, list(lambda = 2, alpha = 1e-9)),
#'           penalties7::penalty_value(pr, beta, list(lambda = 2)))
#'
#' # alpha held at the halfway mixture, lambda still estimated.
#' term_hyper(enet(~ x1 + x2, alpha = 0.5))
#'
#' # The grid is 25 by 5 by default; cyclic sweeps one axis at a time.
#' term_grid(enet(~ x1 + x2))
#' term_search(enet(~ x1 + x2, search = "cyclic"))
#'
#' # The open interval is enforced: an end is one of the other penalties.
#' try(enet(~ x1 + x2, alpha = 1))
#'
#'
#' # Fitted. The data are simulated from a known truth, so the
#' # estimates below can be read against it.
#' if (requireNamespace("statmodels7", quietly = TRUE)) {
#'   set.seed(11)
#'   XX <- matrix(rnorm(150 * 8), 150, 8)
#'   fd <- data.frame(y = as.numeric(XX %*% c(2, -1.5, 1, rep(0, 5))) +
#'                      rnorm(150, sd = 0.4))
#'   fd$X <- XX
#'   cf <- coef(statmodels7::statmod(y ~ enet(X, n_alpha = 3),
#'                                   distributions7::gaussian1_distrib(), fd))$mu
#'   # truth: the first three columns carry 2, -1.5 and 1, the other five
#'   # nothing. Both hyperparameters are chosen by BIC, over a
#'   # shortened alpha axis so that the product grid stays quick.
#'   c(round(cf[2:5], 2), kept = sum(cf[-1] != 0))
#' }
#' @references
#' Zou, H. and Hastie, T. (2005). Regularization and variable selection via
#' the elastic net. *Journal of the Royal Statistical Society, Series B*
#' 67, 301--320.
#'
#' @seealso [penalized_terms()] for what the five share,
#'   [ridge()], [lasso()], [scad()],
#'   [mcp()], [penalties7::elasticnet_penalty()]
#' @export
enet <- function(x, label = "enet", standardize = FALSE,
                 lambda = NULL, alpha = NULL,
                 n_lambda = 25, n_alpha = 5, min_ratio = 1e-4,
                 search = "grid", sparse = NULL, ...) {
  .penalized_spec(x, substitute(x), label, standardize,
                  function(k, map = NULL)
                    penalties7::elasticnet_penalty(map = map, n_coef = k),
                  list(lambda = lambda, alpha = alpha),
                  list(...), list(lambda = n_lambda, alpha = n_alpha),
                  min_ratio, search, sparse)
}


#' SCAD Penalty on a Block of Coefficients
#'
#' @description
#' Smoothly clipped absolute deviation: the lasso's kink at zero, so the
#' term selects, and a penalty that flattens beyond a threshold, so a large
#' coefficient is not shrunk at all.
#'
#' @details
#' It is defined by its derivative, for \eqn{t = \lvert\beta_j\rvert \ge 0},
#' \deqn{\rho'(t) = \lambda\min\!\left\{1,
#'   \frac{(a\lambda - t)_+}{(a-1)\lambda}\right\},}
#' summed over the coefficients. It rises like the lasso near zero, bends
#' from \eqn{t = \lambda}, and is flat past \eqn{t = a\lambda}, which removes
#' the bias the lasso puts on a large coefficient. Being improper it
#' carries no normalizing constant, and is therefore not a log prior and not
#' reachable by a marginal criterion.
#'
#' **Hyperparameters.** `lambda` on \eqn{(0, \infty)}, swept over
#' `n_lambda` values by kink size; `a` on \eqn{(2, \infty)}, below which the
#' penalty is not what its definition intends, swept over
#' `n_a` values on a geometric grid above that bound, since the shape
#' leaves the kink at zero unchanged and no kink-size path can reach it.
#'
#' **The shape is HELD at \eqn{a = 3.7} by default**, the value of Fan and
#' Li and the one \pkg{ncvreg} holds, so only \eqn{\lambda} is searched.
#' `a = NULL` estimates it over the grid instead, and `n_a` sizes that
#' grid. Measured over eight data configurations, estimating it chose the
#' grid's LOWER ENDPOINT every time and changed no model: the columns kept
#' were identical and the error against the truth agreed to four decimals.
#'
#' ⚠️ **The shape is not scale free**, which is why the literature's
#' value belongs to standardized data. See [mcp()], whose page carries the
#' measurement: rescaling the response reproduces the fit only when the
#' shape moves with the square of the scale.
#'
#' @inheritParams penalized_terms
#' @param label A single non-empty character string prefixed to the
#'   coefficient names as `label.name`, `"scad"` by default, so a block
#'   over `x1` reads `scad.x1`. Two penalized terms in one formula stay
#'   apart by their labels.
#' @param lambda The scale of the penalty. One number holds it, several are
#'   the grid the path visits as they stand, and `NULL`, the default, has
#'   the path build one. Must lie in \eqn{(0, \infty)}.
#' @param a The shape, in the same three states and settled independently
#'   of `lambda`, defaulting to the literature's `3.7` rather than to
#'   `NULL`: one number holds it, several are the grid the path visits as
#'   they stand, and `NULL` has the path build one. Must lie in
#'   \eqn{(2, \infty)}.
#' @param n_lambda,n_a How many values the path visits for each, at least 2.
#'   They differ because the axes do: \eqn{\lambda} descends the size of the
#'   kink over four decades and wants that many points, while \eqn{a} spans
#'   the shape's useful range and does not.
#' @param search `"grid"` to visit every combination of \eqn{\lambda}
#'   and \eqn{a}, `"cyclic"` to sweep one at a time with the other held.
#'   See [term_search()].
#'
#' @param min_ratio How far down the path reaches, as a fraction of the
#'   kink that empties the block: smaller reaches a denser fit, larger
#'   stops sooner. A single number in \eqn{(0, 1)}, `1e-4` by default.
#'   Only the sweep by kink size reads it.
#' @return An unbuilt [PenalizedTerm()]: a specification, with `X`,
#'   `coef_names`, `blueprint` and `penalty` empty until [term_build()]
#'   fills them, and the penalty attached there over as many coefficients
#'   as the block turns out to have.
#'
#' @examples
#' set.seed(3)
#' dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
#' p <- term_penalty(term_build(scad(~ x1 + x2), dd))
#' p@params
#' p@params_bounds
#'
#' # The derivative is the definition: lasso-like to lambda, then bending,
#' # then flat past a * lambda.
#' rho1 <- function(t, lambda, a)
#'   lambda * pmin(1, pmax(a * lambda - t, 0) / ((a - 1) * lambda))
#' rho1(c(0.5, 1.5, 4), lambda = 1, a = 3.7)
#'
#' # penalty_kinks() reports every point where some derivative breaks:
#' # the origin, and the two pairs where the second derivative changes
#' # branch, at plus and minus lambda and a * lambda.
#' penalties7::penalty_kinks(p, list(lambda = 1, a = 3))
#'
#' # The literature's shape, held.
#' term_hyper(scad(~ x1 + x2, a = 3.7))
#'
#' # Below 2 the shape is refused.
#' try(scad(~ x1 + x2, a = 2))
#'
#'
#' # Fitted. The data are simulated from a known truth, so the
#' # estimates below can be read against it.
#' if (requireNamespace("statmodels7", quietly = TRUE)) {
#'   set.seed(11)
#'   XX <- matrix(rnorm(150 * 8), 150, 8)
#'   fd <- data.frame(y = as.numeric(XX %*% c(2, -1.5, 1, rep(0, 5))) +
#'                      rnorm(150, sd = 0.4))
#'   fd$X <- XX
#'   cf <- coef(statmodels7::statmod(y ~ scad(X),
#'                                   distributions7::gaussian1_distrib(), fd))$mu
#'   # truth: the first three columns carry 2, -1.5 and 1, the other five
#'   # nothing. Only lambda is searched: the shape is held at the
#'   # literature's 3.7, and `a = NULL` would estimate it instead.
#'   c(round(cf[2:5], 2), kept = sum(cf[-1] != 0))
#' }
#' @references
#' Fan, J. and Li, R. (2001). Variable selection via nonconcave penalized
#' likelihood and its oracle properties. *Journal of the American
#' Statistical Association* 96, 1348--1360.
#'
#' @seealso [penalized_terms()] for what the five share,
#'   [ridge()], [lasso()], [enet()],
#'   [mcp()], [penalties7::scad_penalty()]
#' @export
scad <- function(x, label = "scad", standardize = FALSE,
                 lambda = NULL, a = 3.7,
                 n_lambda = 25, n_a = 5, min_ratio = 1e-4,
                 search = "grid", sparse = NULL, ...) {
  .penalized_spec(x, substitute(x), label, standardize,
                  function(k, map = NULL) penalties7::scad_penalty(map = map,
                                                                   n_coef = k),
                  list(lambda = lambda, a = a),
                  list(...), list(lambda = n_lambda, a = n_a), min_ratio,
                  search, sparse)
}


#' MCP Penalty on a Block of Coefficients
#'
#' @description
#' The minimax concave penalty: like SCAD it selects and then flattens, and
#' it begins to flatten immediately, where SCAD waits for a first threshold.
#'
#' @details
#' Defined by its derivative, for \eqn{t = \lvert\beta_j\rvert \ge 0},
#' \deqn{\rho'(t) = \left(\lambda - \frac{t}{\gamma}\right)_+ ,}
#' summed over the coefficients: it starts at \eqn{\lambda}, falls linearly,
#' and is flat past \eqn{t = \gamma\lambda}. Improper by construction, so it
#' carries no normalizing constant and is not reachable by a marginal
#' criterion.
#'
#' **Hyperparameters.** `lambda` on \eqn{(0, \infty)}, swept over
#' `n_lambda` values by kink size; `gamma` on \eqn{(1, \infty)}, at or below
#' which the penalized objective need not be convex even for an orthogonal
#' design, swept over `n_gamma` values on a geometric grid above that bound.
#'
#' **The shape is HELD at \eqn{\gamma = 3} by default**, the value of
#' Zhang, so only \eqn{\lambda} is searched. `gamma = NULL` estimates it
#' over the grid instead, and `n_gamma` sizes that grid. Measured over eight
#' data configurations, estimating it chose the grid's LOWER ENDPOINT every
#' time and changed no model: the columns kept were identical and the error
#' against the truth agreed to four decimals. The endpoint is the floor plus
#' 0.25, so what was reported as an estimate was the grid's own edge, and it
#' did not move with `n_gamma`, which refines the interior and cannot touch
#' either end.
#'
#' ⚠️ **The shape is not scale free**, which is why the literature's
#' value belongs to standardized data. Rescaling the response by \eqn{k}
#' sends \eqn{\beta \to k\beta}, and the fit is reproduced at
#' \eqn{(\lambda/k,\ \gamma k^2)}, verified to 1.7e-13 at
#' \eqn{k = 2, 5, 10}: the shape carries the units of a proximal step, and
#' the floor it is swept above grows with them. On a response of much larger
#' scale a small held shape falls outside the convex region the compiled
#' coordinate descent needs and the fit takes the general route, measured at
#' 7.4 s against 0.4 s at a hundredfold response.
#'
#' @inheritParams penalized_terms
#' @param label A single non-empty character string prefixed to the
#'   coefficient names as `label.name`, `"mcp"` by default, so a block
#'   over `x1` reads `mcp.x1`. Two penalized terms in one formula stay
#'   apart by their labels.
#' @param lambda The scale of the penalty. One number holds it, several are
#'   the grid the path visits as they stand, and `NULL`, the default, has
#'   the path build one. Must lie in \eqn{(0, \infty)}.
#' @param gamma The shape, in the same three states and settled
#'   independently of `lambda`, defaulting to the literature's `3` rather
#'   than to `NULL`: one number holds it, several are the grid the path
#'   visits as they stand, and `NULL` has the path build one. Must lie in
#'   \eqn{(1, \infty)}.
#' @param n_lambda,n_gamma How many values the path visits for each, at
#'   least 2. They differ because the axes do: \eqn{\lambda} descends the
#'   size of the kink over four decades and wants that many points, while
#'   \eqn{\gamma} spans the shape's useful range and does not.
#' @param search `"grid"` to visit every combination of \eqn{\lambda}
#'   and \eqn{\gamma}, `"cyclic"` to sweep one at a time with the other
#'   held. See [term_search()].
#'
#' @param min_ratio How far down the path reaches, as a fraction of the
#'   kink that empties the block: smaller reaches a denser fit, larger
#'   stops sooner. A single number in \eqn{(0, 1)}, `1e-4` by default.
#'   Only the sweep by kink size reads it.
#' @return An unbuilt [PenalizedTerm()]: a specification, with `X`,
#'   `coef_names`, `blueprint` and `penalty` empty until [term_build()]
#'   fills them, and the penalty attached there over as many coefficients
#'   as the block turns out to have.
#'
#' @examples
#' set.seed(3)
#' dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
#' p <- term_penalty(term_build(mcp(~ x1 + x2), dd))
#' p@params
#' p@params_bounds
#'
#' # The derivative starts at lambda and falls linearly to zero at
#' # gamma * lambda, where SCAD would still be on its first segment.
#' rho1 <- function(t, lambda, gamma) pmax(lambda - t / gamma, 0)
#' rho1(c(0, 1, 2, 3), lambda = 1, gamma = 2)
#'
#' # The kink at zero, and the pair where the second derivative breaks.
#' penalties7::penalty_kinks(p, list(lambda = 1, gamma = 2))
#'
#' # The literature's shape, held.
#' term_hyper(mcp(~ x1 + x2, gamma = 3))
#'
#' # At or below 1 the objective need not be convex, so it is refused.
#' try(mcp(~ x1 + x2, gamma = 1))
#'
#'
#' # Fitted. The data are simulated from a known truth, so the
#' # estimates below can be read against it.
#' if (requireNamespace("statmodels7", quietly = TRUE)) {
#'   set.seed(11)
#'   XX <- matrix(rnorm(150 * 8), 150, 8)
#'   fd <- data.frame(y = as.numeric(XX %*% c(2, -1.5, 1, rep(0, 5))) +
#'                      rnorm(150, sd = 0.4))
#'   fd$X <- XX
#'   cf <- coef(statmodels7::statmod(y ~ mcp(X),
#'                                   distributions7::gaussian1_distrib(), fd))$mu
#'   # truth: the first three columns carry 2, -1.5 and 1, the other five
#'   # nothing. Only lambda is searched: the shape is held at the
#'   # literature's 3, and `gamma = NULL` would estimate it instead.
#'   c(round(cf[2:5], 2), kept = sum(cf[-1] != 0))
#' }
#' @references
#' Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax
#' concave penalty. *The Annals of Statistics* 38, 894--942.
#'
#' @seealso [penalized_terms()] for what the five share,
#'   [ridge()], [lasso()], [enet()],
#'   [scad()], [penalties7::mcp_penalty()]
#' @export
mcp <- function(x, label = "mcp", standardize = FALSE,
                lambda = NULL, gamma = 3,
                n_lambda = 25, n_gamma = 5, min_ratio = 1e-4,
                search = "grid", sparse = NULL, ...) {
  .penalized_spec(x, substitute(x), label, standardize,
                  function(k, map = NULL) penalties7::mcp_penalty(map = map,
                                                                  n_coef = k),
                  list(lambda = lambda, gamma = gamma),
                  list(...), list(lambda = n_lambda, gamma = n_gamma),
                  min_ratio, search, sparse)
}
