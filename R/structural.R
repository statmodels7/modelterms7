#' @include term_classes.R generics.R
NULL

#' @title Parameters of a Structural Term
#'
#' @description
#' Returns the names of a structural term's own parameters, in the order that
#' its filter, its likelihood and their derivatives use. A structural term
#' contributes no design block, so these are not coefficients: they are
#' estimated beside the distribution's, on the unconstrained scale that
#' [term_links()] defines, and [term_npar()] counts them.
#'
#' @details
#' The names are the term's own vocabulary. [gas()] gives `omega` for the
#' level, `kappa1` ... `kappap` for the score loadings and `pacf1` ... `pacfq`
#' for the persistence. [regime()] gives `level1`, the gaps `gap2` ...
#' `gapk` between consecutive levels, and the free entries `alr*.*` of its
#' transition matrix. A marginal break-point term gives the names listed in
#' [MarginalBreakTerm()]. The names index everything else about the
#' term: [term_start()] returns one value per name, [term_readable()] maps
#' them onto the quantities that a reader reads, and the `index` of a
#' [term_penalties()] entry gives positions in this vector.
#'
#' A subformula expands in place the parameter that it develops, so
#' `gas(p = 1, q = 1)` has three parameters and
#' `gas(p = 1, q = 1, omega ~ z)` has four: `omega.(Intercept)` and `omega.z`
#' where the level was.
#'
#' The method on [structural_term()] signals an error that names the class.
#' A structural class supplies its own method.
#'
#' @param term An object inheriting from [structural_term()]. A class that
#'   does not implement the generic signals the error
#'   `"the term class 'X' does not implement term_params()."`.
#' @param ... Passed to methods. The shipped methods do not use these
#'   arguments.
#'
#' @return A character vector, one name per parameter, of length
#'   [term_npar()].
#'
#' @seealso [term_links()] for the chart on which each parameter is carried,
#'   [term_start()] for where the parameters begin, [term_readable()] for the
#'   quantities to which they map, [term_coef_names()] for the additive
#'   branch's equivalent.
#'
#' @examples
#' # The score-driven vocabulary consists of a level, a loading and a
#' # persistence.
#' term_params(gas(p = 1, q = 1))
#' term_params(gas(p = 2, q = 2))
#'
#' # A subformula expands the parameter that it develops, in place.
#' set.seed(1)
#' d <- data.frame(y = rnorm(30), z = rnorm(30), t = 1:30)
#' term_params(term_build(gas(p = 1, q = 1, omega ~ z, time = t), d))
#'
#' # On this branch term_npar() counts the names that term_params() returns.
#' g <- gas(p = 1, q = 2)
#' c(npar = term_npar(g), names = length(term_params(g)))
#'
#' @export
#' @aliases term_params.structural_term
term_params <- S7::new_generic("term_params", "term",
  function(term, ...) S7::S7_dispatch())

#' @title Links of a Structural Term's Parameters
#'
#' @description
#' Returns one \pkg{linkfunctions7} link per parameter of [term_params()],
#' carrying that parameter from its own admissible set onto the whole real
#' line. A fitting layer optimizes on that unconstrained scale, so the
#' optimizer needs no information about the constraint.
#'
#' @details
#' The charts are chosen so that a proposal from anywhere lands somewhere
#' admissible. [gas()] puts its level on the identity, its score loadings on
#' the log so a loading stays positive, and its persistence coordinates on
#' the rhobit so each partial autocorrelation stays inside \eqn{(-1, 1)}
#' and the filter stays stationary. [regime()] puts its first level and the
#' free entries of its transition matrix on the identity and the gaps between
#' levels on the log, so the levels are ordered by construction. The `links`
#' argument of [gas()] overrides the links per parameter.
#'
#' The persistence uses a different kind of chart from the other parameters.
#' The stationary region in the autoregressive coefficients is not a box, so
#' a collection of scalar links cannot cover it; the partial
#' autocorrelations are each in \eqn{(-1, 1)} independently, and
#' Levinson-Durbin carries them onto the coefficients. The function
#' [term_readable()] exists because the coordinate and the quantity that a
#' reader reads are different things.
#'
#' Where a parameter carries a subformula the link is applied inside the
#' development, \eqn{\psi_{j,t} = g_j^{-1}(z_t'\gamma_j)}, so the parameter
#' stays admissible at every observation and the coefficients \eqn{\gamma_j}
#' are unconstrained.
#'
#' The method on [structural_term()] signals an error that names the class.
#'
#' @param term An object inheriting from [structural_term()]. A class that
#'   does not implement the generic signals the error
#'   `"the term class 'X' does not implement term_links()."`.
#' @param ... Passed to methods. The shipped methods do not use these
#'   arguments.
#'
#' @return A named list of \pkg{linkfunctions7} link objects, one per entry of
#'   [term_params()] and named by it.
#'
#' @seealso [term_params()] for the names, [term_start()] for the point on
#'   this scale at which a fit begins, [term_readable()] for the quantities
#'   to which the coordinates map.
#'
#' @examples
#' # The level uses the identity link, the loading uses the log link and the
#' # persistence uses the rhobit link.
#' vapply(term_links(gas(p = 1, q = 2)), function(l) l@link_name, character(1))
#'
#' # Each carries its parameter's own set onto the whole line.
#' lk <- term_links(gas(p = 1, q = 1))
#' vapply(lk, function(l) paste(l@link_bounds, collapse = ", "), character(1))
#'
#' # So any coordinate at all gives an admissible parameter. At a
#' # coordinate of 40 the persistence prints as 1 and is not: the
#' # inverse link is clamped strictly inside its bounds.
#' rho <- linkfunctions7::linkinv(lk$pacf1, c(-40, 0, 40))
#' all(rho > -1 & rho < 1)
#' 1 - rho[3]
#' linkfunctions7::linkinv(lk$kappa1, c(-40, 0, 40))
#'
#' @export
#' @aliases term_links.structural_term
term_links <- S7::new_generic("term_links", "term",
  function(term, ...) S7::S7_dispatch())

#' @title Where a Term's Own Parameters Start
#'
#' @description
#' Returns the starting values of a structural term's parameters, on the
#' unconstrained scale that [term_links()] defines, one per name that
#' [term_params()] gives. A fitting layer reads it to begin a search.
#'
#' @details
#' The start is supplied by the term because the meaning of a coordinate of
#' zero depends on the chart that the term uses for each parameter. The base
#' method returns zero everywhere, which is each link's own natural point, and
#' that is appropriate wherever zero means "the model without this term".
#'
#' Zero does not always have that meaning. [gas()] overrides the base method:
#' its score loadings ride a log chart, so a coordinate of zero is a loading of
#' one, which is a strongly driven filter and a poor place to begin. The
#' method places each loading at 0.1 through the chart, \eqn{\log 0.1 = -2.303},
#' and leaves every other coordinate at zero. When a loading is developed by a
#' subformula, the constant start is projected onto the design of the
#' development.
#'
#' A term that needs the response to place its start reads it from `target`,
#' which reaches the method through the dots. [regime()] places its levels at
#' quantiles of `target`, and a marginal break-point term computes the profile
#' of its break-points from `target` and otherwise uses quantiles of the
#' covariate. Both are computed when the generic is called.
#'
#' @param term A built structural term.
#' @param ... Passed to methods. The methods for [regime()] and for a marginal
#'   break-point term read `target` here, the response on the scale of the
#'   predictor.
#'
#' @return A named numeric vector on the unconstrained scale, of length
#'   [term_npar()] and named as [term_params()].
#'
#' @seealso [term_links()] for the scale on which the start is expressed,
#'   [term_params()] for the names, [term_coef_start()] for the additive
#'   branch's equivalent.
#'
#' @examples
#' # The level and the persistence start at zero; the loading does not.
#' term_start(gas(p = 1, q = 1))
#'
#' # The loading starts at 0.1 on the parameter scale, because zero on a log
#' # chart is a loading of one.
#' lk <- term_links(gas(p = 1, q = 1))
#' linkfunctions7::linkinv(lk$kappa1, term_start(gas(p = 1, q = 1))[["kappa1"]])
#'
#' # The result has one value per parameter, whatever the order.
#' g <- gas(p = 2, q = 2)
#' identical(names(term_start(g)), term_params(g))
#'
#' @export
#' @aliases term_start.structural_term
term_start <- S7::new_generic("term_start", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_start, structural_term) <- function(term, ...) {
  stats::setNames(numeric(length(term_params(term))), term_params(term))
}

#' @title Drawing a Term's Own Parameters
#'
#' @description
#' Draws a point from the region that a structural term's parameters
#' plausibly occupy. The point consists of the values that [term_start()]
#' gives, displaced by normal noise on the unconstrained scale that
#' [term_links()] defines.
#'
#' @details
#' Simulating from a model that carries state requires values for the term's
#' own parameters. A score-driven filter whose level is developed over a
#' hundred groups has a hundred level coordinates in addition to its loading
#' and persistence coordinates, and choosing each of them by hand is
#' impractical.
#'
#' The rule is the start plus noise. The start is used instead of zero because
#' it lies in the region of the term's own charts where the term is well
#' behaved. A score loading rides a log chart, so zero there is a loading of
#' one and a filter strong enough to destabilize its own recursion;
#' [term_start()] places it at 0.1 and the draw spreads from that. A
#' constraint that a chart carries is respected by every draw: [regime()]
#' writes its levels as a first level and positive gaps, so any draw on the
#' unconstrained scale is ordered by construction, and a persistence on a
#' partial-autocorrelation chart is stationary at any coordinate.
#'
#' # The width
#'
#' The standard deviation of the draw is `sd / 2`. At the full width of `sd`
#' a score-driven filter with a Poisson response can produce non-finite
#' series, and a persistence on the partial-autocorrelation chart covers
#' nearly the whole stationary region, which describes noise and not a random
#' effect.
#'
#' # Coordinates covered by a penalty
#'
#' These coordinates are drawn as well, and a caller that has access to the
#' penalty is expected to overwrite them. The deviation of a development
#' follows the prior that its sub-term declares (a Gaussian random effect or a
#' Laplace one). This method does not receive that prior, because its scale is
#' a hyperparameter that the caller chooses. [term_penalties()] identifies
#' those coordinates, and \pkg{penalties7}'s `penalty_draw()` draws them.
#'
#' @param term A built structural term.
#' @param sd The width of the draw, `1` by default. The standard deviation
#'   applied is `sd / 2`, for the reason given in the details.
#' @param ... Passed to methods. The shipped methods do not use these
#'   arguments.
#'
#' @return A named numeric vector on the unconstrained scale, of length
#'   [term_npar()] and named as [term_params()], exactly as [term_start()]
#'   returns.
#'
#' @seealso [term_start()] for the point from which the draw spreads,
#'   [term_links()] for the scale on which the draw is expressed,
#'   [term_simulate()] for the response drawn once these are chosen.
#'
#' @examples
#' set.seed(1)
#' g <- gas(p = 1, q = 1)
#'
#' # The draw has the same names and the same scale as term_start().
#' identical(names(term_draw(g)), term_params(g))
#'
#' # Read on the parameters' own charts, a loading stays positive and a
#' # persistence stays stationary, whatever comes out.
#' lk <- term_links(g)
#' z <- term_draw(g)
#' c(kappa1 = linkfunctions7::linkinv(lk$kappa1, z[["kappa1"]]),
#'   pacf1 = linkfunctions7::linkinv(lk$pacf1, z[["pacf1"]]))
#'
#' # A larger sd widens the draws, and with them the range of models.
#' stats::sd(replicate(200, term_draw(g, sd = 2)[["omega"]]))
#'
#' @export
#' @aliases term_draw.structural_term
term_draw <- S7::new_generic("term_draw", "term",
  function(term, sd = 1, ...) S7::S7_dispatch())

S7::method(term_draw, structural_term) <- function(term, sd = 1, ...) {
  if (!is.numeric(sd) || length(sd) != 1L || !is.finite(sd) || sd < 0) {
    stop("'sd' must be one finite non-negative number.", call. = FALSE)
  }
  z <- term_start(term)
  # HALF the width, measured: at the full one a Poisson filter leaves the
  # doubles and a persistence drawn on a rhobit chart covers the whole
  # stationary region, which is noise rather than a random effect.
  z + stats::rnorm(length(z), 0, sd / 2)
}

#' @title Where a Fit of a Term's Own Parameters Starts, Several Times
#'
#' @description
#' Returns the starting points that a fitting layer tries for a structural
#' term's own parameters, as a list of vectors of the kind that
#' [term_start()] returns. The layer fits from each and keeps the best.
#'
#' @details
#' The base method returns the single start that [term_start()] gives, so a
#' term without a method of its own is fitted once. [regime()] returns
#' `n_start` starts, because its likelihood has several maxima and the starts
#' cover different regions of its parameters. In each start after the first,
#' the free entries of the transition matrix are drawn from a normal
#' distribution of standard deviation 2 and normal draws of standard deviation
#' 0.7 are added to the log gaps. The seed of the start is fixed, and the state
#' of the random number generator of the caller is restored.
#'
#' The first element is always [term_start()]'s, so a fit with one start is
#' the fit without this generic.
#'
#' @param term A built structural term.
#' @param ... Passed to [term_start()].
#' @param target The response on the scale of the predictor, or `NULL`, as
#'   [term_start()] reads it.
#'
#' @return A non-empty list of named numeric vectors on the unconstrained
#'   scale, each of length [term_npar()] and named as [term_params()].
#'
#' @seealso [term_start()] for the first start, [regime()] for a term that
#'   returns several.
#'
#' @examples
#' length(term_starts(gas(p = 1, q = 1)))
#' st <- term_starts(regime(k = 2, n_start = 3))
#' length(st)
#' identical(st[[1]], term_start(regime(k = 2)))
#'
#' @export
#' @aliases term_starts.structural_term
term_starts <- S7::new_generic("term_starts", "term",
  function(term, ..., target = NULL) S7::S7_dispatch())

S7::method(term_starts, structural_term) <- function(term, ...,
                                                     target = NULL) {
  list(term_start(term, ..., target = target))
}


#' @title Apply a Structural Term to a Linear Predictor
#'
#' @description
#' Runs the term's recursion over the data and returns the predictor that it
#' produces, together with the derivative of that predictor with respect
#' to the term's own parameters. This operation defines a
#' filter-shaped structural term such as [gas()]: the predictor at one
#' observation depends on the others, so it cannot be written as a block of
#' columns. The structural terms of the likelihood shape, [regime()] and the
#' marginal break-point term, have no filter and implement [term_loglik()]
#' instead.
#'
#' @details
#' Writing \eqn{\eta_t^{0}} for the static predictor supplied in `eta`
#' and \eqn{\psi} for the term's parameters, the filter returns the pair
#'
#' \deqn{\eta_t = \eta_t^{0} + f_t(\psi),
#'   \qquad J_{tj} = \frac{\partial \eta_t}{\partial \psi_j},}
#'
#' where \eqn{f_t} is the term's own recursion, driven by
#' `score` and `curvature` evaluated at the predictor already
#' produced. Both are read at \eqn{\eta_t}, so `curvature` is the second
#' derivative \eqn{\partial^{2} \ell_t / \partial \eta^{2}} and is negative at
#' an ordinary observation. The sign of `curvature` matters: the predictor
#' does not depend on `curvature`, but the Jacobian does, so passing the
#' information, which is the negative of the curvature, returns the correct
#' predictor together with a Jacobian that is not its derivative, and the
#' filter does not signal an error.
#'
#' The filter applies no scaling of its own. The callbacks supplied are the
#' driving sequence of the recursion, and a fitting layer that uses the
#' `scaling` of [gas()] builds callbacks that return the scaled score and its
#' derivative.
#'
#' The derivative is returned because the recursion is the only place in which
#' it can be computed. Propagating it beside the state costs one extra vector
#' per parameter and is exact; a model layer differencing the filter instead
#' would pay one pass per parameter and inherit the error of the difference.
#'
#' @param term A built structural term.
#' @param eta The static part of the linear predictor, one value per
#'   observation.
#' @param y The response.
#' @param score A function `function(e, i)` of a predictor value `e` and an
#'   observation index `i`, returning the single value that drives the
#'   recursion at that observation: the derivative of the log-likelihood with
#'   respect to the predictor, multiplied by the information scaling when the
#'   `scaling` of [gas()] is not zero.
#' @param curvature A function `function(e, i)` returning the derivative of
#'   that driving value with respect to `e`, which is the second derivative of
#'   the log-likelihood when the scaling is zero.
#' @param psi The term's parameters, on the parameter scale, named as
#'   [term_params()].
#' @param ... Passed to methods.
#'
#' @return A list with `eta`, the predictor that the term produces,
#'   `jacobian`, an `n` by `length(psi)` matrix of its
#'   derivatives with respect to `psi`, and `curv`, the value of
#'   `curvature` at each predictor. That last one is returned because the
#'   recursion evaluates it anyway, so a consumer running a second pass at
#'   the same point, as [term_adjoint()] does, reads it instead of
#'   evaluating the callback again.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:20, y = rnorm(20))
#' term <- term_build(gas(p = 1, q = 1, time = t), dd)
#'
#' # The callbacks return the score and the curvature that a Gaussian mean
#' # would supply.
#' out <- term_filter(term, eta = rep(0, 20), y = dd$y,
#'                    score = function(e, i) dd$y[i] - e,
#'                    curvature = function(e, i) -1,
#'                    psi = list(omega = 0.1, kappa1 = 0.3, pacf1 = 0.5))
#' head(out$eta, 3)
#' dim(out$jacobian)
#'
#' @seealso [gas()] for the filter shape and [regime()] for the likelihood
#'   shape,
#'   [term_static_deriv()] for the same recursion in the static predictor,
#'   [term_adjoint()] for the reverse pass, [term_continue()] for the
#'   recursion past the series.
#' @export
#' @aliases term_filter.structural_term
term_filter <- S7::new_generic("term_filter", "term",
  function(term, eta, y, score, curvature, psi, ...) S7::S7_dispatch())

#' @title The Derivative of a Filtered Predictor in the Static One
#'
#' @description
#' Returns how the predictor that a structural term produces moves when the
#' static part of the predictor moves, one row per observation and one column
#' per direction that the caller supplies.
#'
#' @details
#' A score-driven term's level is driven by scores read at the predictor that
#' the recursion is producing, so a coefficient in the same equation reaches
#' the level as well as the static part: writing \eqn{f_t} for the level and
#' \eqn{x_t} for a row of the design,
#' \deqn{\frac{\partial \eta_t}{\partial \beta} = x_t +
#'   \frac{\partial f_t}{\partial \beta},}
#' and the second piece obeys the recursion that the filter already runs,
#' \deqn{\frac{\partial f_t}{\partial \beta} =
#'   \sum_i \kappa_i \, \ell''_{t-i}
#'     \frac{\partial \eta_{t-i}}{\partial \beta} +
#'   \sum_j \xi_j \frac{\partial f_{t-j}}{\partial \beta},}
#' with \eqn{\kappa_i} the score loadings, \eqn{\xi_j} the autoregressive
#' coefficients and \eqn{\ell''_t} the curvature. The curvature is the one
#' that [term_filter()] returns, so this pass does not evaluate any callback
#' and consists only of arithmetic.
#'
#' Without this derivative, the standard error of the predictor accounts for
#' the static part only and leaves out the part that runs through the filter.
#'
#' The base method returns `NULL`. A term that is not a filter carries no
#' state, so the derivative is the design row itself and the caller needs
#' nothing from the term.
#'
#' @param term A built structural term.
#' @param curv The curvature at each predictor, as `term_filter`
#'   returns it.
#' @param X The directions to propagate, one column each. Ordinarily `X` is
#'   the equation's design, so the result has one row per observation and one
#'   column per coefficient.
#' @param psi The term's parameters, on the parameter scale, named as
#'   [term_params()].
#' @param ... Passed to methods.
#'
#' @return A matrix of `X`'s dimensions, holding the full derivative
#'   \eqn{\partial \eta_t/\partial \beta = x_t + \partial f_t/\partial \beta}
#'   for each column of `X`, so that the design is already included. The
#'   value is `NULL` for a term that carries no state.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:20, y = rnorm(20), x = rnorm(20))
#' term <- term_build(gas(p = 1, q = 1, time = t), dd)
#' psi <- list(omega = 0.1, kappa1 = 0.3, pacf1 = 0.5)
#' out <- term_filter(term, eta = rep(0, 20), y = dd$y,
#'                    score = function(e, i) dd$y[i] - e,
#'                    curvature = function(e, i) -1, psi = psi)
#' D <- term_static_deriv(term, out$curv, cbind(1, dd$x), psi)
#' dim(D)
#'
#' @seealso [term_filter()], [term_adjoint()]
#' @export
term_static_deriv <- S7::new_generic("term_static_deriv", "term",
  function(term, curv, X, psi, ...) S7::S7_dispatch())

#' @title Continuing a Structural Term Past the Observed Series
#'
#' @description
#' Returns the contribution that a structural term makes at rows that come
#' after the rows on which it was built, continuing its recursion rather than
#' restarting it.
#'
#' @details
#' A structural term's contribution at one observation is not a function of
#' that observation: it is the state that a recursion has reached, so
#' predicting past the series means carrying the state forward. This is
#' possible without simulation because the quantity driving the recursion has
#' zero conditional mean, the score itself for a score-driven term, so beyond
#' the data the recursion is deterministic.
#'
#' Among the shipped classes, [gas()] implements the generic. Its method needs
#' `time`, rejects rows that fall inside the observed series and rows of a
#' group that the sample does not contain, and counts a gap between two times
#' as one step. A lag that falls before the observed series takes the starting
#' level of the recursion.
#'
#' The base method signals an error. A term with state whose next value is
#' undefined cannot give a prediction past the series, and returning zero
#' there would read as a term with no effect at all.
#'
#' @param term A built structural term.
#' @param psi The term's parameters, named as [term_params()].
#' @param f_past The term's contribution at each observed row.
#' @param s_past The driving quantity at each observed row.
#' @param newdata The rows onto which the recursion continues.
#' @param ... Passed to methods.
#'
#' @return A numeric vector of `nrow(newdata)` contributions.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:20, y = rnorm(20))
#' term <- term_build(gas(p = 1, q = 1, time = t), dd)
#' psi <- list(omega = 0.1, kappa1 = 0.3, pacf1 = 0.5)
#' out <- term_filter(term, eta = rep(0, 20), y = dd$y,
#'                    score = function(e, i) dd$y[i] - e,
#'                    curvature = function(e, i) -1, psi = psi)
#' sc <- dd$y - out$eta
#' term_continue(term, psi, out$eta, sc, data.frame(t = 21:23))
#'
#' @seealso [term_filter()]
#' @export
term_continue <- S7::new_generic("term_continue", "term",
  function(term, psi, f_past, s_past, newdata, ...) S7::S7_dispatch())

#' @title Drawing a Response From a Structural Term
#'
#' @description
#' Returns the predictor that a structural term produces when the response is
#' being generated rather than read, together with whatever latent quantity
#' the term drew on the way.
#'
#' @details
#' Simulating from a model that carries state differs from fitting one in the
#' direction in which the response moves. A term whose contribution does not
#' read the response can report that contribution and leave the drawing to the
#' caller: the levels of a latent chain ([regime()]) and the break-point of a
#' group drawn from its prior (a marginal break-point term) are both of this
#' kind. A score-driven term cannot do so. Its level at one time is driven by
#' the score of the response at the time before, so the response has to be
#' drawn as the recursion runs.
#'
#' The same calling convention covers both. The caller supplies `draw`, a
#' function of a predictor and a row index returning one response value, and
#' the method returns the predictor that it produced. A method that drew returns
#' the responses as well, and a method that did not returns `NULL` there,
#' leaving the caller to draw at the predictor. [regime()] draws one chain per
#' group, in time order and started from the stationary distribution, and
#' returns the state of each row as `latent`.
#'
#' @param term A built structural term.
#' @param psi The term's parameters, named as [term_params()].
#' @param eta The static part of the predictor, one value per observation.
#' @param draw A function `(e, i)` returning one response value drawn at
#'   predictor `e` for observation `i`.
#' @param ... Passed to methods.
#'
#' @return A list with `eta`, the predictor; `y`, the responses
#'   drawn or `NULL`; and `latent`, whatever the term drew.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:30)
#' term <- term_build(gas(p = 1, q = 1, time = t), dd)
#' out <- term_simulate(term, list(omega = 0.5, kappa1 = 0.3, pacf1 = 0.6),
#'                      rep(0, 30),
#'                      draw = function(e, i) stats::rnorm(1, e, 1))
#' head(out$y, 3)
#'
#' @seealso [term_filter()], [term_continue()]
#' @export
term_simulate <- S7::new_generic("term_simulate", "term",
  function(term, psi, eta, draw, ...) S7::S7_dispatch())

#' @name term_simulate.model_term
#' @title Base Method of term_simulate()
#' @description
#' The base method signals an error: an ordinary term contributes its block
#' times its coefficients and has nothing of its own to draw.
#' @param term A term.
#' @param psi,eta,draw Ignored.
#' @param ... Ignored.
#' @return Nothing; the method always signals an error.
#' @seealso [term_simulate()]
#' @keywords internal
S7::method(term_simulate, model_term) <- function(term, psi, eta, draw, ...) {
  stop(sprintf("'%s' does not implement term_simulate().",
               class(term)[[1L]]), call. = FALSE)
}

#' @name term_continue.model_term
#' @title Base Method of term_continue()
#' @description
#' The base method signals an error.
#' @param term A term.
#' @param psi,f_past,s_past,newdata Ignored.
#' @param ... Ignored.
#' @return Nothing; the method always signals an error.
#' @seealso [term_continue()]
#' @keywords internal
S7::method(term_continue, model_term) <- function(term, psi, f_past, s_past,
                                                  newdata, ...) {
  stop(sprintf("'%s' defines no rule by which its contribution continues past the series.",
               class(term)[[1L]]), call. = FALSE)
}

#' @name term_static_deriv.model_term
#' @title Base Method of term_static_deriv()
#' @description
#' The base method returns `NULL`: an ordinary term's contribution at
#' one observation reads only that observation.
#' @param term A term.
#' @param curv,X,psi Ignored.
#' @param ... Ignored.
#' @return `NULL`.
#' @seealso [term_static_deriv()]
#' @keywords internal
S7::method(term_static_deriv, model_term) <- function(term, curv, X, psi,
                                                      ...) {
  NULL
}

#' @title Which of a Term's Parameters Acts as an Intercept
#'
#' @description
#' Returns the name of the parameter, if any, that shifts the term's
#' contribution by a constant and is therefore the same direction as an
#' intercept in the equation in which the term sits.
#'
#' @details
#' It exists so that a fitting layer can resolve the confounding instead of
#' rejecting the model. A score-driven level \eqn{\omega} and a regime's first
#' level both add a constant to their equation's predictor: with an
#' intercept there too, adding \eqn{c} to one and subtracting the matching
#' amount from the other leaves every predictor unchanged, and the
#' likelihood is flat along that direction. Which of the two is dropped is
#' determined in the layer, because only the layer has the full list of terms
#' in the equation.
#'
#' The base method returns `character(0)`, which means that the term does not
#' shift its contribution by a constant and the question does not arise. The
#' method for [regime()] returns `"level1"`. The method for [gas()] returns
#' `"omega"` when the level is not developed and otherwise every coordinate of
#' the development whose design column is constant, which can be more than
#' one.
#'
#' @param term A term.
#' @param ... Passed to methods.
#'
#' @return A character vector: `character(0)` for a term that shifts nothing
#'   by a constant, `"level1"` for a regime, `"omega"` for an undeveloped
#'   [gas()] level, and otherwise the names of the level's coordinates whose
#'   design column is constant.
#'
#' @examples
#' term_level_param(gas(p = 1, q = 1))
#' term_level_param(regime(k = 2))
#' term_level_param(linpar(~x))
#'
#' @seealso [term_level_design()] for the same question with the level
#'   developed, [term_params()] for the names from which the answer is taken.
#' @export
#' @aliases term_level_param.model_term
term_level_param <- S7::new_generic("term_level_param", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_level_param, model_term) <- function(term, ...) character(0)

#' @title The Design of a Term's Level Development
#'
#' @description
#' Returns, where the parameter that [term_level_param()] names is developed
#' with covariates, the design of that development, with one column per
#' coordinate named as [term_params()] names it. The value is `NULL`
#' for a scalar level and for every other term.
#'
#' @details
#' It exists for the subspace form of the confounding question.
#' `term_level_param` covers the case of a constant: a coordinate whose
#' column is constant shifts the equation's predictor exactly as an
#' intercept does. If the level is developed, a direction of the
#' development's span that also lies in the span of the equation's design
#' raises the same question for that direction. The question can be settled
#' only in a fitting layer, which holds both designs, and this generic
#' supplies the design of the development, which the layer cannot compute
#' from the term alone.
#'
#' @param term A built term.
#' @param ... Passed to methods.
#'
#' @return A numeric matrix with named columns, or `NULL`.
#'
#' @examples
#' is.null(term_level_design(linpar(~x)))
#'
#' @seealso [term_level_param()] for the scalar case, [term_params()] for
#'   the names that its columns carry.
#' @export
#' @aliases term_level_design.model_term
term_level_design <- S7::new_generic("term_level_design", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_level_design, model_term) <- function(term, ...) NULL


#' @title The Quantities That a Fitted Term Reports
#'
#' @description
#' Returns the quantities that a reader reads, with the Jacobian from the
#' term's own parameters on the unconstrained scale, so that a caller who
#' holds their variance matrix can carry it across by the delta method.
#'
#' @details
#' A term's parameters are the coordinates on which it is estimated, chosen so
#' that a search runs unconstrained, and they are not always the quantities
#' of interest in the model. The clearest case is a score-driven persistence,
#' which rides a partial autocorrelation because the stationary region is
#' not a box: the autoregressive coefficient \eqn{\xi_1} is a function of the
#' whole chart and coincides with the coordinate only at \eqn{q = 1}.
#' Reporting the coordinate under the coefficient's name would label one
#' quantity with the name of another. [gas()] therefore reports `xi1`, ..., and
#' [regime()] reports the transition probabilities `p1.1`, `p1.2`, ... row by
#' row.
#'
#' Each row gives a value and the row of \eqn{\partial(\text{value}) /
#' \partial\zeta} at the current parameters, so a standard error is
#' \eqn{\sqrt{J V J^\top}} and an interval is built on whichever scale keeps
#' the quantity in its own set, exactly as \pkg{parameters7}'s
#' `param_readable()` does for a matrix parameter.
#'
#' The base method reports the parameters themselves on the parameter scale,
#' with the diagonal Jacobian of their links. That is appropriate for every
#' term whose coordinates are already its quantities.
#'
#' @param term A built term.
#' @param zeta The term's parameters on the unconstrained scale, named as
#'   [term_params()].
#' @param ... Passed to methods.
#'
#' @return A list with `name`, `value`, `jacobian` (one row
#'   per quantity and one column per parameter) and `scale`, the link
#'   on which the interval for each quantity is built.
#'
#' @examples
#' term_readable(gas(p = 1, q = 1), c(omega = 0.3, kappa1 = 0.4, pacf1 = 0.8))
#'
#' @seealso [term_params()] for the coordinates, [term_links()] for their
#'   charts, [parameters7::param_readable()] for the same shape applied to a
#'   matrix parameter.
#' @export
#' @aliases term_readable.model_term
term_readable <- S7::new_generic("term_readable", "term",
  function(term, zeta, ...) S7::S7_dispatch())

S7::method(term_readable, model_term) <- function(term, zeta, ...) {
  nm <- term_params(term)
  links <- term_links(term)
  z <- unlist(zeta[nm])
  val <- vapply(nm, function(j)
    linkfunctions7::linkinv(links[[j]], z[[j]]), numeric(1))
  J <- diag(vapply(nm, function(j)
    linkfunctions7::dlinkinv(links[[j]], z[[j]]), numeric(1)),
    nrow = length(nm))
  dimnames(J) <- list(nm, nm)
  list(name = nm, value = unname(val), jacobian = J, scale = links[nm])
}


#' @title Differentiate a Structural Term Backwards
#'
#' @description
#' Returns the derivative of a caller's objective with respect to the static
#' predictor that the term was handed, and with respect to the sequence of
#' scores that it was given. Both derivatives account for the recursion.
#'
#' @details
#' [term_filter()] returns the derivative of the predictor in the
#' term's own parameters, which is the piece that is needed to estimate those
#' parameters. Estimating the coefficients of the same equation needs a
#' different derivative: the level at one time is driven by the scores at
#' earlier times, which are read at predictors that those coefficients also
#' enter, so the derivative of the predictor in a coefficient carries a term
#' that the block does not. Propagating that forward would cost one derivative
#' array per coefficient; the reverse recursion here costs one pass whatever
#' their number, and is exact.
#'
#' Two derivatives are returned because the score that the caller supplies
#' depends on more than the predictor at which it is read. Writing
#' \eqn{s_t} for that score and \eqn{\bar{s}_t} for
#' `dscore`, the derivative of the objective in any quantity on which the
#' score depends is
#'
#' \deqn{\frac{\partial L}{\partial \theta}
#'   = \left.\frac{\partial L}{\partial \theta}\right|_{\mathrm{direct}}
#'     + \sum_t \bar{s}_t \frac{\partial s_t}{\partial \theta},}
#'
#' so a model layer whose score is the derivative of its log-likelihood in
#' one distribution parameter obtains the derivative in the predictor of
#' another by multiplying `dscore` by the mixed second derivative of
#' that log-likelihood. `deta` is that formula applied to the term's own
#' equation, where the second factor is the curvature.
#'
#' @param term A built structural term.
#' @param eta The static part of the predictor, as
#'   [term_filter()] takes it.
#' @param y The response.
#' @param score,curvature The callbacks of [term_filter()].
#' @param psi The term's parameters, named as [term_params()].
#' @param g The direct derivative of the objective in the predictor that the
#'   term produced, one value per observation.
#' @param ... Passed to methods.
#'
#' @return A list of two numeric vectors, each of one value per observation:
#'   `deta`, the derivative of the objective in the static predictor the term
#'   was handed, and `dscore`, its derivative in the sequence of scores. A
#'   caller multiplies `dscore` by a mixed second derivative of its
#'   log-density to obtain the derivative in another equation's predictor.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:20, y = rnorm(20))
#' term <- term_build(gas(p = 1, q = 1, time = t), dd)
#' out <- term_adjoint(term, eta = rep(0, 20), y = dd$y,
#'                     score = function(e, i) dd$y[i] - e,
#'                     curvature = function(e, i) -1,
#'                     psi = list(omega = 0.1, kappa1 = 0.3, pacf1 = 0.5),
#'                     g = rep(1, 20))
#' head(out$deta, 3)
#'
#' @seealso [term_filter()] for the forward pass and the callbacks,
#'   [term_curvature()] for the second order, [term_static_deriv()] for the
#'   forward derivative in the static predictor.
#' @export
#' @aliases term_adjoint.structural_term
term_adjoint <- S7::new_generic("term_adjoint", "term",
  function(term, eta, y, score, curvature, psi, g, ...) S7::S7_dispatch())

S7::method(term_adjoint, structural_term) <- function(term, eta, y, score,
                                                      curvature, psi, g, ...) {
  stop(sprintf("the term class '%s' does not implement term_adjoint().",
               attr(S7::S7_class(term), "name")), call. = FALSE)
}


#' @title Second Derivatives of a Structural Term's Predictor
#'
#' @description
#' Returns the Jacobian of the predictor that the term produces with respect
#' to a caller's unknowns, and the second derivative of that predictor
#' contracted against a caller's weights. This second derivative is the piece of
#' the Hessian of the log-likelihood that [term_adjoint()] does not give.
#'
#' @details
#' The gradient of a model carrying a filter needs no second derivative and
#' no forward Jacobian: the reverse recursion gives it in one pass. The
#' curvature needs both. Writing \eqn{u} for the caller's unknowns, \eqn{p}
#' for the distribution parameter whose predictor carries the filter,
#' \eqn{D_t = \partial e_t/\partial u} for the derivative of that predictor,
#' \eqn{V_{q,t}} for the derivative of the predictor of parameter \eqn{q}
#' (with \eqn{V_{p,t} = D_t}) and \eqn{E_t = \partial^2
#' e_t/\partial u \partial u^\top}, the Hessian of the log-likelihood in
#' \eqn{u} is
#'
#' \deqn{\sum_t w_t \sum_{q,r} \ell_{qr,t} V_{q,t}^\top V_{r,t}
#'   + \sum_t w_t \ell_{p,t} E_t,}
#'
#' and the observed information is its negative. The first sum is a weighted
#' crossproduct that the caller assembles. The second sum is the value that
#' this method returns, for the weights \eqn{g_t = w_t\ell_{p,t}} that the
#' caller supplies.
#'
#' `seed` is \eqn{\partial \eta^{0}/\partial u}, one row per
#' observation: the caller supplies how its unknowns reach the static
#' predictor, and the term holds no other information about them. The model's
#' own derivatives enter through `blocks`, since the score that drives the
#' recursion depends on every equation, beyond the predictor at which it is
#' read. The function `blocks` is called once per observation, with the
#' predictor there, the
#' index, the Jacobian row that the recursion has reached (restricted to the
#' active coordinates of the group) and the indices of those coordinates, and
#' it returns the two quantities that seed the first and second derivatives of
#' that score,
#'
#' \deqn{\dot S_t = \ell_{pp,t}D_t + \texttt{cross}, \qquad
#'   \ddot S_t = \ell_{pp,t}E_t + \texttt{M},}
#'
#' with `cross` \eqn{= \sum_{q\ne p}\ell_{pq,t}V_{q,t}} and `M`
#' \eqn{= \sum_{r,r'}\ell_{prr',t}V_{r,t}^\top V_{r',t}}, the third
#' derivatives of the log-density in the predictors. A model of one equation
#' has `cross` zero and `M` equal to \eqn{\ell_{ppp}D^\top D}.
#'
#' @param term A built structural term.
#' @param eta The static part of the predictor.
#' @param y The response.
#' @param score,curvature The callbacks of [term_filter()].
#' @param psi The term's parameters, named as [term_params()].
#' @param g The weights against which the second derivative is contracted,
#'   one per observation.
#' @param seed The derivative of the static predictor in the caller's
#'   unknowns, one row per observation.
#' @param blocks A function `function(e, i, D, act)` of the predictor, the
#'   observation index, the current Jacobian row restricted to the active set
#'   and that set, returning `cross` and `M` on the active set. A function of
#'   the first three arguments only is accepted and is given the full row.
#' @param ... Passed to methods.
#'
#' @return A list of two: `jacobian`, an `n` by `ncol(seed)` matrix holding
#'   the derivative of the predictor in the caller's unknowns, and
#'   `curvature`, the `ncol(seed)` by `ncol(seed)` second derivative
#'   contracted against `g`. Both are exact.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:20, y = rnorm(20))
#' term <- term_build(gas(p = 1, q = 1, time = t), dd)
#'
#' # The model is a gaussian mean of unit variance with one equation, so it
#' # contributes no cross term and no third derivative.
#' m <- 3L
#' out <- term_curvature(
#'   term, eta = rep(0, 20), y = dd$y,
#'   score = function(e, i) dd$y[i] - e,
#'   curvature = function(e, i) -1,
#'   psi = list(omega = 0.1, kappa1 = 0.3, pacf1 = 0.5),
#'   g = rep(1, 20), seed = matrix(0, 20, m),
#'   blocks = function(e, i, D) list(cross = numeric(m),
#'                                   M = matrix(0, m, m)))
#' dim(out$jacobian)
#' dim(out$curvature)
#'
#' @seealso [term_adjoint()] for the first order, [term_third()] for the
#'   third, [term_filter()] for the callbacks.
#' @export
#' @aliases term_curvature.structural_term
term_curvature <- S7::new_generic("term_curvature", "term",
  function(term, eta, y, score, curvature, psi, g, seed, blocks, ...)
    S7::S7_dispatch())

S7::method(term_curvature, structural_term) <- function(term, eta, y, score,
                                                        curvature, psi, g,
                                                        seed, blocks, ...) {
  stop(sprintf("the term class '%s' does not implement term_curvature().",
               attr(S7::S7_class(term), "name")), call. = FALSE)
}

#' @title Third Derivatives of a Structural Term's Predictor, in One Direction
#'
#' @description
#' Returns the second derivative of [term_curvature()] differentiated once
#' more along a single direction, and the derivative of the term's Jacobian
#' along that same direction. The exact gradient of a marginal criterion
#' needs this when a penalty covers the term's own parameters.
#'
#' @details
#' A marginal criterion carries \eqn{-\frac{1}{2}\log|K|} at the penalized
#' mode, so its gradient in a hyperparameter needs \eqn{\mathrm{tr}(M\,\partial
#' K/\partial u\,[v])}, with \eqn{v} the direction in which the mode moves.
#' Where the predictor is \eqn{X\beta} that contraction is the family's third
#' derivative against the design and nothing else. Where a filter produces the
#' predictor, \eqn{K} carries \eqn{\sum_t w_t \ell_{p,t}E_t} as well, and
#' differentiating it requires \eqn{\partial^3 e_t/\partial u^3}.
#'
#' The derivative is only ever contracted. The full third derivative is an
#' \eqn{m^3} array per observation and is never formed: the criterion needs it
#' only along the one direction in which the mode moves, so the method
#' propagates a matrix per observation, which has the same size as the
#' curvature and therefore the same \eqn{O(nm^2)}. A caller that needs several
#' hyperparameters calls
#' this once per direction, which is cheaper than one \eqn{O(nm^3)} assembly
#' whenever the hyperparameters are fewer than the unknowns.
#'
#' Each order of differentiation of the predictor involves one more order of
#' the family, because the recursion is driven by a score read at the
#' predictor that it produces: the curvature's `M` is built from third
#' derivatives, and this method needs a fourth. `blocks` therefore returns two
#' quantities beyond the curvature's, at an observation, from the current
#' Jacobian row \eqn{D} and the direction \eqn{v}:
#'
#' \deqn{\texttt{dcurv} = \sum_r \ell_{ppr}V_r, \qquad
#'   \texttt{N} = \sum_{r,r'}\Big(\sum_{r''}\ell_{prr'r''}(V_{r''}\cdot v)\Big)
#'     V_r^\top V_{r'},}
#'
#' where `dcurv` serves for both the derivative of the curvature along
#' \eqn{v}, which is \eqn{\texttt{dcurv}\cdot v}, and the two terms that
#' arise from differentiating \eqn{M}'s own \eqn{V_p}.
#'
#' The base method returns zeros, which is correct for a term whose predictor
#' is a block of columns times coefficients, since the second derivative of
#' such a predictor in the coefficients is zero. A term that bends the
#' predictor has to implement the method itself, because the base method would
#' report zeros for it. A structural term signals an error instead: it has a
#' second derivative, so a zero third derivative would misstate it.
#'
#' @param term A built term.
#' @param eta The static part of the predictor.
#' @param y The response.
#' @param score,curvature The callbacks of [term_filter()].
#' @param psi The term's parameters, named as [term_params()].
#' @param g The weights against which the third derivative is contracted,
#'   one per observation.
#' @param seed The derivative of the static predictor in the caller's
#'   unknowns.
#' @param blocks A function of the predictor, the index, the current Jacobian
#'   row and the active set, returning `cross`, `M`, `dcurv`
#'   and `N`.
#' @param direction The direction \eqn{v}, in the caller's unknowns.
#' @param ... Passed to methods.
#'
#' @return A list with `jacobian`, the derivative of the predictor in the
#'   caller's unknowns; `dphi`, the second derivative contracted against
#'   `direction`, one row per observation; and `curvature`, the
#'   third derivative contracted against both `g` and `direction`.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:20, y = rnorm(20))
#' term <- term_build(linpar(~t), dd)
#' # an additive term bends no predictor, so every order above the first is
#' # zero, as the base method returns
#' out <- term_third(term, rep(0, 20), dd$y,
#'                   score = function(e, i) dd$y[i] - e,
#'                   curvature = function(e, i) -1, psi = list(),
#'                   g = rep(1, 20), seed = matrix(0, 20, 2),
#'                   blocks = function(e, i, D, act) NULL,
#'                   direction = c(1, 0))
#' all(out$curvature == 0)
#'
#' @seealso [term_curvature()] for the second order this differentiates,
#'   [term_adjoint()] for the first, [statmodels7::reml()] for the criterion
#'   whose gradient requires it.
#' @export
#' @aliases term_third.model_term term_third.structural_term
term_third <- S7::new_generic("term_third", "term",
  function(term, eta, y, score, curvature, psi, g, seed, blocks, direction,
           ...) S7::S7_dispatch())

S7::method(term_third, model_term) <- function(term, eta, y, score, curvature,
                                               psi, g, seed, blocks,
                                               direction, ...) {
  seed <- as.matrix(seed)
  m <- ncol(seed)
  list(jacobian = seed, dphi = matrix(0, nrow(seed), m),
       curvature = matrix(0, m, m))
}

# A structural term REFUSES rather than inheriting the zero above. The base
# method is right for a term whose predictor is a block of columns, whose
# second derivative is already zero; a term that bends the predictor and has
# not written its third derivative would otherwise report zero, which a
# caller cannot tell from a term that genuinely has none.
S7::method(term_third, structural_term) <- function(term, eta, y, score,
                                                    curvature, psi, g, seed,
                                                    blocks, direction, ...) {
  stop(sprintf("the term class '%s' does not implement term_third().",
               attr(S7::S7_class(term), "name")), call. = FALSE)
}

#' @title Fourth Derivatives of a Structural Term's Predictor, in Two Directions
#'
#' @description
#' Returns [term_third()] differentiated once more along a second direction.
#' The exact Hessian of a marginal criterion needs this when the model
#' carries a term that bends the predictor.
#'
#' @details
#' A marginal criterion carries \eqn{-\frac{1}{2}\log|K|} at the penalized
#' mode, so its gradient in one hyperparameter reads \eqn{\partial K/\partial
#' u} along the direction in which the mode moves, which is [term_third()], and
#' its second derivative in a pair of hyperparameters reads the same object
#' along two directions. Where the predictor is \eqn{X\beta} that is the
#' family's derivatives against the design and nothing else; where a filter
#' produces it, \eqn{K} carries \eqn{\sum_t w_t \ell_{p,t}E_t} as well, and
#' differentiating that twice requires \eqn{\partial^4 e_t/\partial u^4}.
#'
#' The derivative is contracted in two directions only. The full fourth
#' derivative is an \eqn{m^4} array per observation and is never formed: the
#' method propagates a matrix per observation, which has the same size as the
#' curvature and therefore the same \eqn{O(nm^2)}. A caller that needs a
#' Hessian over \eqn{n_h} hyperparameters calls this once per pair.
#'
#' The recursion carries five states, \eqn{F_t}, \eqn{\Phi_t},
#' \eqn{\Psi^{(v)}_t}, \eqn{\Psi^{(w)}_t} and \eqn{\Omega_t}. Both third
#' derivatives are carried, because the product rule at the fourth order pairs
#' each direction's third derivative with the other direction, so a recursion
#' holding one contraction cannot reach this order.
#'
#' Each order of differentiation involves one more order of the family, so
#' beyond the third order's `dcurv` and `N` the `blocks` callback returns, at
#' an observation,
#'
#' \deqn{\texttt{Q} = \sum_{r,r'}\ell_{pprr'}V_r^\top V_{r'}, \qquad
#'   \texttt{P} = \sum_{r,r'}\Big(\sum_{s,s'}\ell_{prr'ss'}
#'     (V_s\cdot v)(V_{s'}\cdot w)\Big)V_r^\top V_{r'},}
#'
#' with `cppp` \eqn{= \ell_{ppp}} beside them, and `N` a list of
#' two, one per direction. The family's fifth derivative enters only through
#' `P`. The scalar `cppp` is not recovered by any contraction;
#' it is the amount by which the level's own curvature moves.
#'
#' The base method returns zeros, so an additive term is covered and its
#' author writes no method. A structural term signals an error, for the reason
#' that [term_third()] gives: a zero fourth derivative would misstate it.
#'
#' @param term A built term.
#' @param eta The static part of the predictor.
#' @param y The response.
#' @param score,curvature The callbacks of [term_filter()].
#' @param psi The term's parameters, named as [term_params()].
#' @param g The weights against which the fourth derivative is contracted,
#'   one per observation.
#' @param seed The derivative of the static predictor in the caller's
#'   unknowns.
#' @param blocks A function of the predictor, the index, the current Jacobian
#'   row and the active set, returning `cross`, `M`,
#'   `dcurv`, `N`, `Q`, `P` and `cppp`.
#' @param directions A list of two numeric vectors, the directions \eqn{v}
#'   and \eqn{w} in the caller's unknowns. The result is symmetric in the
#'   two.
#' @param ... Passed to methods.
#'
#' @return A list with `jacobian`, the derivative of the predictor in
#'   the caller's unknowns; `dphi`, a list of two, the second derivative
#'   contracted against each direction, one row per observation;
#'   `dpsi`, the third derivative contracted against both, one row per
#'   observation; and `curvature`, the fourth derivative contracted
#'   against `g` and both directions.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:20, y = rnorm(20))
#' term <- term_build(linpar(~t), dd)
#' # an additive term bends no predictor, so every order above the first is
#' # zero, as the base method returns
#' out <- term_fourth(term, rep(0, 20), dd$y,
#'                    score = function(e, i) dd$y[i] - e,
#'                    curvature = function(e, i) -1, psi = list(),
#'                    g = rep(1, 20), seed = matrix(0, 20, 2),
#'                    blocks = function(e, i, D, act) NULL,
#'                    directions = list(c(1, 0), c(0, 1)))
#' all(out$curvature == 0)
#'
#' @seealso [term_third()] for the third order this differentiates,
#'   [term_curvature()] for the second, [statmodels7::reml()] for the
#'   criterion whose Hessian requires it.
#' @export
#' @aliases term_fourth.model_term term_fourth.structural_term
term_fourth <- S7::new_generic("term_fourth", "term",
  function(term, eta, y, score, curvature, psi, g, seed, blocks, directions,
           ...) S7::S7_dispatch())

S7::method(term_fourth, model_term) <- function(term, eta, y, score,
                                                curvature, psi, g, seed,
                                                blocks, directions, ...) {
  seed <- as.matrix(seed)
  m <- ncol(seed)
  z <- matrix(0, nrow(seed), m)
  list(jacobian = seed, dphi = list(z, z), dpsi = z,
       curvature = matrix(0, m, m))
}

S7::method(term_fourth, structural_term) <- function(term, eta, y, score,
                                                     curvature, psi, g, seed,
                                                     blocks, directions, ...) {
  stop(sprintf("the term class '%s' does not implement term_fourth().",
               attr(S7::S7_class(term), "name")), call. = FALSE)
}

S7::method(term_params, structural_term) <- function(term, ...) {
  stop(sprintf("the term class '%s' does not implement term_params().",
               attr(S7::S7_class(term), "name")), call. = FALSE)
}

S7::method(term_links, structural_term) <- function(term, ...) {
  stop(sprintf("the term class '%s' does not implement term_links().",
               attr(S7::S7_class(term), "name")), call. = FALSE)
}

S7::method(term_filter, structural_term) <- function(term, eta, y, score,
                                                     curvature, psi, ...) {
  stop(sprintf("the term class '%s' does not implement term_filter().",
               attr(S7::S7_class(term), "name")), call. = FALSE)
}

# the observation weights of a likelihood-shape term: one number for every
# row, or one per row
.structural_weights <- function(weights, n) {
  if (is.null(weights)) return(rep(1, n))
  w <- as.numeric(weights)
  if (!(length(w) %in% c(1L, n)) || anyNA(w)) {
    stop(sprintf("'weights' must be one number or %d non-missing numbers.", n),
         call. = FALSE)
  }
  rep_len(w, n)
}
