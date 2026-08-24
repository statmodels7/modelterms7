#' @include term_classes.R generics.R
NULL

#' @title Parameters of a Structural Term
#'
#' @description
#' The names of a structural term's own parameters, in the order its filter
#' and its derivative recursions expect them. A structural term contributes no
#' design block, so these are not coefficients: they are estimated beside the
#' distribution's, on the unconstrained scale [term_links()] defines, and
#' [term_npar()] counts them.
#'
#' @details
#' The names are the term's own vocabulary. [gas()] answers `omega` for the
#' level, `alpha1` ... `alphap` for the score loadings and `pacf1` ... `pacfq`
#' for the persistence; [regime()] answers its levels and the free entries of
#' its transition matrix. They are what indexes everything else about the
#' term: [term_start()] returns one value per name, [term_readable()] carries
#' them onto the quantities a reader reads, and a [term_penalties()] entry's
#' `index` gives positions in this vector.
#'
#' A **subformula** expands the parameter it develops in place, so
#' `gas(p = 1, q = 1)` has three parameters and
#' `gas(p = 1, q = 1, omega ~ z)` has four: `omega.(Intercept)` and `omega.z`
#' where the level was.
#'
#' The method on [structural_term()] throws, naming the class: a structural
#' class supplies this itself.
#'
#' @param term An object inheriting from [structural_term()]. A class that
#'   does not implement the generic throws
#'   `"the term class 'X' does not implement term_params()."`.
#' @param ... Passed to methods. No shipped method reads anything here.
#'
#' @return A character vector, one name per parameter, of length
#'   [term_npar()].
#'
#' @seealso [term_links()] for the chart each rides, [term_start()] for where
#'   they begin, [term_readable()] for the quantities they map to,
#'   [term_coef_names()] for the additive branch's equivalent.
#'
#' @examples
#' # The score-driven vocabulary: a level, a loading, a persistence.
#' term_params(gas(p = 1, q = 1))
#' term_params(gas(p = 2, q = 2))
#'
#' # A subformula expands the parameter it develops, in place.
#' set.seed(1)
#' d <- data.frame(y = rnorm(30), z = rnorm(30), t = 1:30)
#' term_params(term_build(gas(p = 1, q = 1, omega ~ z, time = t), d))
#'
#' # It is what term_npar() counts on this branch.
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
#' One \pkg{linkfunctions7} link per parameter of [term_params()], carrying
#' that parameter from its own admissible set onto the whole real line. A
#' fitting layer optimizes on that unconstrained scale, so an optimizer never
#' has to be told about the constraint.
#'
#' @details
#' The charts are chosen so that a proposal from anywhere lands somewhere
#' admissible. [gas()] puts its level on the identity, its score loadings on
#' the **log** so a loading stays positive, and its persistence coordinates on
#' the **rhobit** so each partial autocorrelation stays inside \eqn{(-1, 1)}
#' and the filter stays stationary. The `links` argument of a constructor
#' overrides them per parameter.
#'
#' The persistence is the case worth understanding. The stationary region in
#' the autoregressive coefficients is not a box, so no collection of scalar
#' links covers it; the partial autocorrelations are each in \eqn{(-1, 1)}
#' independently, and Levinson-Durbin carries them onto the coefficients.
#' That is why [term_readable()] exists: the coordinate and the quantity a
#' reader reads are different things.
#'
#' Where a parameter carries a subformula the link is applied **inside** the
#' development, \eqn{\psi_{j,t} = g_j^{-1}(z_t'\gamma_j)}, so the parameter
#' stays admissible at every observation and the coefficients \eqn{\gamma_j}
#' are unconstrained.
#'
#' The method on [structural_term()] throws, naming the class.
#'
#' @param term An object inheriting from [structural_term()]. A class that
#'   does not implement the generic throws
#'   `"the term class 'X' does not implement term_links()."`.
#' @param ... Passed to methods. No shipped method reads anything here.
#'
#' @return A named list of \pkg{linkfunctions7} link objects, one per entry of
#'   [term_params()] and named by it.
#'
#' @seealso [term_params()] for the names, [term_start()] for the point on
#'   this scale a fit begins at, [term_readable()] for the quantities the
#'   coordinates map to.
#'
#' @examples
#' # Identity for the level, log for the loading, rhobit for the persistence.
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
#' linkfunctions7::linkinv(lk$alpha1, c(-40, 0, 40))
#'
#' @export
#' @aliases term_links.structural_term
term_links <- S7::new_generic("term_links", "term",
  function(term, ...) S7::S7_dispatch())

#' @title Where a Term's Own Parameters Start
#'
#' @description
#' The starting values of a structural term's parameters, on the unconstrained
#' scale [term_links()] defines, one per name [term_params()] gives. A fitting
#' layer reads it to begin a search.
#'
#' @details
#' The start belongs to the term because only the term knows what a coordinate
#' of zero means on each of its charts. The base method returns zero
#' everywhere, which is each link's own natural point, and that is right
#' wherever zero means "the model without this term".
#'
#' It does not always. [gas()] overrides it: its score loadings ride a log
#' chart, so a coordinate of zero is a loading of **one**, which is a strongly
#' driven filter and a poor place to begin. It starts them at 0.1 through the
#' chart, \eqn{\log 0.1 = -2.303}, and leaves every other coordinate at zero.
#'
#' A term that needs the data to place its start, as a marginal break-point
#' term does, computes it at [term_build()] and returns it from here.
#'
#' @param term A built structural term.
#' @param ... Passed to methods. No shipped method reads anything here.
#'
#' @return A named numeric vector on the unconstrained scale, of length
#'   [term_npar()] and named as [term_params()].
#'
#' @seealso [term_links()] for the scale it is on, [term_params()] for the
#'   names, [term_coef_start()] for the additive branch's equivalent.
#'
#' @examples
#' # The level and the persistence start at zero; the loading does not.
#' term_start(gas(p = 1, q = 1))
#'
#' # Because zero on a log chart is a loading of one.
#' lk <- term_links(gas(p = 1, q = 1))
#' linkfunctions7::linkinv(lk$alpha1, term_start(gas(p = 1, q = 1))[["alpha1"]])
#'
#' # One value per parameter, whatever the order.
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

#' @title Apply a Structural Term to a Linear Predictor
#'
#' @description
#' Runs the term's recursion over the data and returns the predictor it
#' produces, together with the derivative of that predictor with respect
#' to the term's own parameters. This is the operation that makes a
#' structural term structural: the predictor at one observation depends on
#' the others, so it cannot be written as a block of columns.
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
#' an ordinary observation. **Its sign is load-bearing**: passing the
#' information, which is its negative, returns a predictor and a Jacobian
#' that are internally consistent and wrong, with no error to say so.
#'
#' The derivative is returned because the recursion is the only place it can
#' be computed. Propagating it beside the state costs one extra vector per
#' parameter and is exact; a model layer differencing the filter instead
#' would pay one pass per parameter and inherit the error of the difference.
#'
#' @param term A built structural term.
#' @param eta The static part of the linear predictor, one value per
#'   observation.
#' @param y The response.
#' @param score A function of the predictor returning the derivative of
#'   the log-likelihood with respect to it, one value per observation.
#' @param curvature A function of the predictor returning the second
#'   derivative of the log-likelihood with respect to it.
#' @param psi The term's parameters, on the parameter scale, named as
#'   [term_params()].
#' @param ... Passed to methods.
#'
#' @return A list with `eta`, the predictor the term produces,
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
#' # the score and curvature a Gaussian mean would supply
#' out <- term_filter(term, eta = rep(0, 20), y = dd$y,
#'                    score = function(e, i) dd$y[i] - e,
#'                    curvature = function(e, i) -1,
#'                    psi = list(omega = 0.1, alpha1 = 0.3, pacf1 = 0.5))
#' head(out$eta, 3)
#' dim(out$jacobian)
#'
#' @seealso [gas()] and [regime()] for the two structural shapes,
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
#' How the predictor a structural term produces moves when the static part
#' of the predictor moves, one row per observation and one column per
#' direction the caller supplies.
#'
#' @details
#' A score-driven term's level is driven by scores read AT the predictor the
#' recursion is producing, so a coefficient in the same equation reaches the
#' level as well as the static part: writing \eqn{f_t} for the level and
#' \eqn{x_t} for a row of the design,
#' \deqn{\frac{\partial \eta_t}{\partial \beta} = x_t +
#'   \frac{\partial f_t}{\partial \beta},}
#' and the second piece obeys the recursion the filter already runs,
#' \deqn{\frac{\partial f_t}{\partial \beta} =
#'   \sum_i \alpha_i \, \ell''_{t-i}
#'     \frac{\partial \eta_{t-i}}{\partial \beta} +
#'   \sum_j \beta_j \frac{\partial f_{t-j}}{\partial \beta}.}
#' The curvature it needs is the one [term_filter()] returns, so
#' no callback is evaluated here and the pass is arithmetic alone.
#'
#' Without it a standard error of the predictor counts the static part only.
#' Measured on a score-driven mean with one covariate beside it, that
#' understates the standard error by about a quarter.
#'
#' The base method returns `NULL`. A term that is not a filter carries no
#' state, so the derivative is the design row itself and the caller needs
#' nothing from the term.
#'
#' @param term A built structural term.
#' @param curv The curvature at each predictor, as `term_filter`
#'   returns it.
#' @param X The directions to propagate, one column each. Ordinarily the
#'   equation's design, so the result is one row per observation and one
#'   column per coefficient.
#' @param psi The term's parameters, on the parameter scale, named as
#'   [term_params()].
#' @param ... Passed to methods.
#'
#' @return A matrix of `X`'s dimensions, holding \eqn{\partial f_t/\partial
#'   \beta} for each column of `X`, so the full derivative of the predictor
#'   is `X + result`. `NULL` for a term that carries no state.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:20, y = rnorm(20), x = rnorm(20))
#' term <- term_build(gas(p = 1, q = 1, time = t), dd)
#' psi <- list(omega = 0.1, alpha1 = 0.3, pacf1 = 0.5)
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
#' The contribution a structural term makes at rows that come after the ones
#' it was built on, continuing its recursion rather than restarting it.
#'
#' @details
#' A structural term's contribution at one observation is not a function of
#' that observation: it is the state a recursion has reached, so predicting
#' past the series means carrying the state forward. What makes it possible
#' without simulation is that the quantity driving the recursion has zero
#' conditional mean, the score itself for a score-driven term, so beyond the
#' data the recursion is deterministic.
#'
#' The base method signals an error. A term with state that cannot say what
#' its state does next has nothing to offer a prediction, and returning zero
#' there would read as a term with no effect at all.
#'
#' @param term A built structural term.
#' @param psi The term's parameters, named as [term_params()].
#' @param f_past The term's contribution at each observed row.
#' @param s_past The driving quantity at each observed row.
#' @param newdata The rows to continue onto.
#' @param ... Passed to methods.
#'
#' @return A numeric vector of `nrow(newdata)` contributions.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:20, y = rnorm(20))
#' term <- term_build(gas(p = 1, q = 1, time = t), dd)
#' psi <- list(omega = 0.1, alpha1 = 0.3, pacf1 = 0.5)
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
#' The predictor a structural term produces when the response is being
#' GENERATED rather than read, together with whatever latent quantity the
#' term drew on the way.
#'
#' @details
#' Simulating from a model that carries state is not the same operation as
#' fitting one, and the difference is which direction the response moves in.
#' A term whose contribution does not read the response can report that
#' contribution and leave the drawing to the caller: a latent chain's levels
#' and a group's break-point drawn from its prior are both like that. A
#' score-driven term cannot. Its level at one time is driven by the score of
#' the response at the time before, so the response has to be drawn AS the
#' recursion runs.
#'
#' One contract covers both. The caller supplies `draw`, a function of
#' a predictor and a row index returning one response value, and the method
#' returns the predictor it produced; a method that drew returns the
#' responses as well and one that did not returns `NULL` there, leaving
#' the caller to draw at the predictor.
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
#' out <- term_simulate(term, list(omega = 0.5, alpha1 = 0.3, pacf1 = 0.6),
#'                      rep(0, 30),
#'                      draw = function(e, i) stats::rnorm(1, e, 1))
#' head(out$y, 3)
#'
#' @seealso [term_filter()], [term_continue()]
#' @export
term_simulate <- S7::new_generic("term_simulate", "term",
  function(term, psi, eta, draw, ...) S7::S7_dispatch())

#' @name term_simulate.model_term
#' @title A Term Without State Draws Nothing
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
  stop(sprintf("'%s' does not say how a response is drawn from it.",
               class(term)[[1L]]), call. = FALSE)
}

#' @name term_continue.model_term
#' @title A Term Without State Is Not Continued
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
  stop(sprintf("'%s' does not say how its contribution continues past the series.",
               class(term)[[1L]]), call. = FALSE)
}

#' @name term_static_deriv.model_term
#' @title No State, No Propagation
#' @description
#' The base method returns `NULL`: an ordinary term's contribution at
#' one observation reads that observation alone.
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
#' The name of the parameter, if any, that shifts the term's contribution by
#' a constant, and is therefore the same direction as an intercept in the
#' equation the term sits in.
#'
#' @details
#' It exists so that a fitting layer can resolve the confounding instead of
#' refusing the model. A score-driven level \eqn{\omega} and a regime's first
#' level both add a constant to their equation's predictor: with an
#' intercept there too, adding \eqn{c} to one and subtracting the matching
#' amount from the other leaves every predictor unchanged, and the
#' likelihood is flat along that direction. Which of the two is dropped is
#' the layer's decision and not the term's, since only the layer knows what
#' else the equation carries.
#'
#' The base method returns `character(0)`, meaning the term shifts
#' nothing by a constant and no such question arises.
#'
#' @param term A term.
#' @param ... Passed to methods.
#'
#' @return A single string, or `character(0)`.
#'
#' @examples
#' term_level_param(gas(p = 1, q = 1))
#' term_level_param(linpar(~x))
#'
#' @seealso [term_level_design()] for the same question with the level
#'   developed, [term_params()] for the names it answers among.
#' @export
#' @aliases term_level_param.model_term
term_level_param <- S7::new_generic("term_level_param", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_level_param, model_term) <- function(term, ...) character(0)

#' @title The Design of a Term's Level Development
#'
#' @description
#' Where the parameter [term_level_param()] names is developed
#' with covariates, the design of that development, with one column per
#' coordinate named as [term_params()] names it. `NULL`
#' for a scalar level and for every other term.
#'
#' @details
#' It exists for the subspace form of the confounding question.
#' `term_level_param` answers for the constant: a coordinate whose
#' column is constant shifts the equation's predictor exactly as an
#' intercept does. With the level developed, a direction of the
#' development's span that also lies in the span of the equation's design
#' raises the same question for that direction, and only a fitting layer,
#' which holds both designs, can ask it. This generic hands it the one
#' half it cannot see.
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
#'   the names its columns carry.
#' @export
#' @aliases term_level_design.model_term
term_level_design <- S7::new_generic("term_level_design", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_level_design, model_term) <- function(term, ...) NULL


#' @title The Quantities a Fitted Term Reports
#'
#' @description
#' What a reader reads, with the Jacobian from the term's own parameters on
#' the unconstrained scale, so that a caller holding their variance matrix
#' can carry it across by the delta method.
#'
#' @details
#' A term's parameters are the coordinates it is ESTIMATED on, chosen so
#' that a search runs unconstrained, and they are not always the quantities
#' the model is about. The clearest case is a score-driven persistence,
#' which rides a partial autocorrelation because the stationary region is
#' not a box: what the literature writes as \eqn{\beta_1} is the
#' autoregressive coefficient, which is a function of the whole chart and
#' coincides with it only at \eqn{q = 1}. Reporting the coordinate under the
#' coefficient's name would promise one quantity and print another.
#'
#' Each row gives a value and the row of \eqn{\partial(\text{value}) /
#' \partial\zeta} at the current parameters, so a standard error is
#' \eqn{\sqrt{J V J^\top}} and an interval is built on whichever scale keeps
#' the quantity in its own set, exactly as \pkg{parameters7}'s
#' `param_readable()` does for a matrix parameter.
#'
#' The base method reports the parameters themselves on the PARAMETER scale,
#' with the diagonal Jacobian of their links, which is what every term whose
#' coordinates are already its quantities wants.
#'
#' @param term A built term.
#' @param zeta The term's parameters on the unconstrained scale, named as
#'   [term_params()].
#' @param ... Passed to methods.
#'
#' @return A list with `name`, `value`, `jacobian` (one row
#'   per quantity and one column per parameter) and `scale`, the link
#'   an interval for each quantity is built on.
#'
#' @examples
#' term_readable(gas(p = 1, q = 1), c(omega = 0.3, alpha1 = 0.4, pacf1 = 0.8))
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
#' The derivative of a caller's objective with respect to the static
#' predictor the term was handed, and with respect to the sequence of scores
#' it was given, both accounting for the recursion.
#'
#' @details
#' [term_filter()] returns the derivative of the predictor in the
#' term's OWN parameters, which is what estimating those needs. It is not
#' what estimating the coefficients of the same equation needs: the level at
#' one time is driven by the scores at earlier ones, which are read at
#' predictors those coefficients also enter, so the derivative of the
#' predictor in a coefficient carries a term the block does not. Propagating
#' that forward would cost one derivative array per coefficient; the reverse
#' recursion here costs one pass whatever their number, and is exact.
#'
#' Two derivatives are returned rather than one because the score the caller
#' supplies depends on more than the predictor it is read at. Writing
#' \eqn{s_t} for that score and \eqn{\bar{s}_t} for
#' `dscore`, the derivative of the objective in anything the score
#' depends on is
#'
#' \deqn{\frac{\partial L}{\partial \theta}
#'   = \left.\frac{\partial L}{\partial \theta}\right|_{\mathrm{direct}}
#'     + \sum_t \bar{s}_t \frac{\partial s_t}{\partial \theta},}
#'
#' so a model layer whose score is the derivative of its log-likelihood in
#' one distribution parameter obtains the derivative in the predictor of
#' ANOTHER by multiplying `dscore` by the mixed second derivative of
#' that log-likelihood. `deta` is that formula applied to the term's own
#' equation, where the second factor is the curvature.
#'
#' @param term A built structural term.
#' @param eta The static part of the predictor, as
#'   [term_filter()] takes it.
#' @param y The response.
#' @param score,curvature The callbacks of [term_filter()].
#' @param psi The term's parameters, named as [term_params()].
#' @param g The direct derivative of the objective in the predictor the term
#'   produced, one value per observation.
#' @param ... Passed to methods.
#'
#' @return A list with `deta` and `dscore`, each one value per
#'   observation.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:20, y = rnorm(20))
#' term <- term_build(gas(p = 1, q = 1, time = t), dd)
#' out <- term_adjoint(term, eta = rep(0, 20), y = dd$y,
#'                     score = function(e, i) dd$y[i] - e,
#'                     curvature = function(e, i) -1,
#'                     psi = list(omega = 0.1, alpha1 = 0.3, pacf1 = 0.5),
#'                     g = rep(1, 20))
#' head(out$deta, 3)
#'
#' @seealso [term_filter()]
#' @export
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
#' The Jacobian of the predictor the term produces with respect to a
#' caller's unknowns, and the second derivative of that predictor contracted
#' against a caller's weights. It is what an observed information needs and
#' [term_adjoint()] does not give.
#'
#' @details
#' The gradient of a model carrying a filter needs no second derivative and
#' no forward Jacobian: the reverse recursion answers it in one pass. The
#' curvature needs both. Writing \eqn{u} for the caller's unknowns,
#' \eqn{D_t = \partial e_t/\partial u} and \eqn{E_t = \partial^2
#' e_t/\partial u \partial u^\top}, the observed information of the model is
#'
#' \deqn{\sum_t w_t \sum_{q,r} \ell_{qr,t} V_{q,t}^\top V_{r,t}
#'   + \sum_t w_t \ell_{p,t} E_t,}
#'
#' whose first sum is a weighted crossproduct the caller assembles and whose
#' second is what this returns, the weights \eqn{g_t = w_t\ell_{p,t}} being
#' supplied.
#'
#' `seed` is \eqn{\partial \eta^{0}/\partial u}, one row per
#' observation: the caller says how its unknowns reach the static predictor
#' and the term knows nothing else about them. `blocks` is where the
#' model's own derivatives enter, since the score the recursion is driven by
#' depends on every equation and not only on the predictor it is read at. It
#' is called once per observation with the predictor there and the Jacobian
#' the recursion has reached, and returns the two quantities that seed the
#' first and second derivatives of that score,
#'
#' \deqn{\dot S_t = \ell_{pp,t}D_t + \texttt{cross}, \qquad
#'   \ddot S_t = \ell_{pp,t}E_t + \texttt{M},}
#'
#' with `cross` \eqn{= \sum_{q\ne p}\ell_{pq,t}C_{q,t}} and `M`
#' \eqn{= \sum_{r,r'}\ell_{prr',t}V_{r,t}^\top V_{r',t}}, the third
#' derivatives of the log-density in the predictors. A model of one equation
#' has `cross` zero and `M` equal to \eqn{\ell_{ppp}D^\top D}.
#'
#' @param term A built structural term.
#' @param eta The static part of the predictor.
#' @param y The response.
#' @param score,curvature The callbacks of [term_filter()].
#' @param psi The term's parameters, named as [term_params()].
#' @param g The weights the second derivative is contracted against, one per
#'   observation.
#' @param seed The derivative of the static predictor in the caller's
#'   unknowns, one row per observation.
#' @param blocks A function of the predictor, the index and the current
#'   Jacobian row, returning `cross` and `M`.
#' @param ... Passed to methods.
#'
#' @return A list with `jacobian`, the derivative of the predictor in
#'   the caller's unknowns, and `curvature`, the contracted second
#'   derivative.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:20, y = rnorm(20))
#' term <- term_build(gas(p = 1, q = 1, time = t), dd)
#'
#' # a gaussian mean of unit variance and one equation, so the model
#' # contributes no cross term and no third derivative
#' m <- 3L
#' out <- term_curvature(
#'   term, eta = rep(0, 20), y = dd$y,
#'   score = function(e, i) dd$y[i] - e,
#'   curvature = function(e, i) -1,
#'   psi = list(omega = 0.1, alpha1 = 0.3, pacf1 = 0.5),
#'   g = rep(1, 20), seed = matrix(0, 20, m),
#'   blocks = function(e, i, D) list(cross = numeric(m),
#'                                   M = matrix(0, m, m)))
#' dim(out$jacobian)
#' dim(out$curvature)
#'
#' @seealso [term_adjoint()], [term_filter()]
#' @export
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
#' The second derivative of [term_curvature()] differentiated once
#' more along a single direction, and the derivative of the term's Jacobian
#' along that same direction. It is what the exact gradient of a marginal
#' criterion needs when a penalty covers the term's own parameters.
#'
#' @details
#' A marginal criterion carries \eqn{-\frac{1}{2}\log|K|} at the penalized
#' mode, so its gradient in a hyperparameter needs \eqn{\mathrm{tr}(M\,\partial
#' K/\partial u\,[v])}, with \eqn{v} the direction the mode moves in. Where
#' the predictor is \eqn{X\beta} that contraction is the family's third
#' derivative against the design and nothing else. Where a filter produces the
#' predictor, \eqn{K} carries \eqn{\sum_t w_t \ell_{p,t}E_t} as well, and
#' differentiating it asks for \eqn{\partial^3 e_t/\partial u^3}.
#'
#' **Only contracted.** The full third derivative is an \eqn{m^3} array
#' per observation and is never formed: the criterion asks for it along the
#' one direction the mode moves in, so what is propagated is a matrix per
#' observation, the same size as the curvature and therefore the same
#' \eqn{O(nm^2)}. A caller wanting several hyperparameters calls this once
#' per direction, which is cheaper than one \eqn{O(nm^3)} assembly whenever
#' the hyperparameters are fewer than the unknowns.
#'
#' **What the model supplies.** Each order of differentiation of the
#' predictor pulls in one more order of the family, the recursion being driven
#' by a score read at the predictor it produces: the curvature's `M` is
#' built from third derivatives, and this needs a fourth. `blocks`
#' therefore returns two quantities beyond the curvature's, at an observation
#' with the current Jacobian row \eqn{D} and its directional derivative:
#'
#' \deqn{\texttt{dcurv} = \sum_r \ell_{ppr}V_r, \qquad
#'   \texttt{N} = \sum_{r,r'}\Big(\sum_{r''}\ell_{prr'r''}(V_{r''}\cdot v)\Big)
#'     V_r^\top V_{r'},}
#'
#' `dcurv` serving both the derivative of the curvature along \eqn{v},
#' which is \eqn{\texttt{dcurv}\cdot v}, and the two terms differentiating
#' \eqn{M}'s own \eqn{V_p}.
#'
#' The base method returns zeros, so an additive term -- whose second
#' derivative is already zero -- is covered without writing anything, and a
#' term written later that does not implement this reports no third
#' derivative rather than a wrong one.
#'
#' @param term A built term.
#' @param eta The static part of the predictor.
#' @param y The response.
#' @param score,curvature The callbacks of [term_filter()].
#' @param psi The term's parameters, named as [term_params()].
#' @param g The weights the third derivative is contracted against, one per
#'   observation.
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
#' # zero and the base method says so
#' out <- term_third(term, rep(0, 20), dd$y,
#'                   score = function(e, i) dd$y[i] - e,
#'                   curvature = function(e, i) -1, psi = list(),
#'                   g = rep(1, 20), seed = matrix(0, 20, 2),
#'                   blocks = function(e, i, D, act) NULL,
#'                   direction = c(1, 0))
#' all(out$curvature == 0)
#'
#' @seealso [term_curvature()], [term_adjoint()]
#' @export
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
