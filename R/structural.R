#' @include term_classes.R generics.R
NULL

#' @title Parameters of a Structural Term
#'
#' @description
#' The names of a structural term's own parameters, in the order its
#' filter expects them. A structural term contributes no design block, so
#' its parameters are not coefficients: they are estimated alongside the
#' distribution's, on the unconstrained scale its links define.
#'
#' @param term An object inheriting from \code{\link{structural_term}}.
#' @param ... Passed to methods.
#'
#' @return A character vector.
#'
#' @examples
#' term_params(gas(p = 1, q = 1))
#'
#' @seealso \code{\link{term_links}}, \code{\link{term_filter}}
#' @export
term_params <- S7::new_generic("term_params", "term",
  function(term, ...) S7::S7_dispatch())

#' @title Links of a Structural Term's Parameters
#'
#' @description
#' One \pkg{linkfunctions7} link per parameter of
#' \code{\link{term_params}}, carrying it to the unconstrained scale the
#' model layer optimizes on.
#'
#' @param term An object inheriting from \code{\link{structural_term}}.
#' @param ... Passed to methods.
#'
#' @return A named list of link objects.
#'
#' @examples
#' vapply(term_links(gas(p = 1, q = 1)), function(l) l@link_name, character(1))
#'
#' @seealso \code{\link{term_params}}
#' @export
term_links <- S7::new_generic("term_links", "term",
  function(term, ...) S7::S7_dispatch())

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
#' Writing \eqn{\eta_t^{0}} for the static predictor supplied in \code{eta}
#' and \eqn{\psi} for the term's parameters, the filter returns the pair
#'
#' \deqn{\eta_t = \eta_t^{0} + f_t(\psi),
#'   \qquad J_{tj} = \frac{\partial \eta_t}{\partial \psi_j},}
#'
#' where \eqn{f_t} is the term's own recursion, driven by
#' \code{score} and \code{curvature} evaluated at the predictor already
#' produced. Both are read at \eqn{\eta_t}, so \code{curvature} is the
#' second derivative \eqn{\partial^{2} \ell_t / \partial \eta^{2}} and is
#' negative at an ordinary observation; passing the information, its
#' negative, gives a predictor and a Jacobian that are internally
#' consistent and wrong.
#'
#' The derivative is returned because the recursion is the only place it
#' can be computed. A model layer differencing the filter would pay one
#' pass per parameter and inherit the error of the difference; propagating
#' the derivative alongside the state costs one extra vector per parameter
#' and is exact.
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
#'   \code{\link{term_params}}.
#' @param ... Passed to methods.
#'
#' @return A list with \code{eta}, the predictor the term produces, and
#'   \code{jacobian}, an \code{n} by \code{length(psi)} matrix of its
#'   derivatives with respect to \code{psi}.
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
#' @seealso \code{\link{gas}}
#' @export
term_filter <- S7::new_generic("term_filter", "term",
  function(term, eta, y, score, curvature, psi, ...) S7::S7_dispatch())

#' @title Which of a Term's Parameters Acts as an Intercept
#'
#' @description
#' The name of the parameter, if any, that shifts the term's contribution by
#' a constant, and is therefore the same direction as an intercept in the
#' equation the term sits in.
#'
#' @details
#' It exists so that a fitting layer can resolve the confounding rather than
#' refuse the model. A score-driven level \eqn{\omega} and a regime's first
#' level both add a constant to their equation's predictor: with an
#' intercept there too, adding \eqn{c} to one and subtracting the matching
#' amount from the other leaves every predictor unchanged, and the
#' likelihood is flat along that direction. Which of the two is dropped is
#' the layer's decision and not the term's, since only the layer knows what
#' else the equation carries.
#'
#' The base method returns \code{character(0)}, meaning the term shifts
#' nothing by a constant and no such question arises.
#'
#' @param term A term.
#' @param ... Passed to methods.
#'
#' @return A single string, or \code{character(0)}.
#'
#' @examples
#' term_level_param(gas(p = 1, q = 1))
#' term_level_param(linpar(~x))
#'
#' @seealso \code{\link{term_params}}
#' @export
term_level_param <- S7::new_generic("term_level_param", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_level_param, model_term) <- function(term, ...) character(0)


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
#' \code{param_readable()} does for a matrix parameter.
#'
#' The base method reports the parameters themselves on the PARAMETER scale,
#' with the diagonal Jacobian of their links, which is what every term whose
#' coordinates are already its quantities wants.
#'
#' @param term A built term.
#' @param zeta The term's parameters on the unconstrained scale, named as
#'   \code{\link{term_params}}.
#' @param ... Passed to methods.
#'
#' @return A list with \code{name}, \code{value}, \code{jacobian} (one row
#'   per quantity and one column per parameter) and \code{scale}, the link
#'   an interval for each quantity is built on.
#'
#' @examples
#' term_readable(gas(p = 1, q = 1), c(omega = 0.3, alpha1 = 0.4, pacf1 = 0.8))
#'
#' @seealso \code{\link{term_params}}, \code{\link{term_links}}
#' @export
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
#' \code{\link{term_filter}} returns the derivative of the predictor in the
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
#' \code{dscore}, the derivative of the objective in anything the score
#' depends on is
#'
#' \deqn{\frac{\partial L}{\partial \theta}
#'   = \left.\frac{\partial L}{\partial \theta}\right|_{\mathrm{direct}}
#'     + \sum_t \bar{s}_t \frac{\partial s_t}{\partial \theta},}
#'
#' so a model layer whose score is the derivative of its log-likelihood in
#' one distribution parameter obtains the derivative in the predictor of
#' ANOTHER by multiplying \code{dscore} by the mixed second derivative of
#' that log-likelihood. \code{deta} is that formula applied to the term's own
#' equation, where the second factor is the curvature.
#'
#' @param term A built structural term.
#' @param eta The static part of the predictor, as
#'   \code{\link{term_filter}} takes it.
#' @param y The response.
#' @param score,curvature The callbacks of \code{\link{term_filter}}.
#' @param psi The term's parameters, named as \code{\link{term_params}}.
#' @param g The direct derivative of the objective in the predictor the term
#'   produced, one value per observation.
#' @param ... Passed to methods.
#'
#' @return A list with \code{deta} and \code{dscore}, each one value per
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
#' @seealso \code{\link{term_filter}}
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
#' \code{\link{term_adjoint}} does not give.
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
#' \code{seed} is \eqn{\partial \eta^{0}/\partial u}, one row per
#' observation: the caller says how its unknowns reach the static predictor
#' and the term knows nothing else about them. \code{blocks} is where the
#' model's own derivatives enter, since the score the recursion is driven by
#' depends on every equation and not only on the predictor it is read at. It
#' is called once per observation with the predictor there and the Jacobian
#' the recursion has reached, and returns the two quantities that seed the
#' first and second derivatives of that score,
#'
#' \deqn{\dot S_t = \ell_{pp,t}D_t + \texttt{cross}, \qquad
#'   \ddot S_t = \ell_{pp,t}E_t + \texttt{M},}
#'
#' with \code{cross} \eqn{= \sum_{q\ne p}\ell_{pq,t}C_{q,t}} and \code{M}
#' \eqn{= \sum_{r,r'}\ell_{prr',t}V_{r,t}^\top V_{r',t}}, the third
#' derivatives of the log-density in the predictors. A model of one equation
#' has \code{cross} zero and \code{M} equal to \eqn{\ell_{ppp}D^\top D}.
#'
#' @param term A built structural term.
#' @param eta The static part of the predictor.
#' @param y The response.
#' @param score,curvature The callbacks of \code{\link{term_filter}}.
#' @param psi The term's parameters, named as \code{\link{term_params}}.
#' @param g The weights the second derivative is contracted against, one per
#'   observation.
#' @param seed The derivative of the static predictor in the caller's
#'   unknowns, one row per observation.
#' @param blocks A function of the predictor, the index and the current
#'   Jacobian row, returning \code{cross} and \code{M}.
#' @param ... Passed to methods.
#'
#' @return A list with \code{jacobian}, the derivative of the predictor in
#'   the caller's unknowns, and \code{curvature}, the contracted second
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
#' @seealso \code{\link{term_adjoint}}, \code{\link{term_filter}}
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
