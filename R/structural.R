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
#'                    psi = list(omega = 0.1, a1 = 0.3, pacf1 = 0.5))
#' head(out$eta, 3)
#' dim(out$jacobian)
#'
#' @seealso \code{\link{gas}}
#' @export
term_filter <- S7::new_generic("term_filter", "term",
  function(term, eta, y, score, curvature, psi, ...) S7::S7_dispatch())

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
