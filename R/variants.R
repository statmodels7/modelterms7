#' @include term_classes.R generics.R
NULL

#' @title The Term That Stands in for a Term During a Pilot Fit
#'
#' @description
#' Returns the unbuilt term with which a fitting layer replaces `term` in a
#' pilot fit, or `NULL` when the term needs no pilot. A term returns one when
#' it carries a discrete choice that the data must settle (the index of the
#' mode of [shape_unimodal()]), which the layer makes with
#' [term_variants()] from the contribution that the stand-in fits.
#'
#' @param term A term, built or not.
#' @param ... Passed to methods.
#'
#' @return `NULL`, or an unbuilt term carrying the label of `term`.
#'
#' @seealso [term_variants()], which reads the pilot.
#'
#' @examples
#' term_pilot(s(x))
#' class(term_pilot(scs(x, shape_unimodal())))[1]
#'
#' # A mode given in the shape needs no pilot.
#' term_pilot(scs(x, shape_unimodal(mode = 0.4)))
#'
#' @export
term_pilot <- S7::new_generic("term_pilot", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_pilot, model_term) <- function(term, ...) NULL


#' @title The Versions of a Term That a Fitting Layer Tries
#'
#' @description
#' Returns the unbuilt versions of `term` among which a fitting layer chooses
#' by its criterion, given the contribution of the term's pilot (see
#' [term_pilot()]) at the rows on which `term` was built. The layer fits the
#' model with each and keeps the best.
#'
#' @details
#' The base method returns `term` alone. For [scs()] with [shape_unimodal()]
#' and no mode, the pilot's maximum gives the Greville abscissa
#' \eqn{\xi_{M_0}} nearest to it, and the versions are the shape with the mode
#' at the indices \eqn{M_0 - 2, \dots, M_0 + 2} that exist, named `"M = j"`.
#'
#' @param term A built term.
#' @param pilot A numeric vector, the contribution of the pilot at the rows on
#'   which `term` was built, or `NULL`.
#' @param ... Passed to methods.
#'
#' @return A non-empty named list of unbuilt terms carrying the label of
#'   `term`.
#'
#' @seealso [term_pilot()].
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = runif(100))
#' tm <- term_build(scs(x, shape_unimodal(), basis7::bspline_smooth(k = 10)), d)
#' names(term_variants(tm, pilot = -(d$x - 0.3)^2))
#'
#' @export
term_variants <- S7::new_generic("term_variants", "term",
  function(term, pilot = NULL, ...) S7::S7_dispatch())

S7::method(term_variants, model_term) <- function(term, pilot = NULL, ...) {
  list(term)
}
