#' @include structural.R regime.R nonlinear.R
NULL

#' Charted Coordinates of a Term
#'
#' @description
#' Names the free values of a term that are coordinates of a chart onto a
#' bounded set, so that a coordinate running to infinity reaches the boundary of
#' the quantity that it describes. A statistic read at such a point requires a
#' check that the free scale cannot provide: the derivative on the free scale
#' equals the derivative on the bounded scale multiplied by the derivative of
#' the chart, and the derivative of the chart tends to zero at the edge, so the
#' derivative on the free scale vanishes there whatever the derivative on the
#' bounded scale is.
#'
#' @details
#' For a structural term the result is a character vector of names from
#' [term_params()]: every parameter whose link in [term_links()] is not the
#' identity, and for [regime()] also the additive log-ratios of the
#' transition matrix, whose chart is the one that
#' [parameters7::transition_matrix()] provides and which [term_links()]
#' therefore reports as the identity. For [nl()] it is an integer vector of
#' positions in the term's block, named by parameter: the scalar parameters
#' (those with no subformula) whose link is not the identity. For every other
#' term the result is an empty vector.
#'
#' @param term A model term. For [nl()] it must have been built.
#' @param ... Unused.
#'
#' @return A character vector for a structural term, an integer vector for
#'   [nl()] named by parameter (`integer(0)` when the term has no charted
#'   parameter), and `character(0)` otherwise.
#'
#' @seealso [term_links()] for the charts, [statmodels7::statmod_certificate()]
#'   for the check that reads this.
#'
#' @examples
#' term_charted(regime(2))
#' term_charted(gas(p = 1, q = 1))
#'
#' # A scalar parameter of nl() on a log link is charted; the position is
#' # its column in the term's block.
#' nb <- term_build(nl(~ a * exp(-r * x), start = list(a = 1, r = 1),
#'                     links = list(r = linkfunctions7::log_link())),
#'                  data.frame(x = seq(0.1, 2, length.out = 10)))
#' term_charted(nb)
#' @export
#' @aliases term_charted.model_term term_charted.structural_term
#'   term_charted.RegimeTerm term_charted.NlTerm
term_charted <- S7::new_generic("term_charted", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_charted, model_term) <- function(term, ...) character(0)

S7::method(term_charted, structural_term) <- function(term, ...) {
  lk <- term_links(term)
  names(lk)[!vapply(lk, function(l) identical(l@link_name, "identity"),
                    logical(1))]
}

S7::method(term_charted, RegimeTerm) <- function(term, ...) {
  lk <- term_links(term)
  nm <- names(lk)[!vapply(lk, function(l) identical(l@link_name, "identity"),
                          logical(1))]
  c(nm, term@chain@free_names)
}

S7::method(term_charted, NlTerm) <- function(term, ...) {
  bp <- term@build_info
  if (!length(bp)) {
    stop("the term has not been built; call term_build(term, data) first.",
         call. = FALSE)
  }
  out <- integer(0)
  for (p in bp$params) {
    if (!is.null(bp$subs[[p]])) next
    if (identical(.nl_link(bp$links, p)@link_name, "identity")) next
    ix <- bp$index[[p]]
    if (length(ix) == 1L) out[[p]] <- as.integer(ix)
  }
  out
}
