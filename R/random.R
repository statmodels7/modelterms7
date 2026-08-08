#' @include term_classes.R generics.R
NULL

#' @title S7 Class for Grouped Random-Effect Terms
#' @name RandomTerm
#'
#' @description
#' A subclass of \code{\link{additive_term}} for grouped coefficients with
#' a distribution on the effects: the indicator block of a grouping factor,
#' one coefficient per level, with the penalty carrying the effects'
#' distribution. Constructed by \code{\link{random}}.
#'
#' @inheritParams additive_term
#' @param formula The bar formula, e.g. \code{~ 1 | g}.
#' @param precision A \pkg{parameters7} matrix parameter for the effects'
#'   precision, or \code{NULL}.
#' @param distrib A univariate \pkg{distributions7} object for the
#'   effects, or \code{NULL}.
#' @param kinks The declared kink set of \code{distrib}'s log-density.
#'
#' @return An object of class \code{RandomTerm}.
#'
#' @seealso \code{\link{random}}
#' @examples
#' S7::S7_inherits(random(~ 1 | g), RandomTerm)
#' @export
RandomTerm <- S7::new_class(
  name = "RandomTerm",
  parent = additive_term,
  properties = list(
    formula = S7::class_any,
    precision = S7::class_any,
    distrib = S7::class_any,
    kinks = S7::class_numeric
  )
)

#' Grouped Random-Effect Term
#'
#' @description
#' Random intercepts for a grouping factor: \code{random(~ 1 | g)} builds
#' the indicator block of \code{g}, one coefficient per level, and
#' attaches the distribution of the effects as the penalty on those
#' coefficients -- which is what a random effect is under penalized
#' likelihood.
#'
#' @details
#' Three choices of effect distribution are available, mutually exclusive.
#' With neither argument the effects are independent Gaussians with one
#' free scale (the effects' standard deviation), through
#' \code{\link[penalties7]{ridge_penalty}}. With \code{precision}, a
#' \pkg{parameters7} matrix parameter of dimension equal to the number of
#' levels enters as the precision of a joint Gaussian on the effects,
#' through \code{\link[penalties7]{structured_penalty}}: the
#' hyperparameters are the structure's free values, so correlated effects
#' (an AR(1) over ordered groups, a compound symmetry) come from the
#' structure and not from the term. With \code{distrib}, a univariate
#' \pkg{distributions7} object -- holding its own location, typically
#' through \code{\link[distributions7]{fixed}} at zero -- is applied
#' coordinatewise through \code{\link[penalties7]{distrib_penalty}}. A
#' distribution used as a penalty gives joint-mode (penalized likelihood)
#' estimation of the effects; the marginal likelihood, which integrates
#' them out, is not provided here.
#'
#' Only random intercepts are implemented. A bar formula with covariates
#' on its left side (\code{~ x | g}) needs the block-diagonal composition
#' of per-group covariance structures, which \pkg{parameters7} does not
#' provide yet, and is rejected with that reason.
#'
#' Prediction maps new data onto the levels seen at build time; a level
#' the term has not seen is rejected.
#'
#' @param formula A bar formula, \code{~ 1 | g}, with \code{g} evaluating
#'   to the grouping variable in the data.
#' @param precision A \pkg{parameters7} matrix parameter of dimension
#'   equal to the number of levels, or \code{NULL}.
#' @param distrib A univariate \pkg{distributions7} object for the
#'   effects, or \code{NULL}.
#' @param kinks The kink set of \code{distrib}'s log-density in its
#'   argument, passed to \code{\link[penalties7]{distrib_penalty}}.
#' @param label A single non-empty string prefixed to the coefficient
#'   names.
#'
#' @return An object of class \code{\link{RandomTerm}} (a specification;
#'   see \code{\link{term_build}}).
#'
#' @examples
#' dd <- data.frame(y = rnorm(9), g = factor(rep(c("a", "b", "c"), 3)))
#' built <- term_build(random(~ 1 | g), dd)
#' term_coef_names(built)
#' term_penalty(built)@params
#'
#' @export
random <- function(formula, precision = NULL, distrib = NULL,
                   kinks = numeric(0), label = "random") {
  if (!inherits(formula, "formula") || length(formula) != 2L) {
    stop("'formula' must be a one-sided bar formula, e.g. ~ 1 | g.",
         call. = FALSE)
  }
  e <- formula[[2L]]
  if (!is.call(e) || !identical(e[[1L]], as.name("|"))) {
    stop("'formula' must contain a grouping bar, e.g. ~ 1 | g.",
         call. = FALSE)
  }
  if (!identical(e[[2L]], 1) && !identical(e[[2L]], 1L)) {
    stop(paste("only random intercepts (~ 1 | g) are implemented: a",
               "random slope needs the block-diagonal composition of",
               "per-group structures, which parameters7 does not provide",
               "yet."), call. = FALSE)
  }
  if (!is.null(precision) && !is.null(distrib)) {
    stop("'precision' and 'distrib' are mutually exclusive.", call. = FALSE)
  }
  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !nzchar(label)) {
    stop("'label' must be a single non-empty character string.",
         call. = FALSE)
  }
  RandomTerm(label = label, formula = formula,
             precision = precision, distrib = distrib, kinks = kinks,
             X = NULL, coef_names = character(0),
             blueprint = list(), penalty = NULL)
}

.random_group <- function(expr, data, levels = NULL) {
  v <- eval(expr, data, baseenv())
  if (is.null(levels)) return(factor(v))
  f <- factor(v, levels = levels)
  if (any(is.na(f) & !is.na(v))) {
    bad <- unique(as.character(v)[is.na(f) & !is.na(v)])
    stop(sprintf("grouping level '%s' was not present at build time.",
                 bad[1L]), call. = FALSE)
  }
  f
}

S7::method(term_build, RandomTerm) <- function(term, data, ...) {
  gexpr <- term@formula[[2L]][[3L]]
  g <- .random_group(gexpr, data)
  m <- nlevels(g)
  if (m < 2L) {
    stop("the grouping variable must have at least two levels.",
         call. = FALSE)
  }
  mf <- stats::model.frame(~ 0 + g, data.frame(g = g),
                           na.action = stats::na.pass)
  Z <- stats::model.matrix(~ 0 + g, mf)
  attr(Z, "assign") <- NULL
  attr(Z, "contrasts") <- NULL
  cn <- paste(term@label, levels(g), sep = ".")
  colnames(Z) <- cn
  rownames(Z) <- NULL

  pen <- if (!is.null(term@precision)) {
    p <- penalties7::structured_penalty(term@precision)
    if (p@n_coef != m) {
      stop(sprintf("'precision' has dimension %d and the grouping has %d levels.",
                   p@n_coef, m), call. = FALSE)
    }
    p
  } else if (!is.null(term@distrib)) {
    penalties7::distrib_penalty(term@distrib, n_coef = m,
                                kinks = term@kinks)
  } else {
    penalties7::ridge_penalty(n_coef = m)
  }

  term@X <- Z
  term@coef_names <- cn
  term@blueprint <- list(gexpr = gexpr, glevels = levels(g))
  term@penalty <- pen
  term
}

S7::method(term_predict, RandomTerm) <- function(term, newdata, ...) {
  .assert_built(term)
  bp <- term@blueprint
  g <- .random_group(bp$gexpr, newdata, levels = bp$glevels)
  mf <- stats::model.frame(~ 0 + g, data.frame(g = g),
                           na.action = stats::na.pass)
  Z <- stats::model.matrix(~ 0 + g, mf)
  attr(Z, "assign") <- NULL
  attr(Z, "contrasts") <- NULL
  colnames(Z) <- term@coef_names
  rownames(Z) <- NULL
  Z
}
