#' @title S7 Base Class for Model Terms
#' @name model_term
#'
#' @description
#' The abstract root of the term hierarchy. A model term records what a
#' formula names: the recipe that turns a data frame into a contribution to
#' the model, together with the metadata a fit reads. A term as written in a
#' formula is a specification; \code{\link{term_build}} turns it into a
#' built term carrying its design block.
#'
#' The formula interpreter (\code{\link{interpret_formula}}) recognizes a
#' term by this class: any call in a formula that evaluates to an object
#' inheriting from \code{model_term} is treated as a term, so a term class
#' defined outside the package works in a formula without registration.
#'
#' @param label A character string prefixed to the term's coefficient names
#'   when non-empty.
#' @param hyper The hyperparameters of the term's penalty that the
#'   caller HELD, as a named list. Empty, the default, means every one
#'   of them is estimated. See \code{\link{term_hyper}}.
#' @param grid How many values a path visits for each of the term's
#'   hyperparameters, as a named list. Empty, the default, leaves it to the
#'   criterion. See \code{\link{term_grid}}.
#' @param min_ratio How far down the path over the size of the kink
#'   reaches, as a fraction of the value that empties the block, or
#'   \code{numeric(0)} for the criterion's own. See
#'   \code{\link{term_path_min}}.
#'
#' @return An object inheriting from class \code{model_term}.
#'
#' @seealso \code{\link{additive_term}}, \code{\link{structural_term}},
#'   \code{\link{term_build}}
#' @examples
#' S7::S7_inherits(linpar(~1), model_term)
#' @export
model_term <- S7::new_class(
  name = "model_term",
  abstract = TRUE,
  properties = list(
    label = S7::class_character,
    # Which of its penalty's hyperparameters the caller has HELD, and at
    # what. Empty is the default and means every one of them is estimated:
    # which are and which are not is a property of the term, since the term
    # is where the penalty is named, and not of whatever criterion the fit
    # happens to run. See term_hyper().
    hyper = S7::new_property(S7::class_list, default = quote(list())),
    # How many values a PATH visits for each of them. A hyperparameter with
    # a kink is chosen by sweeping its own values, and how finely is a
    # property of the term for the same reason as which ones are estimated:
    # a block of four columns and one of four hundred want different grids,
    # and the criterion does not know which it is looking at.
    grid = S7::new_property(S7::class_list, default = quote(list())),
    # How far DOWN the path over the size of the kink reaches, as a
    # fraction of the value that empties the block. It belongs beside the
    # grid size and for the same reason, and it is one number rather than
    # one per hyperparameter because only that path uses it: a bounded
    # hyperparameter is swept over its own interval and a shape over a
    # geometric grid above its lower bound, where a ratio means nothing.
    min_ratio = S7::new_property(S7::class_numeric,
                                 default = quote(numeric(0)))
  )
)

#' @title S7 Class for Additive Terms
#' @name additive_term
#'
#' @description
#' The branch of \code{\link{model_term}} for terms that contribute a
#' linear block \eqn{X_j \beta_j} to the linear predictor of one
#' distribution parameter. A built additive term carries the design block,
#' the coefficient names, the blueprint that reproduces the mapping on new
#' data, and the penalty attached to its coefficients (\code{NULL} when the
#' term is unpenalized).
#'
#' @param label A character string prefixed to the coefficient names when
#'   non-empty.
#' @param hyper The hyperparameters of the term's penalty that the
#'   caller HELD, as a named list. Empty, the default, means every one
#'   of them is estimated. See \code{\link{term_hyper}}.
#' @param grid How many values a path visits for each of the term's
#'   hyperparameters, as a named list. Empty, the default, leaves it to the
#'   criterion. See \code{\link{term_grid}}.
#' @param min_ratio How far down the path over the size of the kink
#'   reaches, as a fraction of the value that empties the block, or
#'   \code{numeric(0)} for the criterion's own. See
#'   \code{\link{term_path_min}}.
#' @param X The design block, filled by \code{\link{term_build}}.
#' @param coef_names The coefficient names, filled by
#'   \code{\link{term_build}}.
#' @param blueprint The information needed to reproduce the mapping on new
#'   data, filled by \code{\link{term_build}}.
#' @param penalty The penalty on the term's coefficients, or \code{NULL}.
#'
#' @return An object inheriting from class \code{additive_term}.
#'
#' @seealso \code{\link{linpar}}, \code{\link{term_build}},
#'   \code{\link{term_predict}}
#' @examples
#' S7::S7_inherits(linpar(~1), additive_term)
#' @export
additive_term <- S7::new_class(
  name = "additive_term",
  parent = model_term,
  abstract = TRUE,
  properties = list(
    X = S7::class_any,
    coef_names = S7::class_character,
    blueprint = S7::class_list,
    penalty = S7::class_any
  )
)

#' @title S7 Class for Structural Terms
#' @name structural_term
#'
#' @description
#' The branch of \code{\link{model_term}} for terms that rewrite the
#' likelihood contributions rather than adding a design block:
#' \code{\link{gas}}, whose predictor is a recursion, and
#' \code{\link{regime}}, whose contribution is a likelihood mixed over
#' latent states. A term on this branch reports its own parameters through
#' \code{\link{term_params}} rather than coefficients, so
#' \code{\link{term_npar}} counts those and \code{\link{term_penalties}}
#' indexes into them.
#'
#' @param label A character string naming the term.
#' @param hyper The hyperparameters of the term's penalty that the
#'   caller HELD, as a named list. Empty, the default, means every one
#'   of them is estimated. See \code{\link{term_hyper}}.
#' @param grid How many values a path visits for each of the term's
#'   hyperparameters, as a named list. Empty, the default, leaves it to the
#'   criterion. See \code{\link{term_grid}}.
#' @param min_ratio How far down the path over the size of the kink
#'   reaches, as a fraction of the value that empties the block, or
#'   \code{numeric(0)} for the criterion's own. See
#'   \code{\link{term_path_min}}.
#'
#' @return An object inheriting from class \code{structural_term}.
#'
#' @seealso \code{\link{model_term}}
#' @examples
#' S7::S7_inherits(linpar(~1), structural_term)
#' @export
structural_term <- S7::new_class(
  name = "structural_term",
  parent = model_term,
  abstract = TRUE
)

#' @title S7 Class for the Unpenalized Parametric Term
#' @name LinparTerm
#'
#' @description
#' A subclass of \code{\link{additive_term}} for the unpenalized parametric
#' block built from a one-sided formula through \code{\link[stats]{model.matrix}}.
#' Constructed by \code{\link{linpar}}, and implicitly by
#' \code{\link{interpret_formula}}, which collects the bare covariates of a
#' model formula into one term of this class.
#'
#' @inheritParams additive_term
#' @param formula The one-sided formula defining the block.
#'
#' @return An object of class \code{LinparTerm}.
#'
#' @seealso \code{\link{linpar}}
#' @examples
#' S7::S7_inherits(linpar(~1), LinparTerm)
#' @export
LinparTerm <- S7::new_class(
  name = "LinparTerm",
  parent = additive_term,
  properties = list(
    formula = S7::class_any
  )
)
