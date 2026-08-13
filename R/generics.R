#' @include term_classes.R
NULL

#' @title Build a Term on Data
#'
#' @description
#' Turns a term specification into a built term: the design block is
#' computed from the data, the coefficient names are assigned, and the
#' blueprint that reproduces the mapping on new data is recorded. The
#' returned object is a copy of the specification with those properties
#' filled; the specification itself is unchanged.
#'
#' @details
#' An additive term contributes to the linear predictor through a design
#' block and, when it is penalized, a penalty on the coefficients of that
#' block:
#'
#' \deqn{\eta = \sum_{t} X_t \beta_t,
#'   \qquad \text{penalized objective} \quad
#'   -\ell(\beta) + \sum_{t} \rho_t(\beta_t; \theta_t),}
#'
#' and building the term is what produces \eqn{X_t} from the data and
#' attaches \eqn{\rho_t}. \code{\link{term_matrix}} reads the block,
#' \code{\link{term_penalty}} the penalty and \code{\link{term_predict}}
#' reproduces \eqn{X_t} on new rows through the blueprint. A structural
#' term is the exception: its contribution cannot be written as a block of
#' columns, and it reports itself through \code{\link{term_filter}} or
#' \code{\link{term_loglik}} instead.
#'
#' @param term An object inheriting from class \code{\link{model_term}}.
#' @param data A data frame.
#' @param ... Passed to methods.
#'
#' @return A built term of the same class as \code{term}.
#'
#' @examples
#' built <- term_build(linpar(~x), data.frame(x = 1:4))
#' term_matrix(built)
#'
#' @seealso \code{\link{term_predict}}, \code{\link{term_refresh}}, \code{\link{term_matrix}}, \code{\link{term_coef_names}}, \code{\link{term_npar}}, \code{\link{term_is_built}}
#' @export
term_build <- S7::new_generic("term_build", "term",
  function(term, data, ...) {
    if (!is.data.frame(data)) {
      stop("'data' must be a data frame.", call. = FALSE)
    }
    S7::S7_dispatch()
  })

S7::method(term_build, model_term) <- function(term, data, ...) {
  stop(sprintf("the term class '%s' does not implement term_build().",
               attr(S7::S7_class(term), "name")), call. = FALSE)
}

S7::method(term_build, structural_term) <- function(term, data, ...) {
  stop("structural terms are reserved for a later release; none is implemented yet.",
       call. = FALSE)
}

#' @title Whether a Term Has Been Built
#'
#' @description
#' \code{TRUE} for a term returned by \code{\link{term_build}} and
#' \code{FALSE} for a bare specification. The accessors
#' \code{\link{term_matrix}}, \code{\link{term_npar}},
#' \code{\link{term_coef_names}} and \code{\link{term_predict}} reject a
#' specification, and this predicate is the test they use.
#'
#' @param term An object inheriting from class \code{\link{model_term}}.
#'
#' @return A logical scalar.
#'
#' @examples
#' term_is_built(linpar(~x))
#' term_is_built(term_build(linpar(~x), data.frame(x = 1:4)))
#'
#' @seealso \code{\link{term_build}}, \code{\link{term_predict}}, \code{\link{term_refresh}}, \code{\link{term_matrix}}, \code{\link{term_coef_names}}, \code{\link{term_npar}}
#' @export
term_is_built <- function(term) {
  if (!S7::S7_inherits(term, model_term)) {
    stop("'term' must inherit from 'model_term'.", call. = FALSE)
  }
  S7::S7_inherits(term, additive_term) && length(term@coef_names) > 0L
}

.assert_built <- function(term) {
  if (!term_is_built(term)) {
    stop("the term has not been built; call term_build(term, data) first.",
         call. = FALSE)
  }
  invisible(term)
}

#' @title Design Block of a Built Term
#'
#' @description
#' The \eqn{n \times k} design block of a built additive term, with the
#' term's coefficient names as column names.
#'
#' @param term A built term (see \code{\link{term_build}}).
#' @param ... Passed to methods.
#'
#' @return A numeric matrix.
#'
#' @examples
#' term_matrix(term_build(linpar(~x), data.frame(x = 1:4)))
#'
#' @seealso \code{\link{term_build}}, \code{\link{term_predict}}, \code{\link{term_refresh}}, \code{\link{term_coef_names}}, \code{\link{term_npar}}, \code{\link{term_is_built}}
#' @export
term_matrix <- S7::new_generic("term_matrix", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_matrix, additive_term) <- function(term, ...) {
  .assert_built(term)
  term@X
}

#' @title Penalty of a Term
#'
#' @description
#' The penalty attached to the whole of the term's coefficients, or
#' \code{NULL} when there is none. The hyperparameters, their bounds and
#' links, and every derivative in the coefficients and the hyperparameters
#' are the penalty object's, not the term's.
#'
#' @details
#' A term whose penalty reaches only part of its parameters returns
#' \code{NULL} here and declares that penalty through
#' \code{\link{term_penalties}}, which names the parameters it covers:
#' \code{\link{seg}} penalizes the changes and not the linear effect or the
#' break-points, and the developments of \code{\link{nl}} and
#' \code{\link{gas}} carry their sub-terms' penalties. Reading
#' a partial penalty here would say that it covers the block, so the
#' question this generic asks is answered only where the answer is the whole
#' of it.
#'
#' @param term An object inheriting from class \code{\link{additive_term}}.
#' @param ... Passed to methods.
#'
#' @return A penalty object, or \code{NULL}.
#'
#' @examples
#' term_penalty(linpar(~x))
#'
#' @seealso \code{\link{term_penalties}}, \code{\link{term_smooth}}, \code{\link{edf}}
#' @export
term_penalty <- S7::new_generic("term_penalty", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_penalty, additive_term) <- function(term, ...) {
  term@penalty
}

#' @title Every Penalty a Term Carries
#'
#' @description
#' What a term declares it wants penalized: a list of entries, each naming a
#' subset of the term's own parameters and the penalty over them.
#'
#' @details
#' \code{\link{term_penalty}} answers for the common case, one penalty over the
#' whole of a term's design block, and this generalizes it in two directions a
#' model layer needs.
#'
#' A term may carry \strong{more than one} penalty, over different parameters of
#' its own. A panel model with a population value and a departure per group
#' wants the population value free and the departures shrunk, which is one
#' penalty over part of the parameters and none over the rest.
#'
#' The parameters need \strong{not be coefficients of a design block}. The
#' persistence of a score-driven term, the nonlinear parameters of
#' \code{\link{nl}}, the break-point of \code{\link{seg}} are parameters of the
#' term and nothing else, and everything a penalty needs from them is a vector
#' of numbers and their positions.
#'
#' The base method answers from \code{term_penalty()}, so a term that carries
#' one penalty over its whole block -- every term shipped here -- needs no
#' method of its own and behaves exactly as before. Its single entry is named
#' with the empty string, meaning the whole term, so a caller that keys the
#' hyperparameters by term name keys them exactly as it did.
#'
#' @param term A built term.
#' @param ... Passed to methods.
#'
#' @return A list, possibly empty. Each entry has \code{name} (a label unique
#'   WITHIN the term, empty for a penalty over the whole of it), \code{index}
#'   (positions among the term's parameters) and \code{penalty} (a
#'   \pkg{penalties7} object). The name is not the term's: two \code{ridge()}
#'   terms in one formula are two terms with their own hyperparameters, and it
#'   is the caller that knows what it called each one.
#'
#' @examples
#' term_penalties(term_build(ridge(~x), data.frame(x = rnorm(20))))
#' term_penalties(term_build(linpar(~x), data.frame(x = rnorm(20))))
#'
#' @seealso \code{\link{term_penalty}}, \code{\link{term_npar}}
#' @export
term_penalties <- S7::new_generic("term_penalties", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_penalties, model_term) <- function(term, ...) {
  pen <- tryCatch(term_penalty(term), error = function(e) NULL)
  if (is.null(pen)) return(list())
  list(list(name = "", index = seq_len(term_npar(term)), penalty = pen))
}

#' @title Number of Parameters of a Built Term
#'
#' @description
#' How many parameters of its own a built term carries: the columns of the
#' design block for an additive term, and the entries of
#' \code{\link{term_params}} for a structural one, which contributes no
#' block. It is the length of the vector \code{\link{term_penalties}}
#' indexes into.
#'
#' @param term A built term (see \code{\link{term_build}}).
#' @param ... Passed to methods.
#'
#' @return An integer.
#'
#' @examples
#' term_npar(term_build(linpar(~x), data.frame(x = 1:4)))
#'
#' @seealso \code{\link{term_build}}, \code{\link{term_predict}}, \code{\link{term_refresh}}, \code{\link{term_matrix}}, \code{\link{term_coef_names}}, \code{\link{term_is_built}}
#' @export
term_npar <- S7::new_generic("term_npar", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_npar, additive_term) <- function(term, ...) {
  .assert_built(term)
  ncol(term@X)
}

S7::method(term_npar, structural_term) <- function(term, ...) {
  length(term_params(term))
}

#' @title Coefficient Names of a Built Term
#'
#' @description
#' The names of the term's coefficients, prefixed by the term's label when
#' the label is non-empty.
#'
#' @param term A built term (see \code{\link{term_build}}).
#' @param ... Passed to methods.
#'
#' @return A character vector.
#'
#' @examples
#' term_coef_names(term_build(linpar(~x), data.frame(x = 1:4)))
#'
#' @seealso \code{\link{term_build}}, \code{\link{term_predict}}, \code{\link{term_refresh}}, \code{\link{term_matrix}}, \code{\link{term_npar}}, \code{\link{term_is_built}}
#' @export
term_coef_names <- S7::new_generic("term_coef_names", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_coef_names, additive_term) <- function(term, ...) {
  .assert_built(term)
  term@coef_names
}

#' @title Whether a Term's Penalized Objective Is Smooth
#'
#' @description
#' \code{TRUE} when the term's contribution to the penalized objective is
#' differentiable in the coefficients. The answer is read from the penalties
#' rather than declared by the term: an unpenalized term is smooth, and a
#' penalized one is smooth exactly when no penalty it carries declares a
#' kink, so a term cannot disagree with its own penalties. The model layer
#' uses this flag to split the coefficient vector into the block the
#' classical optimizers handle and the block that needs non-smooth
#' strategies.
#'
#' @details
#' The enumeration is \code{\link{term_penalties}}, so a term carrying one
#' penalty over part of its parameters and none over the rest answers for
#' the part: \code{seg(x, penalty = penalties7::lasso_penalty)} is not
#' smooth, its slope
#' changes sitting at a kink, although its linear effect and its
#' break-points are unpenalized.
#'
#' @param term An object inheriting from class \code{\link{model_term}}.
#' @param ... Passed to methods.
#'
#' @return A logical scalar.
#'
#' @examples
#' term_smooth(linpar(~x))
#'
#' @seealso \code{\link{term_penalties}}, \code{\link{term_penalty}}, \code{\link{edf}}
#' @export
term_smooth <- S7::new_generic("term_smooth", "term",
  function(term, ...) S7::S7_dispatch())

# The penalty a term's `penalty` argument asks for, over a given number of
# coordinates: no map, so that the separable branch of penalties7 applies and
# a fitting layer keeps its proximal step and its coordinate descent.
#
# Two spellings. A penalties7 penalty is used as it stands, so anything that
# package offers -- an elastic net, a heavy-tailed prior, a structured
# precision -- reaches a term without a name having to be invented for it
# here; and a function of the coordinate count is what a penalty whose SIZE
# is not known until the data are seen needs, a panel's deviations existing
# only once the groups are counted. A function is called by the NAME
# `n_coef` when it has that formal -- which is what lets a penalties7
# constructor be passed bare, `penalty = penalties7::ridge_penalty`, whose
# FIRST formal is the map and would otherwise receive the count --
# and positionally otherwise, which is what `function(k) ...` wants.
.penalty_factory <- function(kind) {
  if (is.function(kind)) {
    by_name <- "n_coef" %in% names(formals(kind))
    return(function(k) {
      pen <- if (by_name) kind(n_coef = k) else kind(k)
      .penalty_check(pen, k, "the function")
    })
  }
  if (S7::S7_inherits(kind, penalties7::penalty)) {
    return(function(k) .penalty_check(kind, k, "the penalty"))
  }
  stop(sprintf("unknown penalty '%s'.", as.character(kind)[1L]),
       call. = FALSE)
}

# what a supplied penalty has to be, checked where the coordinate count is
# finally known: a penalty of the wrong width would be evaluated at a
# coefficient vector of another length and R would recycle it in silence
.penalty_check <- function(pen, k, what) {
  if (!S7::S7_inherits(pen, penalties7::penalty)) {
    stop(sprintf("%s must give a penalties7 penalty; it gave a %s.",
                 what, paste(class(pen), collapse = "/")), call. = FALSE)
  }
  if (!identical(as.integer(pen@n_coef), as.integer(k))) {
    stop(sprintf(paste("%s covers %d coefficients and the term has %d.",
                       "Pass a function of the count where the term's size",
                       "is not known in advance."),
                 what, as.integer(pen@n_coef), as.integer(k)), call. = FALSE)
  }
  pen
}

# is a `penalty` argument the absence of one? NULL is the default of every
# term that takes the argument
.penalty_is_none <- function(kind) is.null(kind)

# What a `penalty` argument may be, checked in the CONSTRUCTOR so that a
# mistake is reported where it is written rather than at term_build(): a
# penalties7 penalty, or a function of the number of coefficients returning
# one. A string is rejected -- a name was a second spelling of an object the
# toolkit already has, and `penalty = penalties7::ridge_penalty` is exactly
# as short.
.penalty_arg <- function(kind) {
  if (is.null(kind) || is.function(kind) ||
      S7::S7_inherits(kind, penalties7::penalty)) {
    return(kind)
  }
  stop(paste("'penalty' must be NULL, a penalties7 penalty, or a function",
             "of the number of coefficients returning one",
             "(e.g. penalty = penalties7::lasso_penalty)."), call. = FALSE)
}

# how a penalty argument prints
.penalty_label <- function(kind) {
  if (is.function(kind)) return("a penalty per coefficient count")
  if (S7::S7_inherits(kind, penalties7::penalty)) return(kind@penalty_name)
  as.character(kind)[1L]
}

# a hyperparameter value inside each domain, at which the kink set is asked
# for; the kinks are structural, so any admissible value answers the question
.penalty_probe_theta <- function(pen) {
  stats::setNames(lapply(pen@params, function(p) {
    b <- pen@params_bounds[[p]]
    if (is.finite(b[1L]) && is.finite(b[2L])) return(mean(b))
    if (is.finite(b[1L])) return(b[1L] + 1)
    if (is.finite(b[2L])) return(b[2L] - 1)
    0
  }), pen@params)
}

S7::method(term_smooth, model_term) <- function(term, ...) {
  for (e in term_penalties(term)) {
    pen <- e$penalty
    if (length(penalties7::penalty_kinks(pen, .penalty_probe_theta(pen)))) {
      return(FALSE)
    }
  }
  TRUE
}

#' @title Design Block on New Data
#'
#' @description
#' Applies a built term's mapping to new data, reproducing the block the
#' term would have produced had the new rows been part of the original
#' data: factor levels, contrasts and any constants recorded in the
#' blueprint at build time are reused, never recomputed. New data carrying
#' a factor level unknown to the blueprint is rejected.
#'
#' @details
#' The block returned is \eqn{\tilde{X}_t} such that
#' \eqn{\tilde{\eta} = \tilde{X}_t \beta_t} is the term's contribution at
#' the new rows, evaluated at the coefficients the model already carries.
#' Reapplying the recorded mapping rather than rebuilding it is what makes
#' that identity hold: a rebuilt factor encoding, spline knot placement or
#' basis reparametrization would give a block of the same shape multiplying
#' the same coefficients and meaning something else.
#'
#' @param term A built term (see \code{\link{term_build}}).
#' @param newdata A data frame.
#' @param ... Passed to methods.
#'
#' @return A numeric matrix with \code{nrow(newdata)} rows and one column
#'   per coefficient.
#'
#' @examples
#' built <- term_build(linpar(~x), data.frame(x = 1:4))
#' term_predict(built, data.frame(x = c(0.5, 2.5)))
#'
#' @seealso \code{\link{term_build}}, \code{\link{term_refresh}}, \code{\link{term_matrix}}, \code{\link{term_coef_names}}, \code{\link{term_npar}}, \code{\link{term_is_built}}
#' @export
term_predict <- S7::new_generic("term_predict", "term",
  function(term, newdata, ...) {
    if (!is.data.frame(newdata)) {
      stop("'newdata' must be a data frame.", call. = FALSE)
    }
    S7::S7_dispatch()
  })

# --- printing ---------------------------------------------------------------

S7::method(print, model_term) <- function(x, ...) {
  cls <- attr(S7::S7_class(x), "name")
  lab <- if (nzchar(x@label)) sprintf(" '%s'", x@label) else ""
  if (term_is_built(x)) {
    cat(sprintf("<%s>%s built: %d coefficient%s\n", cls, lab,
                ncol(x@X), if (ncol(x@X) == 1L) "" else "s"))
  } else {
    cat(sprintf("<%s>%s (specification; call term_build() with data)\n",
                cls, lab))
  }
  invisible(x)
}
