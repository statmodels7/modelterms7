#' @include term_classes.R generics.R
NULL

# Which of a term's hyperparameters the caller has held, and at what.
#
# The choice belongs to the term because the term is where the penalty is
# named. A criterion decides HOW an estimated hyperparameter is found; it
# does not decide WHICH ones are estimated, and an argument on the criterion
# saying so was read by nothing when the term disagreed with it.

#' Check a Term's Held Hyperparameters Against Its Penalty
#'
#' @description
#' Validates the values a constructor was given against the names and the
#' bounds of the penalty the term will build, and returns them as a named
#' list with the \code{NULL} entries -- the ones to be estimated -- dropped.
#'
#' @details
#' The check happens at CONSTRUCTION, where the caller can see it, and not at
#' the fit three layers away. A name the penalty does not carry is an error
#' naming the ones it does, which is what turns \code{ridge(lambda = 2)} --
#' the elastic net's spelling applied to a Gaussian prior, whose
#' hyperparameter is \code{sigma} -- into a message instead of an argument
#' that lands in \code{...} and does nothing.
#'
#' The bounds are OPEN, as they are everywhere in the toolkit: a Gaussian
#' prior at \eqn{\sigma = 0} is not a prior, and an elastic net at
#' \eqn{\alpha = 0} has no kink and is a penalty of another kind.
#'
#' @param values A named list of the constructor's arguments, \code{NULL}
#'   where the hyperparameter is to be estimated.
#' @param penalty A \pkg{penalties7} penalty, or a function returning one,
#'   used only to read the names and the bounds.
#' @param what The constructor's name, for the message.
#'
#' @return A named list of the held values.
#'
#' @seealso \code{\link{term_hyper}}, \code{\link{term_penalties}}
#'
#' @keywords internal
check_hyper <- function(values, penalty, what = "this term") {
  values <- values[!vapply(values, is.null, logical(1))]
  if (!length(values)) return(list())
  pen <- if (is.function(penalty)) {
    tryCatch(penalty(1L), error = function(e) NULL)
  } else {
    penalty
  }
  if (is.null(pen)) return(as.list(values))
  ok <- pen@params
  for (nm in names(values)) {
    if (!nm %in% ok) {
      stop(sprintf(paste0("'%s' has no hyperparameter '%s'. It carries: %s.",
                          "\n  Give one of those to hold it, or leave it",
                          " NULL to have it estimated."),
                   what, nm, paste(ok, collapse = ", ")), call. = FALSE)
    }
    v <- values[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v)) {
      stop(sprintf("'%s' in '%s' must be a single finite number, or NULL.",
                   nm, what), call. = FALSE)
    }
    b <- pen@params_bounds[[nm]]
    if (!is.null(b) && (v <= b[1L] || v >= b[2L])) {
      stop(sprintf(paste0("'%s' in '%s' must lie strictly inside (%s, %s);",
                          " it is %s."),
                   nm, what, format(b[1L]), format(b[2L]), format(v)),
           call. = FALSE)
    }
  }
  as.list(values)
}


#' @title The Hyperparameters a Term Holds
#'
#' @description
#' The values the caller fixed in the constructor, by hyperparameter. Those a
#' term does not name are estimated by whatever criterion the fit runs.
#'
#' @details
#' A term carrying several penalties answers per penalty, under the same names
#' \code{\link{term_penalties}} gives its entries, and every entry of that
#' enumeration carries its own held values in the field \code{fixed} -- so a
#' term that copies the entries of its sub-terms, which is what a structural
#' term with subformulas does, propagates them without knowing they exist.
#'
#' @param term A term, built or not.
#' @param ... Passed to methods.
#'
#' @return A named list, one entry per penalty of the term, each a named list
#'   of held values. Empty where the term holds nothing.
#'
#' @examples
#' term_hyper(lasso(~x, lambda = 3))
#' term_hyper(lasso(~x))
#'
#' @seealso \code{\link{term_penalties}}, \code{\link{ridge}}
#' @export
term_hyper <- S7::new_generic("term_hyper", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_hyper, model_term) <- function(term, ...) {
  # A BUILT term answers from its entries, so a term whose penalties come
  # from sub-terms -- a structural one with subformulas -- reports what its
  # sub-terms hold without a method of its own. An unbuilt one has no
  # entries yet and answers from its own property.
  ent <- tryCatch(term_penalties(term), error = function(e) list())
  if (length(ent)) {
    out <- list()
    for (e in ent) {
      f <- e$fixed
      if (is.null(f) || !length(f)) next
      out[[if (is.null(e$name)) "" else e$name]] <- f
    }
    if (length(out)) return(out)
  }
  h <- term@hyper
  if (!length(h)) return(list())
  stats::setNames(list(h), "")
}


#' The Held Smoothing Parameters of a Smooth
#'
#' @description
#' Carries the \code{lambda} a smooth's constructor was given onto the names
#' its penalty uses, and checks it against the positivity every smoothing
#' parameter obeys.
#'
#' @details
#' A one-dimensional smooth and an isotropic tensor product carry one
#' smoothing parameter, named \code{lambda}; an anisotropic tensor product
#' carries one per margin, \code{lambda1}, \code{lambda2} and so on, which is
#' what \code{\link[penalties7]{additive_penalty}} names them. An unnamed
#' vector is read in that order and must be as long as there are margins; a
#' named one may hold some and leave the others to be estimated.
#'
#' @param lambda What the constructor was given, or \code{NULL}.
#' @param names The penalty's own hyperparameter names.
#' @param what The term's label, for the message.
#'
#' @return A named list of held values.
#'
#' @seealso \code{\link{check_hyper}}, \code{\link{s}}, \code{\link{te}}
#'
#' @keywords internal
smooth_hyper <- function(lambda, names, what = "this smooth") {
  if (is.null(lambda)) return(list())
  if (!is.numeric(lambda) || !length(lambda) || anyNA(lambda)) {
    stop(sprintf("'lambda' in '%s' must be a number, or NULL.", what),
         call. = FALSE)
  }
  nm <- base::names(lambda)
  if (is.null(nm)) {
    if (length(lambda) == 1L && length(names) > 1L) {
      # one number for a term with a smoothing parameter per margin is the
      # isotropic reading of an anisotropic penalty, which is a different
      # model rather than a shorthand for this one
      stop(sprintf(paste0("'%s' has %d smoothing parameters (%s), and",
                          " 'lambda' has one.\n  Give one per margin, or",
                          " name the ones to hold, or use anisotropic =",
                          " FALSE."),
                   what, length(names), paste(names, collapse = ", ")),
           call. = FALSE)
    }
    if (length(lambda) != length(names)) {
      stop(sprintf("'lambda' in '%s' must have %d value%s.", what,
                   length(names), if (length(names) == 1L) "" else "s"),
           call. = FALSE)
    }
    lambda <- stats::setNames(lambda, names)
  }
  bad <- setdiff(base::names(lambda), names)
  if (length(bad)) {
    stop(sprintf(paste0("'%s' has no smoothing parameter '%s'. It carries:",
                        " %s."), what, bad[1L], paste(names, collapse = ", ")),
         call. = FALSE)
  }
  if (any(lambda <= 0)) {
    stop(sprintf("every 'lambda' in '%s' must be strictly positive.", what),
         call. = FALSE)
  }
  as.list(lambda)
}


#' Read a Caller's Held Hyperparameters
#'
#' @description
#' Accepts a named vector or a named list and returns a named list, checking
#' only what can be checked before the penalty exists.
#'
#' @details
#' \code{\link{random}} builds one of three penalties depending on what it was
#' given -- a ridge, a structured prior over a \pkg{parameters7} matrix, or a
#' \pkg{distributions7} family used coordinatewise -- so which names there are
#' is not known until the term is built. The shape is checked here and the
#' names against the penalty by \code{\link{check_hyper}} at that point.
#'
#' @param x A named vector, a named list, or \code{NULL}.
#' @param what The term's label, for the message.
#'
#' @return A named list, possibly empty.
#'
#' @keywords internal
as_hyper <- function(x, what = "this term") {
  if (is.null(x) || !length(x)) return(list())
  nm <- names(x)
  if (is.null(nm) || !all(nzchar(nm))) {
    stop(sprintf(paste0("'hyper' in '%s' must be named, one entry per",
                        " hyperparameter to hold."), what), call. = FALSE)
  }
  x <- as.list(x)
  for (h in nm) {
    v <- x[[h]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v)) {
      stop(sprintf("'hyper$%s' in '%s' must be a single finite number.",
                   h, what), call. = FALSE)
    }
  }
  x
}
