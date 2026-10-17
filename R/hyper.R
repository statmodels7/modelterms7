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
#' Validates the values that a constructor was given against the names and the
#' bounds of the penalty that the term will build, and returns them as a named
#' list. The `NULL` entries, which are the hyperparameters to be estimated, are
#' dropped.
#'
#' @details
#' The check is made at construction and not at the fit. A name that the
#' penalty does not carry is an error that lists the names that it does carry.
#' In the penalized constructors, an argument written for another penalty, such
#' as `mcp(x, a = 3)`, which gives SCAD's shape to an MCP whose own shape is
#' `gamma`, is rejected earlier, as an argument that the constructor does not
#' have (see [penalized_terms()]).
#'
#' The bounds are open, as they are everywhere in the toolkit: a ridge at
#' \eqn{\lambda = 0} is no penalty at all, and an elastic net at
#' \eqn{\alpha = 0} has no kink and is a penalty of another kind.
#'
#' One argument carries three states, read for one hyperparameter at a time:
#' `NULL` leaves the hyperparameter to be estimated, one number holds it, and
#' several numbers are the grid itself. This returns the second of the three;
#' [check_values()] returns the third.
#'
#' @param values A named list of the constructor's arguments, `NULL`
#'   where the hyperparameter is to be estimated.
#' @param penalty A \pkg{penalties7} penalty, or a function returning one,
#'   used only to read the names and the bounds.
#' @param what The constructor's name, for the message.
#'
#' @return A named list, keyed by the penalty's own hyperparameter names, that
#'   holds the single value of each hyperparameter that was given exactly one
#'   number.
#'   The list is empty where every hyperparameter was left `NULL` or written
#'   out as a grid.
#'
#' @seealso [check_values()], [term_hyper()],
#'   [term_penalties()]
#'
#' @keywords internal
check_hyper <- function(values, penalty, what = "this term") {
  .hyper_parts(values, penalty, what)$hyper
}


#' Check a Term's Written-Out Grids Against Its Penalty
#'
#' @description
#' Validates the hyperparameter values that a constructor was given as a
#' vector, which a path visits as they stand.
#'
#' @details
#' Several values are neither a held hyperparameter nor a request to build a
#' grid: they are the grid. Each is checked against the penalty's bounds as a
#' held value would be, and the bounds are the only check. The value that
#' empties the block does not cap them and the depth of the path does not
#' apply to them, since the values to visit are given.
#'
#' They are sorted, because a path is walked from the emptiest fit towards the
#' fullest and its warm starts follow that order, and duplicates are dropped.
#' Which end of the order is the sparse one depends on the penalty, so the
#' direction is settled where the path is built.
#'
#' @inheritParams check_hyper
#'
#' @return A named list of numeric vectors, one entry per hyperparameter that
#'   the caller wrote out.
#'
#' @seealso [check_hyper()], [term_values()]
#'
#' @keywords internal
check_values <- function(values, penalty, what = "this term") {
  .hyper_parts(values, penalty, what)$values
}


#' Split a Constructor's Hyperparameter Arguments by Length
#'
#' @description
#' Validates each against the penalty and files it as held or as a grid.
#'
#' @details
#' The validation is one body because the two states differ only in length:
#' the name must be one that the penalty carries and every value must lie
#' strictly inside its bounds, whether there is one of them or twenty.
#'
#' @inheritParams check_hyper
#'
#' @return A list of two named lists: `hyper`, the hyperparameters given one
#'   value each, and `values`, those given several, sorted and deduplicated.
#'   Both are keyed by the penalty's own hyperparameter names, and either may
#'   be empty.
#'
#' @keywords internal
.hyper_parts <- function(values, penalty, what = "this term") {
  none <- list(hyper = list(), values = list())
  values <- values[!vapply(values, is.null, logical(1))]
  if (!length(values)) return(none)
  pen <- if (is.function(penalty)) {
    tryCatch(penalty(1L), error = function(e) NULL)
  } else {
    penalty
  }
  held <- list()
  grids <- list()
  for (nm in names(values)) {
    v <- values[[nm]]
    if (!is.null(pen) && !nm %in% pen@params) {
      stop(sprintf(paste0("'%s' has no hyperparameter '%s'. It carries: %s.",
                          "\n  Give one of those to hold it, or leave it",
                          " NULL to have it estimated."),
                   what, nm, paste(pen@params, collapse = ", ")),
           call. = FALSE)
    }
    if (!is.numeric(v) || !length(v) || any(!is.finite(v))) {
      stop(sprintf(paste0("'%s' in '%s' must be a number to hold it, several",
                          " numbers to sweep\n  exactly those, or NULL to",
                          " have the grid built."), nm, what), call. = FALSE)
    }
    b <- if (is.null(pen)) NULL else pen@params_bounds[[nm]]
    if (!is.null(b) && any(v <= b[1L] | v >= b[2L])) {
      bad <- v[v <= b[1L] | v >= b[2L]][[1L]]
      stop(sprintf(paste0("'%s' in '%s' must lie strictly inside (%s, %s);",
                          " it is %s."),
                   nm, what, format(b[1L]), format(b[2L]), format(bad)),
           call. = FALSE)
    }
    if (length(v) == 1L) {
      held[[nm]] <- as.numeric(v)
    } else {
      grids[[nm]] <- sort(unique(as.numeric(v)))
    }
  }
  list(hyper = held, values = grids)
}


#' @title The Hyperparameters That a Term Holds
#'
#' @description
#' Reports the hyperparameter values that the caller fixed in the constructor,
#' one entry per penalty that the term carries. A hyperparameter left `NULL`,
#' or written out as a grid, is estimated by the fit and does not appear
#' here, so an empty result means that the term holds none.
#'
#' @details
#' # The three states of one argument
#'
#' A constructor's hyperparameter argument carries three meanings, read for
#' one hyperparameter at a time. `NULL` leaves the hyperparameter to be
#' estimated: a path builds a grid for a penalty with a kink, and a marginal
#' criterion estimates the value at the mode for a smooth one. `NULL` is the
#' default for `lambda` and `alpha`, whereas the shape of [scad()] and [mcp()]
#' defaults to a held number, and `NULL` there has the shape estimated.
#' One number holds the hyperparameter, and this function reports that number.
#' Several numbers are a grid that the path visits as they stand, and
#' [term_values()] reports that grid. The hyperparameter is still estimated,
#' because the values to search are given and the estimate is left to the fit.
#'
#' # How the entries are keyed
#'
#' The names are the ones [term_penalties()] gives its entries: `""` for a
#' penalty covering the whole term, and `parameter::label` for a penalty
#' brought in by a subformula. A term whose penalties come from sub-terms, such
#' as [nl()] or [seg()] with subformulas, propagates their held values without
#' a method of its own, because every entry of that enumeration carries its own
#' in the field `fixed`.
#'
#' A built term reports those entries, and an unbuilt term reports its own
#' `hyper` property under the names as they were written. The two answers are
#' the same for the penalized terms, and for [s()] and [random()] where the
#' names are not qualified. A qualified name, such as `lambda.a` for the level
#' `a` of a factor, is reported as written before the build and under its entry
#' after it. For a term whose held values sit in subformulas, such as [nl()] and
#' [seg()], the sub-terms exist only after the build, so an unbuilt term
#' reports an empty list.
#'
#' @param term A term, built or not.
#' @param ... Passed to methods. The shipped methods do not use it.
#'
#' @return A named list, one entry per penalty of the term and named as
#'   [term_penalties()] names it, each entry a named list of single numbers
#'   keyed by the penalty's own hyperparameter names. The list is empty where
#'   the term holds nothing.
#'
#' @seealso [term_values()] for a written-out grid, [term_grid()] for how many
#'   values a built grid visits, [term_penalties()] for the entry names that
#'   key the result.
#'
#' @examples
#' # The first term holds its hyperparameter and the second leaves it free.
#' term_hyper(lasso(~ x, lambda = 3))
#' term_hyper(lasso(~ x))
#'
#' # Several values are a grid, not a held value, so they go elsewhere.
#' term_hyper(lasso(~ x, lambda = c(0.1, 1, 10)))
#' term_values(lasso(~ x, lambda = c(0.1, 1, 10)))
#'
#' # One of the two hyperparameters is held, and the other is estimated.
#' term_hyper(enet(~ x, alpha = 0.5))
#'
#' # The shape of SCAD is held by default; NULL has it estimated.
#' term_hyper(scad(~ x))
#' term_hyper(scad(~ x, a = NULL))
#'
#' # A held value on a sub-term is reported under that entry's name once the
#' # term is built.
#' set.seed(4)
#' d <- data.frame(x = runif(40, 0, 5), g = factor(rep(1:4, 10)))
#' nu <- nl(~ a * exp(-r * x), a ~ 0 + lasso(~ g, lambda = 2), r ~ 1)
#' term_hyper(nu)
#' term_hyper(term_build(nu, d))
#'
#' @export
#' @aliases term_hyper.model_term
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
#' Carries the `hyper` argument of [s()] or [te()] onto the names that its
#' penalty uses, and checks that every held smoothing parameter is finite and
#' strictly positive.
#'
#' @details
#' A one-dimensional smooth and an isotropic tensor product carry one
#' smoothing parameter, named `lambda`; an anisotropic tensor product
#' carries one per margin, `lambda1`, `lambda2` and so on, as
#' [penalties7::additive_penalty()] names them. An unnamed
#' vector is read in that order and must be as long as there are margins; a
#' named one may hold some and leave the others to be estimated.
#'
#' @param hyper What the constructor was given, or `NULL`.
#' @param names The penalty's own hyperparameter names.
#' @param what The term's label, for the message.
#'
#' @return A named list of held values.
#'
#' @seealso [check_hyper()], [s()], [te()]
#'
#' @keywords internal
smooth_hyper <- function(hyper, names, what = "this smooth") {
  if (is.null(hyper)) return(list())
  if (!is.numeric(hyper) || !length(hyper) || !all(is.finite(hyper))) {
    stop(sprintf("'hyper' in '%s' must be a finite number, or NULL.", what),
         call. = FALSE)
  }
  nm <- base::names(hyper)
  if (is.null(nm)) {
    if (length(hyper) == 1L && length(names) > 1L) {
      # one number for a term with a smoothing parameter per margin is the
      # isotropic reading of an anisotropic penalty, which is a different
      # model rather than a shorthand for this one
      stop(sprintf(paste0("'%s' has %d smoothing parameters (%s), and",
                          " 'hyper' has one.\n  Give one per margin, or",
                          " name the ones to hold, or use anisotropic =",
                          " FALSE."),
                   what, length(names), paste(names, collapse = ", ")),
           call. = FALSE)
    }
    if (length(hyper) != length(names)) {
      stop(sprintf("'hyper' in '%s' must have %d value%s.", what,
                   length(names), if (length(names) == 1L) "" else "s"),
           call. = FALSE)
    }
    hyper <- stats::setNames(hyper, names)
  }
  bad <- setdiff(base::names(hyper), names)
  if (length(bad)) {
    stop(sprintf(paste0("'%s' has no smoothing parameter '%s'. It carries:",
                        " %s."), what, bad[1L], paste(names, collapse = ", ")),
         call. = FALSE)
  }
  if (any(hyper <= 0)) {
    stop(sprintf("every 'hyper' in '%s' must be strictly positive.", what),
         call. = FALSE)
  }
  as.list(hyper)
}


#' Read a Caller's Held Hyperparameters
#'
#' @description
#' Accepts a named vector or a named list and returns a named list, checking
#' only what can be checked before the penalty exists.
#'
#' @details
#' [random()] builds its penalty from the distribution of the effects with
#' [penalties7::distrib_penalty()]: a centered Gaussian by default, a
#' multivariate \pkg{distributions7} family for correlated effects, or a
#' univariate family used once for each within-group column. A smooth whose
#' penalty comes from a factory, or that has one smoothing parameter per level
#' of a factor, is read in the same way. The names of the hyperparameters are
#' therefore known only when the term is built. The shape is checked here and
#' the names against the penalty by [check_hyper()] at that point.
#'
#' A vector entry is a written-out grid, exactly as it is in the constructors
#' that name their hyperparameters, so the length is not checked here either;
#' [.hyper_parts()] splits the two once the penalty exists.
#'
#' @param x A named vector, a named list, or `NULL`.
#' @param what The term's label, for the message.
#'
#' @return A named list, one entry per hyperparameter named, each a numeric
#'   vector of length one or more. The list is empty for `NULL` or a
#'   zero-length input.
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
    if (!is.numeric(v) || !length(v) || any(!is.finite(v))) {
      stop(sprintf(paste0("'hyper$%s' in '%s' must be a number to hold it,",
                          " or several to sweep\n  exactly those."),
                   h, what), call. = FALSE)
    }
  }
  x
}


#' Check a Term's Grid Sizes Against Its Penalty
#'
#' @description
#' Validates the number of values that the constructor was given per
#' hyperparameter, and returns them as a named list from which the `NULL`
#' entries are dropped.
#'
#' @details
#' How finely a hyperparameter is swept belongs to the term for the same
#' reason as whether it is swept at all: a penalized block of four columns
#' and one of four hundred need different grids, and a criterion is
#' applied to every term of the model at once. Where a term gives no size,
#' the path of the fitting layer uses its own default.
#'
#' The smallest grid has two values. There is no upper limit; each point of a
#' path is a whole fit, so the cost grows linearly in the number of values on
#' one axis.
#'
#' @param values A named list of the constructor's arguments, `NULL`
#'   where the path's default size is wanted.
#' @param penalty A \pkg{penalties7} penalty, or a function returning one,
#'   used only to read the names.
#' @param what The constructor's name, for the message.
#'
#' @return A named list of grid sizes.
#'
#' @seealso [check_hyper()], [term_grid()]
#'
#' @keywords internal
check_grid <- function(values, penalty, what = "this term") {
  values <- values[!vapply(values, is.null, logical(1))]
  if (!length(values)) return(list())
  pen <- if (is.function(penalty)) {
    tryCatch(penalty(1L), error = function(e) NULL)
  } else {
    penalty
  }
  for (nm in names(values)) {
    if (!is.null(pen) && !nm %in% pen@params) {
      stop(sprintf(paste0("'%s' has no hyperparameter '%s' on which to",
                          " put a grid. It carries: %s."),
                   what, nm, paste(pen@params, collapse = ", ")),
           call. = FALSE)
    }
    v <- values[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) ||
        v != round(v) || v < 2) {
      stop(sprintf(paste0("the grid for '%s' in '%s' must be a whole number",
                          " of at least 2, or NULL."), nm, what),
           call. = FALSE)
    }
    values[[nm]] <- as.integer(v)
  }
  as.list(values)
}


#' @title The Grid Sizes of a Term
#'
#' @description
#' Reports the number of values that a path visits for each hyperparameter of
#' the term that is swept, one entry per penalty. A hyperparameter that the
#' term holds has no grid and is not reported. A swept hyperparameter for
#' which the term records no size is covered at the default size of the fitting
#' layer.
#'
#' @details
#' How finely a hyperparameter is swept belongs to the term because the term is
#' where the penalty is named: a block of four columns and one of four
#' hundred need different grids, and an outer criterion applies to every
#' hyperparameter of the model alike, so a size set there could not
#' distinguish between them.
#'
#' Only a hyperparameter with a path has a grid, and a held
#' hyperparameter has none. [ridge()], [s()] and [te()] with their default
#' penalty, and [random()] with its default Gaussian prior, report nothing,
#' because a marginal criterion estimates their hyperparameters at the mode.
#' A smoother whose penalty factory has a kink is swept by a path but reports
#' nothing either, and the path uses its default size. [lasso()] reports its
#' `n_lambda`, [enet()] reports `n_lambda` and `n_alpha`, and [scad()] and
#' [mcp()] report `n_lambda` together with `n_a` or `n_gamma` where the shape
#' is estimated, since the shape is held by default. These constructors write
#' their sizes into the term, so `lasso(~ x)` reports `lambda = 25` where an
#' unset argument would give an empty list. A hyperparameter written out as a
#' vector has no grid to build and is not reported here; the path visits the
#' written-out values, which [term_values()] reports.
#'
#' The keys are [term_penalties()]'s entry names, `""` for a penalty over the
#' whole term.
#'
#' @param term A term, built or not.
#' @param ... Passed to methods. The shipped methods do not use it.
#'
#' @return A named list, one entry per penalty of the term, each a named list
#'   of single whole numbers keyed by hyperparameter. The list is empty where
#'   the term records no grid size for any hyperparameter.
#'
#' @seealso [term_path_min()] for how far down the path reaches,
#'   [term_search()] for how several of them are combined, [term_values()] for
#'   a grid written out instead of counted.
#'
#' @examples
#' # A path term reports its grid, at the default and when set.
#' term_grid(lasso(~ x))
#' term_grid(lasso(~ x, n_lambda = 50))
#'
#' # The two axes have two sizes, and the kink axis is swept more finely.
#' term_grid(enet(~ x))
#'
#' # A held hyperparameter has no grid, and the held shape of SCAD has none.
#' term_grid(enet(~ x, alpha = 0.5))
#' term_grid(scad(~ x))
#' term_grid(scad(~ x, a = NULL))
#'
#' # A penalty with no kink has no path, so there is nothing to count.
#' term_grid(ridge(~ x))
#' term_grid(s(x, basis7::bspline_smooth(k = 5)))
#'
#' @export
#' @aliases term_grid.model_term
term_grid <- S7::new_generic("term_grid", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_grid, model_term) <- function(term, ...) {
  out <- list()
  for (e in .hyper_entries(term)) {
    g <- e$n_values[.built(e)]
    if (length(g)) out[[e$name]] <- g
  }
  out
}

# The entries the hyperparameter reporters read: term_penalties() of a built
# term, and otherwise one entry over the whole term from its own properties.
.hyper_entries <- function(term) {
  ent <- tryCatch(term_penalties(term), error = function(e) list())
  if (!length(ent)) {
    ent <- list(list(name = "", fixed = term@hyper, n_values = term@grid,
                     values = term@values, min_ratio = term@min_ratio,
                     search = term@search))
  }
  lapply(ent, function(e) {
    if (is.null(e$name)) e$name <- ""
    e
  })
}

# The hyperparameters of an entry that a path sweeps: those with a grid to
# build or one written out, less those held.
.swept <- function(e) {
  setdiff(union(names(e$n_values), names(e$values)), names(e$fixed))
}

# The hyperparameters of an entry whose grid is built: swept, and not
# written out.
.built <- function(e) setdiff(.swept(e), names(e$values))


#' Reject a Written-Out Grid Where There Is No Path
#'
#' @description
#' Signals an error when the caller gave several values for a hyperparameter
#' whose penalty has no kink.
#'
#' @details
#' Several values are a grid for a path to visit, and only a penalty with a
#' kink is swept along one. For any other penalty the hyperparameter is
#' estimated at the mode by a marginal criterion, and a vector given there
#' would have no effect. Whether the penalty has a kink is read from the
#' penalty at a probe value of its own hyperparameters. The kink set is
#' structural, so a ridge and a random effect under a Gaussian prior are
#' covered by the same test as a family added later.
#'
#' @param values The written-out grids, as [check_values()] returns
#'   them.
#' @param pen The penalty that the term will build, or `NULL`.
#' @param what The term's label, for the message.
#'
#' @return `NULL`, invisibly.
#'
#' @seealso [check_values()], [term_values()]
#'
#' @keywords internal
reject_pathless_values <- function(values, pen, what = "this term") {
  if (!length(values) || is.null(pen)) return(invisible(NULL))
  k <- tryCatch(penalties7::penalty_kinks(pen, .penalty_probe_theta(pen)),
                error = function(e) 0)
  if (length(k)) return(invisible(NULL))
  stop(sprintf(paste0("'%s' in '%s' has several values, and this term has no",
                      " path to visit them:\n  its penalty has no kink, so",
                      " the hyperparameter is estimated by the\n  criterion at",
                      " the mode. Give one number to hold it, or NULL to",
                      " estimate it."),
               names(values)[[1L]], what), call. = FALSE)
}


#' @title The Hyperparameter Values That a Term Writes Out
#'
#' @description
#' Reports the values that a path visits, for each hyperparameter that the
#' caller wrote out as a vector instead of leaving to be built. It is the
#' third state of a constructor's hyperparameter argument, beside `NULL` for a
#' built grid and one number for a held value.
#'
#' @details
#' A written-out grid is used as it stands. The value that empties the block
#' does not cap it and [term_path_min()] does not extend it: those two
#' construct a grid, and here there is nothing to construct.
#'
#' A hyperparameter written out is still estimated. The caller supplies the
#' values to search and not the estimate, so a fit reports the hyperparameter as
#' chosen by the criterion, and [term_hyper()] does not carry it.
#'
#' The values are sorted and deduplicated at construction, because a path is
#' walked from the emptiest fit toward the fullest and its warm starts follow
#' that order. Which end of the sorted order is the sparse one depends on the
#' penalty, so the direction is settled where the path is built.
#'
#' The keys are [term_penalties()]'s entry names, `""` for a penalty over the
#' whole term.
#'
#' @param term A term, built or not.
#' @param ... Passed to methods. The shipped methods do not use it.
#'
#' @return A named list, one entry per penalty of the term, each a named list
#'   of numeric vectors keyed by hyperparameter, sorted and without
#'   duplicates. The list is empty where the term wrote nothing out.
#'
#' @seealso [term_hyper()] for a held value, [term_grid()] for a grid to be
#'   built, [term_penalties()] for the entry names that key the result.
#'
#' @examples
#' # A grid that is written out is reported as the grid to visit.
#' term_values(lasso(~ x, lambda = c(0.1, 1, 10)))
#'
#' # The values are sorted and deduplicated at construction.
#' term_values(lasso(~ x, lambda = c(10, 0.1, 1, 1)))
#'
#' # One number is a held value, and belongs to the other reporter.
#' term_values(lasso(~ x, lambda = 3))
#' term_hyper(lasso(~ x, lambda = 3))
#'
#' # A penalty with no kink has no path, so several values are rejected.
#' try(ridge(~ x, lambda = c(1, 2)))
#'
#' @export
#' @aliases term_values.model_term
term_values <- S7::new_generic("term_values", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_values, model_term) <- function(term, ...) {
  ent <- tryCatch(term_penalties(term), error = function(e) list())
  if (length(ent)) {
    out <- list()
    for (e in ent) {
      v <- e$values
      if (is.null(v) || !length(v)) next
      out[[if (is.null(e$name)) "" else e$name]] <- v
    }
    if (length(out)) return(out)
  }
  v <- term@values
  if (!length(v)) return(list())
  stats::setNames(list(v), "")
}


#' Check a Term's Path Depth
#'
#' @description
#' Validates the fraction of the emptying value to which a path descends.
#'
#' @details
#' The fraction is one number per term, because only the path over the size of
#' the kink uses it. A bounded hyperparameter is swept over its own interval,
#' and a shape that does not move the kink over a geometric grid above its
#' lower bound; a fraction of an emptying value has no meaning on either axis.
#'
#' @param v What the constructor was given, or `NULL`.
#' @param what The term's label, for the message.
#'
#' @return A numeric vector of length one, or of length zero.
#'
#' @seealso [check_grid()], [term_path_min()]
#'
#' @keywords internal
check_min_ratio <- function(v, what = "this term") {
  if (is.null(v)) return(numeric(0))
  if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0 || v >= 1) {
    stop(sprintf(paste0("'min_ratio' in '%s' must be a single number in",
                        " (0, 1), or NULL."), what), call. = FALSE)
  }
  as.numeric(v)
}


#' @title How Far Down Its Path a Term Reaches
#'
#' @description
#' Reports the fraction of the emptying value to which the path descends. The
#' path runs from the kink that leaves every coefficient of the block at zero
#' down to that fraction of it, so a smaller number reaches a denser fit and a
#' larger one stops sooner.
#'
#' @details
#' It belongs to the term for the same reason that the grid size does: how far
#' the useful range of a hyperparameter extends is a property of the block, and
#' a criterion is applied to every term of the model at once.
#'
#' It is one number per penalty, because only the path over the size of
#' the kink uses it. A bounded hyperparameter is swept over its own interval,
#' and a shape that does not move the kink over a geometric grid above its
#' lower bound; a fraction of an emptying value has no meaning on either axis.
#'
#' [lasso()], [enet()], [scad()] and [mcp()] write their `min_ratio` default
#' into the term, so each reports it whether or not the caller set one, as
#' long as the path over `lambda` is built. A penalty whose `lambda` is held or
#' written out reports nothing, and the constructors reject an explicit
#' `min_ratio` for it.
#' [ridge()], [s()], [te()] and [random()] have no path and report nothing.
#'
#' The keys are [term_penalties()]'s entry names, `""` for a penalty over the
#' whole term.
#'
#' @param term A term, built or not.
#' @param ... Passed to methods. The shipped methods do not use it.
#'
#' @return A named list, one entry per penalty of the term, each a single
#'   number in \eqn{(0, 1)}. The list is empty where every penalty of the term
#'   lacks a swept hyperparameter or a recorded depth.
#'
#' @seealso [term_grid()] for how many values the path visits,
#'   [term_search()] for how several hyperparameters are combined.
#'
#' @examples
#' # The first call shows the default, and the second a deeper path.
#' term_path_min(lasso(~ x))
#' term_path_min(lasso(~ x, min_ratio = 1e-6))
#'
#' # A term with no path, or with no swept hyperparameter, has nothing to
#' # report.
#' term_path_min(ridge(~ x))
#' term_path_min(lasso(~ x, lambda = 2))
#'
#' @export
#' @aliases term_path_min.model_term
term_path_min <- S7::new_generic("term_path_min", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_path_min, model_term) <- function(term, ...) {
  out <- list()
  for (e in .hyper_entries(term)) {
    m <- e$min_ratio
    # the depth belongs to the path over lambda, the size of the kink, and
    # to the built path of an entry whose hyperparameters have other names
    b <- .built(e)
    uses <- if ("lambda" %in% .swept(e) || "lambda" %in% names(e$fixed))
      "lambda" %in% b else length(b) > 0L
    if (length(m) && uses) out[[e$name]] <- as.numeric(m)
  }
  out
}


#' Check How a Term Covers Its Own Hyperparameters
#'
#' @description
#' Validates the sweep that a constructor was given for the term's own kinked
#' hyperparameters.
#'
#' @details
#' The sweep is one word per term, because it states how the hyperparameters
#' are combined with each other, which is not a property of any single
#' hyperparameter. It belongs to the term because a penalty with a kink is
#' fitted by a scheme of its own, and the way that scheme sweeps its
#' hyperparameters is part of it. A criterion is applied to every
#' hyperparameter of the model, including the smooth ones, and most of them
#' cannot use such an argument.
#'
#' @param v What the constructor was given, or `NULL`.
#' @param what The term's label, for the message.
#'
#' @return A character vector of length one, or of length zero.
#'
#' @seealso [term_search()], [check_min_ratio()]
#'
#' @keywords internal
check_search <- function(v, what = "this term") {
  if (is.null(v)) return(character(0))
  if (!is.character(v) || length(v) != 1L || is.na(v) ||
      !v %in% c("grid", "cyclic")) {
    stop(sprintf(paste0("'search' in '%s' must be \"grid\", \"cyclic\" or",
                        " NULL."), what), call. = FALSE)
  }
  v
}


#' @title How a Term Covers Its Own Kinked Hyperparameters
#'
#' @description
#' Reports whether the term's own hyperparameters are swept as a product,
#' `"grid"`, or one at a time while the others are held, `"cyclic"`. A term
#' that names neither word is covered in the way that the fitting layer covers
#' a term by default.
#'
#' @details
#' It matters only for a term with more than one swept hyperparameter:
#' [enet()] by default, and [scad()] and [mcp()] when their shape is left
#' `NULL` to be estimated, since the shape is held by default. There is
#' nothing to combine otherwise, and [lasso()] does not take the argument.
#' Under `"grid"` the cost is the product of the grids of the swept
#' hyperparameters; under `"cyclic"` it is their sum per pass.
#'
#' It is one word per penalty, and it states how the hyperparameters are
#' combined with each other, which is not a property of any single one of
#' them.
#'
#' Between two terms the sweep alternates, whichever search each term names, so
#' `y ~ lasso(X) + enet(R)` costs the sum of the costs of the two blocks. A
#' term with `search = "grid"` therefore does not raise the cost of the other
#' term.
#'
#' The keys are [term_penalties()]'s entry names, `""` for a penalty over the
#' whole term.
#'
#' @param term A term, built or not.
#' @param ... Passed to methods. The shipped methods do not use it.
#'
#' @return A named list, one entry per penalty of the term, each the single
#'   string `"grid"` or `"cyclic"`. The list is empty where every penalty of
#'   the term has fewer than two swept hyperparameters.
#'
#' @seealso [term_grid()] for the sizes that are combined, [term_path_min()]
#'   for the depth of the kink-size axis.
#'
#' @examples
#' # The default is a product of the two axes; cyclic sweeps one at a time.
#' term_search(enet(~ x))
#' term_search(enet(~ x, search = "cyclic"))
#'
#' # SCAD sweeps its shape only when it is estimated.
#' term_search(scad(~ x))
#' term_search(scad(~ x, a = NULL))
#'
#' # With one hyperparameter there is nothing to combine and no argument to
#' # take.
#' term_search(lasso(~ x))
#' try(lasso(~ x, search = "cyclic"))
#'
#' @export
#' @aliases term_search.model_term
term_search <- S7::new_generic("term_search", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_search, model_term) <- function(term, ...) {
  out <- list()
  for (e in .hyper_entries(term)) {
    s <- e$search
    if (length(s) && length(.swept(e)) > 1L) out[[e$name]] <- as.character(s)
  }
  out
}


#' The Hyperparameter Names That a Penalty Factory Produces
#'
#' @description
#' Builds the penalty at a probe size and reads its `params`, returning `NULL`
#' where it cannot be built.
#'
#' @details
#' A factory is a function of the coefficient count, and the names that it
#' produces do not depend on that count, so one coefficient is enough to read
#' them. The probe is wrapped because a factory may reject a size of one, and
#' the two outcomes are kept apart: `NULL` marks a penalty whose names cannot
#' be read here, and an empty vector marks a penalty with no hyperparameter.
#' Every shipped factory returns a penalty, so `NULL` arises only from a
#' factory written by a caller.
#'
#' @param factory A function of the coefficient count returning a penalty, or
#'   a penalty.
#'
#' @return A character vector, empty where the penalty carries no
#'   hyperparameter, or `NULL` where it could not be built.
#'
#' @keywords internal
.factory_params <- function(factory) {
  pen <- if (is.function(factory)) {
    tryCatch(factory(1L), error = function(e) NULL)
  } else {
    factory
  }
  if (is.null(pen)) NULL else pen@params
}


#' Check the Sharing Labels That a Constructor Was Given
#'
#' @description
#' Validates `id` and returns it as a named character vector, one entry per
#' hyperparameter to share, empty where the constructor was given none.
#'
#' @details
#' Two forms are accepted and both have the same meaning. An unnamed single
#' string shares the hyperparameter, and is valid only where the penalty
#' carries one: it is the short form for the frequent case, a smooth or a
#' ridge. A named vector, `c(alpha = "A")`, names the hyperparameter to share,
#' and a penalty carrying several requires this form, since sharing the
#' mixing of an elastic net while leaving its strength free is a different
#' model from the reverse, and the code does not infer which one was meant.
#'
#' The short form on a penalty carrying several hyperparameters is an error
#' that lists them. The names are checked against the penalty's own, so a
#' misspelling is reported where it was written.
#'
#' Sharing does not merge the penalties. They remain two objects estimated at
#' one value, and the code does not check whether they are comparable: the two
#' may weigh quite different matrices, and whether that is the intended model
#' is left to the caller. A shared value is most meaningful between penalties
#' of the same shape over covariates of comparable scale.
#'
#' @param x The `id` argument as given: `NULL`, a single string, or a named
#'   character vector.
#' @param names The hyperparameter names that the term's penalty carries,
#'   against which the names of the labels are checked, or `NULL` where the
#'   penalty does not exist yet. In that case the shape is checked and the
#'   labels are carried as written, to be resolved by a second call where the
#'   names are known.
#' @param what How to name the term in a message.
#'
#' @return A character vector, possibly empty, named by the hyperparameters.
#'   In the deferred case (`names = NULL`) a single unnamed string is returned
#'   unnamed, to be named by the second call.
#'
#' @seealso [term_ids()], which reports it; [check_hyper()] for the argument
#'   that holds a hyperparameter instead of sharing it.
#'
#' @keywords internal
check_ids <- function(x, names, what = "this term") {
  if (is.null(x) || !length(x)) return(character(0))
  if (!is.character(x) || any(is.na(x)) || !all(nzchar(x))) {
    stop(sprintf(paste0("'id' in '%s' must be a non-empty string, or a named\n",
                        "  vector of them with one entry per hyperparameter",
                        " to share."), what), call. = FALSE)
  }
  avail <- names
  nm <- base::names(x)
  # `NULL` names are the deferred case: the term's penalty is not built yet,
  # as in a random effect, whose hyperparameters are the effects'
  # distribution's and are known only once the group is read. The shape is
  # checked here and the labels are carried as written, for the same call to
  # resolve them where the names exist.
  if (is.null(avail)) {
    if (is.null(nm)) return(as.character(x))
    return(stats::setNames(as.character(x), nm))
  }
  # A term whose penalty carries no hyperparameter has nothing to share, and
  # says so rather than storing a label nothing will ever read.
  if (!length(avail)) {
    stop(sprintf(paste0("'%s' has no hyperparameter to share, so 'id' means",
                        " nothing here."), what), call. = FALSE)
  }
  if (is.null(nm) || !all(nzchar(nm))) {
    if (length(x) != 1L || length(avail) != 1L) {
      stop(sprintf(paste0(
        "'id' in '%s' has no names, and an unnamed 'id' shares the\n",
        "  hyperparameter, which is only meaningful where the penalty carries",
        " one. It carries %d: %s.\n",
        "  Name which to share, as id = c(%s = \"...\")."),
        what, length(avail), paste(avail, collapse = ", "), avail[1L]),
        call. = FALSE)
    }
    return(stats::setNames(as.character(x), avail))
  }
  bad <- setdiff(nm, avail)
  if (length(bad)) {
    stop(sprintf(paste0("'id$%s' in '%s' is not a hyperparameter of its",
                        " penalty.\n  It carries: %s."),
                 bad[1L], what, paste(avail, collapse = ", ")), call. = FALSE)
  }
  stats::setNames(as.character(x), nm)
}
