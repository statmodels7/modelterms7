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
#' list, dropping the `NULL` entries, which are the ones to be estimated.
#'
#' @details
#' The check happens at construction, where the caller can see it, instead of
#' at the fit three layers away. A name the penalty does not carry is an error
#' naming the ones it does. That is what turns `mcp(x, a = 3)`, SCAD's shape
#' written on an MCP whose own shape is `gamma`, into a message instead of an
#' argument that lands in `...` and does nothing.
#'
#' The bounds are open, as they are everywhere in the toolkit: a ridge at
#' \eqn{\lambda = 0} is no penalty at all, and an elastic net at
#' \eqn{\alpha = 0} has no kink and is a penalty of another kind.
#'
#' One argument carries three states, each read for one hyperparameter at a
#' time: `NULL` has the path build the grid, one number holds the
#' hyperparameter, and several are the grid itself. This returns the second of
#' the three; [check_values()] returns the third.
#'
#' @param values A named list of the constructor's arguments, `NULL`
#'   where the hyperparameter is to be estimated.
#' @param penalty A \pkg{penalties7} penalty, or a function returning one,
#'   used only to read the names and the bounds.
#' @param what The constructor's name, for the message.
#'
#' @return A named list, keyed by the penalty's own hyperparameter names, of
#'   the values given exactly one number. Empty where every hyperparameter
#'   was left `NULL` or written out as a grid.
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
#' Validates the hyperparameter values a constructor was given as a vector,
#' which a path visits as they stand.
#'
#' @details
#' Several values are neither a held hyperparameter nor a request to build a
#' grid: they are the grid. Each is checked against the penalty's bounds as a
#' held value would be. Nothing else is applied to them: the value that
#' empties the block does not cap them and the depth of the path does not
#' reach them, the caller having said which values to visit.
#'
#' They are sorted, because a path is walked from the emptiest fit towards the
#' fullest and its warm starts follow that order, and duplicates are dropped.
#' Which end of the order is the sparse one depends on the penalty, so the
#' direction is settled where the path is built.
#'
#' @inheritParams check_hyper
#'
#' @return A named list of numeric vectors, one entry per hyperparameter the
#'   caller wrote out.
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
#' The validation is one body because the two states differ in length alone:
#' the name must be one the penalty carries and every value must lie strictly
#' inside its bounds, whether there is one of them or twenty.
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


#' @title The Hyperparameters a Term Holds
#'
#' @description
#' Reports the hyperparameter values the caller fixed in the constructor, one
#' entry per penalty the term carries. A hyperparameter left `NULL` is
#' **estimated** by whatever criterion the fit runs and does not appear here,
#' so an empty result means every one of them is free.
#'
#' @details
#' # The three states of one argument
#'
#' A constructor's hyperparameter argument carries three meanings, read per
#' hyperparameter at a time. `NULL`, the default, has the path build a
#' grid or a marginal criterion estimate the value at the mode. **One number**
#' holds it, and that is what this reports. **Several numbers** are a grid the
#' path visits as they stand, which [term_values()] reports; the
#' hyperparameter is still estimated there, the caller having said where to
#' look, leaving the answer to the fit.
#'
#' # How the entries are keyed
#'
#' The names are the ones [term_penalties()] gives its entries: `""` for a
#' penalty covering the whole term, and `parameter::label` for one a
#' subformula brought in. A term whose penalties come from sub-terms, which is
#' what a structural term with subformulas is, propagates their held values
#' without a method of its own, every entry of that enumeration carrying its
#' own in the field `fixed`.
#'
#' A built term answers from those entries and an unbuilt one from its own
#' `hyper` property, so the answer is the same before and after a build.
#'
#' @param term A term, built or not.
#' @param ... Passed to methods. No shipped method reads anything here.
#'
#' @return A named list, one entry per penalty of the term and named as
#'   [term_penalties()] names it, each entry a named list of single numbers
#'   keyed by the penalty's own hyperparameter names. Empty where the term
#'   holds nothing.
#'
#' @seealso [term_values()] for a written-out grid, [term_grid()] for how many
#'   values a built grid visits, [term_penalties()] for the entries this is
#'   keyed by.
#'
#' @examples
#' # Held, and free.
#' term_hyper(lasso(~ x, lambda = 3))
#' term_hyper(lasso(~ x))
#'
#' # Several values are a grid, not a held value, so they go elsewhere.
#' term_hyper(lasso(~ x, lambda = c(0.1, 1, 10)))
#' term_values(lasso(~ x, lambda = c(0.1, 1, 10)))
#'
#' # One of two held, the other estimated.
#' term_hyper(enet(~ x, alpha = 0.5))
#'
#' # A held value on a sub-term is reported under that entry's name.
#' set.seed(4)
#' d <- data.frame(x = runif(40, 0, 5), g = factor(rep(1:4, 10)))
#' nb <- term_build(nl(~ a * exp(-r * x), a ~ 0 + lasso(~ g, lambda = 2),
#'                     r ~ 1), d)
#' term_hyper(nb)
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
#' Carries the `lambda` a smooth's constructor was given onto the names
#' its penalty uses, and checks it against the positivity every smoothing
#' parameter obeys.
#'
#' @details
#' A one-dimensional smooth and an isotropic tensor product carry one
#' smoothing parameter, named `lambda`; an anisotropic tensor product
#' carries one per margin, `lambda1`, `lambda2` and so on, which is
#' what [penalties7::additive_penalty()] names them. An unnamed
#' vector is read in that order and must be as long as there are margins; a
#' named one may hold some and leave the others to be estimated.
#'
#' @param lambda What the constructor was given, or `NULL`.
#' @param names The penalty's own hyperparameter names.
#' @param what The term's label, for the message.
#'
#' @return A named list of held values.
#'
#' @seealso [check_hyper()], [s()], [te()]
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
#' [random()] builds one of three penalties depending on what it was given: a
#' ridge, a structured prior over a \pkg{parameters7} matrix, or a
#' \pkg{distributions7} family used coordinatewise. Which names there are is
#' therefore not known until the term is built. The shape is checked here and the
#' names against the penalty by [check_hyper()] at that point.
#'
#' A vector entry is a written-out grid, exactly as it is in the constructors
#' that name their hyperparameters, so the length is not
#' checked here either; [.hyper_parts()] splits the two once the penalty
#' exists.
#'
#' @param x A named vector, a named list, or `NULL`.
#' @param what The term's label, for the message.
#'
#' @return A named list, one entry per hyperparameter named, each a numeric
#'   vector of length one or more. Empty for `NULL` or a zero-length input.
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
#' Validates the number of values the constructor was given per
#' hyperparameter, and returns them as a named list with the `NULL`
#' entries dropped.
#'
#' @details
#' How finely a hyperparameter is swept belongs to the term for the same
#' reason as whether it is swept at all: a penalized block of four columns
#' and one of four hundred want different grids, and a criterion applies to
#' every term of the model at once and cannot know which it is looking at.
#' Where a term says nothing the criterion's own default is used.
#'
#' Two values is the smallest grid that is a grid. There is no upper limit
#' beyond the caller's patience: each point of a path is a whole fit, so the
#' cost is linear in the number asked for.
#'
#' @param values A named list of the constructor's arguments, `NULL`
#'   where the criterion's default is wanted.
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
      stop(sprintf(paste0("'%s' has no hyperparameter '%s' to put a grid",
                          " on. It carries: %s."),
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


#' @title The Grid a Term Asks For
#'
#' @description
#' Reports how many values a path visits for each of the term's
#' hyperparameters, one entry per penalty. A hyperparameter the term names
#' nothing for is swept at the fitting layer's own default.
#'
#' @details
#' How finely a hyperparameter is swept belongs to the term because the term is
#' where the penalty is named: a block of four columns and one of four hundred
#' want different grids, and an outer criterion, which is put to every
#' hyperparameter of the model, does not know which it is looking at.
#'
#' Only a hyperparameter with a **path** has a grid. [ridge()], [s()], [te()]
#' and [random()] report nothing, their hyperparameters being estimated at the
#' mode by a marginal criterion; [lasso()] reports its `n_lambda`, and
#' [enet()], [scad()] and [mcp()] report both of theirs. Those constructors
#' write their defaults into the term, so `lasso(~ x)` reports `lambda = 25`
#' where an unset argument would give an empty list.
#'
#' The keys are [term_penalties()]'s entry names, `""` for a penalty over the
#' whole term.
#'
#' @param term A term, built or not.
#' @param ... Passed to methods. No shipped method reads anything here.
#'
#' @return A named list, one entry per penalty of the term, each a named list
#'   of single whole numbers keyed by hyperparameter. Empty where the term
#'   names no grid at all.
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
#' # Two axes, and their two sizes: the kink axis is swept far more finely.
#' term_grid(enet(~ x))
#'
#' # A penalty with no kink has no path, so nothing to count.
#' term_grid(ridge(~ x))
#' term_grid(s(x, k = 5))
#'
#' @export
#' @aliases term_grid.model_term
term_grid <- S7::new_generic("term_grid", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_grid, model_term) <- function(term, ...) {
  ent <- tryCatch(term_penalties(term), error = function(e) list())
  if (length(ent)) {
    out <- list()
    for (e in ent) {
      g <- e$n_values
      if (is.null(g) || !length(g)) next
      out[[if (is.null(e$name)) "" else e$name]] <- g
    }
    if (length(out)) return(out)
  }
  g <- term@grid
  if (!length(g)) return(list())
  stats::setNames(list(g), "")
}


#' Reject a Written-Out Grid Where There Is No Path
#'
#' @description
#' Signals an error when the caller gave several values for a hyperparameter
#' whose penalty has no kink.
#'
#' @details
#' Several values are a grid for a path to visit, and only a penalty with a
#' kink is swept along one: everything else has its hyperparameter read at the
#' mode by a marginal criterion, which would take the vector and do nothing
#' with it. The question is put to the penalty at a probe value of its own
#' hyperparameters, the kink set being structural, so a ridge and a random
#' effect under a Gaussian prior are covered by the same line as a family
#' added later.
#'
#' @param values The written-out grids, as [check_values()] returns
#'   them.
#' @param pen The penalty the term will build, or `NULL`.
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
                      " path to visit them\n  on: its penalty has no kink, so",
                      " the hyperparameter is estimated by the\n  criterion at",
                      " the mode. Give one number to hold it, or NULL to",
                      " estimate it."),
               names(values)[[1L]], what), call. = FALSE)
}


#' @title The Hyperparameter Values a Term Writes Out
#'
#' @description
#' Reports the values a path visits, for each hyperparameter the caller wrote
#' out as a vector instead of leaving to be built. It is the third state of a
#' constructor's hyperparameter argument, beside `NULL` for a built grid and
#' one number for a held value.
#'
#' @details
#' A written-out grid is used as it stands. The value that empties the block
#' does not cap it and [term_path_min()] does not extend it: those two
#' construct a grid, and here there is nothing to construct.
#'
#' A hyperparameter written out is still **estimated**. What the caller fixed
#' is where to look, not the answer, so a fit reports it as chosen by the
#' criterion, and [term_hyper()] does not carry it.
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
#' @param ... Passed to methods. No shipped method reads anything here.
#'
#' @return A named list, one entry per penalty of the term, each a named list
#'   of numeric vectors keyed by hyperparameter, sorted and without
#'   duplicates. Empty where the term wrote nothing out.
#'
#' @seealso [term_hyper()] for a held value, [term_grid()] for a grid to be
#'   built, [term_penalties()] for the entries this is keyed by.
#'
#' @examples
#' # Written out, and reported as the grid to visit.
#' term_values(lasso(~ x, lambda = c(0.1, 1, 10)))
#'
#' # Sorted and deduplicated at construction.
#' term_values(lasso(~ x, lambda = c(10, 0.1, 1, 1)))
#'
#' # One number is a held value, and belongs to the other reporter.
#' term_values(lasso(~ x, lambda = 3))
#' term_hyper(lasso(~ x, lambda = 3))
#'
#' # A penalty with no kink has no path, so several values are refused.
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
#' Validates the fraction of the emptying value a path descends to.
#'
#' @details
#' One number per term, one per hyperparameter being unnecessary: only the
#' path over the size of the kink uses it: a bounded hyperparameter is swept
#' over its own interval and a shape that does not move the kink over a
#' geometric grid above its lower bound, and a fraction of an emptying value
#' means nothing in either.
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
#' Reports the fraction of the emptying value the path descends to. The path
#' runs from the kink that leaves every coefficient of the block at zero down
#' to that fraction of it, so a smaller number reaches a denser fit and a
#' larger one stops sooner.
#'
#' @details
#' It belongs to the term for the same reason the grid size does: how far the
#' useful range of a hyperparameter extends is a property of the block, and a
#' criterion applies to every term of the model at once.
#'
#' It is **one number per penalty**, one per hyperparameter being unnecessary:
#' only the path over the size of the kink uses it. A bounded hyperparameter
#' is swept over its own interval, and a shape that does not move the kink over
#' a geometric grid above its lower bound; a fraction of an emptying value
#' means nothing in either.
#'
#' [lasso()], [enet()], [scad()] and [mcp()] write their `min_ratio` default of
#' `1e-4` into the term, so they report it whether or not the caller set one.
#' [ridge()], [s()], [te()] and [random()] have no path and report nothing.
#'
#' The keys are [term_penalties()]'s entry names, `""` for a penalty over the
#' whole term.
#'
#' @param term A term, built or not.
#' @param ... Passed to methods. No shipped method reads anything here.
#'
#' @return A named list, one entry per penalty of the term, each a single
#'   number in \eqn{(0, 1)}. Empty where the term names none.
#'
#' @seealso [term_grid()] for how many values the path visits,
#'   [term_search()] for how several hyperparameters are combined.
#'
#' @examples
#' # The default, and a deeper path.
#' term_path_min(lasso(~ x))
#' term_path_min(lasso(~ x, min_ratio = 1e-6))
#'
#' # No path, nothing to report.
#' term_path_min(ridge(~ x))
#'
#' @export
#' @aliases term_path_min.model_term
term_path_min <- S7::new_generic("term_path_min", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_path_min, model_term) <- function(term, ...) {
  ent <- tryCatch(term_penalties(term), error = function(e) list())
  if (length(ent)) {
    out <- list()
    for (e in ent) {
      m <- e$min_ratio
      if (is.null(m) || !length(m)) next
      out[[if (is.null(e$name)) "" else e$name]] <- as.numeric(m)
    }
    if (length(out)) return(out)
  }
  m <- term@min_ratio
  if (!length(m)) return(list())
  stats::setNames(list(as.numeric(m)), "")
}


#' Check How a Term Covers Its Own Hyperparameters
#'
#' @description
#' Validates the sweep a constructor was given for the term's own kinked
#' hyperparameters.
#'
#' @details
#' One word per term, one per hyperparameter being meaningless: it says how the
#' hyperparameters are combined with each other, which is not a property any
#' one of them has. It belongs to the term because a penalty with a kink is
#' fitted by a scheme of its own, and how that scheme sweeps its own
#' hyperparameters is part of the scheme. A criterion applies to every
#' hyperparameter of the model, the smooth ones included, and would be
#' carrying an argument most of them cannot read.
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
#' `"grid"`, or one at a time with the others held, `"cyclic"`. A term that
#' names neither is covered the way the fitting layer covers one by default.
#'
#' @details
#' It matters only for a term carrying **more than one** hyperparameter with a
#' kink, so [enet()], [scad()] and [mcp()]; there is nothing to combine
#' otherwise, and [lasso()] does not take the argument. Under `"grid"` the cost
#' is the product of the term's own grids, 25 by 5 at the defaults; under
#' `"cyclic"` it is their sum per pass.
#'
#' It is one word per penalty. It says how the hyperparameters are combined
#' **with each other**, which is not a property any one of them has.
#'
#' Between two terms the sweep alternates whatever each one names, so
#' `y ~ lasso(X) + enet(R)` costs the two blocks added, so one term asking for
#' a product does not make the other pay for it.
#'
#' The keys are [term_penalties()]'s entry names, `""` for a penalty over the
#' whole term.
#'
#' @param term A term, built or not.
#' @param ... Passed to methods. No shipped method reads anything here.
#'
#' @return A named list, one entry per penalty of the term, each the single
#'   string `"grid"` or `"cyclic"`. Empty where the term names none.
#'
#' @seealso [term_grid()] for the sizes being combined, [term_path_min()] for
#'   the depth of the kink-size axis.
#'
#' @examples
#' # The default is a product of the two axes; cyclic sweeps one at a time.
#' term_search(enet(~ x))
#' term_search(enet(~ x, search = "cyclic"))
#'
#' # One hyperparameter, so nothing to combine and no argument to take.
#' term_search(lasso(~ x))
#' try(lasso(~ x, search = "cyclic"))
#'
#' @export
#' @aliases term_search.model_term
term_search <- S7::new_generic("term_search", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_search, model_term) <- function(term, ...) {
  ent <- tryCatch(term_penalties(term), error = function(e) list())
  if (length(ent)) {
    out <- list()
    for (e in ent) {
      s <- e$search
      if (is.null(s) || !length(s)) next
      out[[if (is.null(e$name)) "" else e$name]] <- as.character(s)
    }
    if (length(out)) return(out)
  }
  s <- term@search
  if (!length(s)) return(list())
  stats::setNames(list(as.character(s)), "")
}
