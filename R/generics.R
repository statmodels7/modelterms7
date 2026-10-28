#' @include term_classes.R
NULL

#' @title Build a Term on Data
#'
#' @description
#' Turns a term specification into a built term: an additive term computes its
#' design block from `data`, assigns the coefficient names and records the
#' `build_info` that will reproduce the mapping on other rows; a structural term
#' records whatever its recursion needs, its grouping and its ordering. The
#' returned object is a copy of the specification in which those properties are
#' filled, and the specification is unchanged.
#'
#' @details
#' # What building produces
#'
#' An additive term contributes to the linear predictor through a block, and
#' through a penalty on that block's coefficients when it is penalized:
#'
#' \deqn{\eta = \sum_{t} X_t \beta_t,
#'   \qquad \text{penalized objective} \quad
#'   -\ell(\beta) + \sum_{t} \rho_t(\beta_t; \theta_t).}
#'
#' Building produces \eqn{X_t} from the data and attaches
#' \eqn{\rho_t}. [term_matrix()] then reads the block, [term_penalties()] the
#' penalties, and [term_predict()] reproduces \eqn{X_t} at other rows through
#' the `build_info`.
#'
#' A structural term has no such block. [gas()] records the group to which each
#' row belongs and its place in that group's series; [regime()] and the marginal
#' break-point terms record the quantities that their forward recursion reads.
#' They report themselves through [term_filter()] or [term_loglik()], and
#' [term_matrix()] has no method for them.
#'
#' # Building twice, and building on other data
#'
#' Building is not idempotent in general: a term built again on new data
#' re-derives its factor levels and its knots from those rows.
#' [term_predict()] exists to avoid this re-derivation, and the subset check of
#' [check_term()] tests for it. A term is therefore built once and applied to
#' other rows with [term_predict()].
#'
#' # The default
#'
#' The default method is registered on [model_term()]. It signals the error
#' `"the term class 'X' does not implement term_build()."`, which names the
#' class of the term. It covers both branches: every shipped structural term
#' registers a method of its own, and a structural class written elsewhere
#' that does not register one receives the same error.
#'
#' @param term An object inheriting from [model_term()], built or not. A built
#'   term is rebuilt.
#' @param data A data frame carrying every variable that the term names. Any
#'   other object is rejected by the generic, before dispatch, with the error
#'   `"'data' must be a data frame."`.
#' @param ... Passed to methods. The shipped methods do not use it.
#'
#' @return A built term of the same class as `term`. For an additive term the
#'   `X`, `coef_names` and `build_info` properties are filled, and `penalty` is
#'   filled where the term is penalized; for a structural term the `build_info`
#'   is filled. [term_is_built()] is `TRUE` either way, because it reads the
#'   one of the two records that the branch fills.
#'
#' @seealso [term_predict()] for the block at other rows, [term_matrix()] and
#'   [term_coef_names()] for what a build filled, [term_refresh()] for a block
#'   that moves with its coefficients, and [check_term()] for validating the
#'   result.
#'
#' @examples
#' d <- data.frame(x = rnorm(20), g = factor(rep(letters[1:4], 5)))
#'
#' # A specification carries no block; building fills it.
#' spec <- linpar(~ x + g)
#' built <- term_build(spec, d)
#' c(spec = term_is_built(spec), built = term_is_built(built))
#' dim(term_matrix(built))
#' term_coef_names(built)
#'
#' # The specification is untouched: building returns a copy.
#' term_is_built(spec)
#'
#' # A structural term records its recursion's bookkeeping and no block.
#' g <- term_build(gas(p = 1, q = 1), data.frame(y = rnorm(30)))
#' term_params(g)
#' try(term_matrix(g))
#'
#' # A class without a method signals an error that names the class.
#' Foo <- S7::new_class("Foo", parent = additive_term)
#' try(term_build(Foo(), d))
#'
#' @export
#' @aliases term_build.model_term term_build.structural_term
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

# No structural default: the `model_term` one above names the class, which is
# the useful message. A structural default saying the branch is unimplemented
# was unreachable for gas(), regime() and the marginal break-point terms, all
# of which register term_build(), and misleading for a class written outside
# the package, which is the only caller that can reach it.

#' @title Whether a Term Has Been Built
#'
#' @description
#' Returns `TRUE` for a term that [term_build()] has filled and `FALSE` for a
#' bare specification. [term_matrix()], [term_npar()], [term_coef_names()],
#' [term_predict()] and [plot()] apply this test before reading a design block,
#' so the result determines which error a caller receives from those
#' functions.
#'
#' @details
#' The two branches record being built in different places, so the predicate
#' reads the record of each branch. An additive term is built when it has
#' coefficient names, `length(term@coef_names) > 0L`; a structural term
#' contributes no design columns and so has none to count, and is built when
#' its `build_info` is filled, `length(term@build_info) > 0L`.
#'
#' A built additive term with no columns would give `FALSE`. The shipped
#' constructors do not produce one: a formula that gives no columns, such as
#' `~ 0`, is rejected when [term_build()] runs.
#'
#' The `build_info` property is read directly. It is declared on
#' [additive_term()] and on the abstract [structural_term()], which the
#' three shipped structural classes inherit, so every term inheriting from
#' [model_term()] has one. A structural class written outside the package
#' inherits an empty `build_info`, which reads as not built until [term_build()]
#' fills it, and the predicate returns a logical for such a term as well.
#'
#' @param term An object inheriting from [model_term()]. Any other object is
#'   rejected with the error `"'term' must inherit from 'model_term'."`.
#'
#' @return A single logical, never `NA`.
#'
#' @seealso [term_build()], which makes it `TRUE`; [term_matrix()] and
#'   [term_npar()], which reject a term for which it is `FALSE`.
#'
#' @examples
#' d <- data.frame(x = 1:4)
#' term_is_built(linpar(~ x))
#' term_is_built(term_build(linpar(~ x), d))
#'
#' # The accessors apply this test, so it predicts their error.
#' try(term_matrix(linpar(~ x)))
#'
#' # A structural term carries no design block, so the predicate reads its
#' # build_info instead.
#' gs <- gas(p = 1, q = 1)
#' gb <- term_build(gs, data.frame(y = rnorm(20)))
#' c(spec = term_is_built(gs), built = term_is_built(gb))
#'
#' @seealso [term_build()], [term_matrix()], [term_coef_names()], [term_npar()]
#' @export
term_is_built <- function(term) {
  if (!S7::S7_inherits(term, model_term)) {
    stop("'term' must inherit from 'model_term'.", call. = FALSE)
  }
  # The two branches record being built in different places: an additive term
  # in its coefficient names, a structural one in its `build_info`, having no
  # design block to name coefficients of.
  #
  # `build_info` is declared on additive_term and on the abstract
  # structural_term, so every model_term has one and it can be read rather
  # than asked for. A class written outside the package inherits an empty
  # list, which reads as not built, rather than raising S7's "Can't find
  # property" from inside a predicate whose guard promises a logical.
  if (S7::S7_inherits(term, structural_term)) {
    return(length(term@build_info) > 0L)
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
#' Returns the \eqn{n \times k} block that a built additive term contributes to
#' the linear predictor, with the term's coefficient names as column names and
#' one row per observation of the data on which it was built. The block is the
#' term's `X` property, returned as it is stored.
#'
#' @details
#' The block is not necessarily a base matrix. [random()] builds a grouping
#' indicator as a `dgCMatrix`, since a row belongs to one group and the
#' fraction of nonzero entries is \eqn{1/m} for \eqn{m} groups, and [linpar()],
#' the penalized terms and a smooth with a factor `by` build a sparse block
#' where their `sparse` argument selects it. Code that reads a block therefore
#' tests `is.matrix(x) && is.numeric(x)` or a two-dimensional S4 object. A test
#' that consisted of `is.matrix()` only would reject a block merely for being
#' stored sparsely, because `is.matrix()` is `FALSE` for every \pkg{Matrix}
#' class.
#'
#' For a term whose block moves with its own coefficients, [nl()] and [seg()],
#' this returns the block at the coefficients last committed by
#' [term_refresh()]. [term_jacobian_block()] reports whether that block is a
#' Jacobian or a frozen working linearization.
#'
#' A structural term contributes no block and registers no method, so
#' `term_matrix()` on one signals S7's method-not-found error.
#'
#' @param term A built additive term (see [term_build()]). A specification
#'   signals the error
#'   `"the term has not been built; call term_build(term, data) first."`.
#' @param ... Passed to methods. The shipped methods do not use it.
#'
#' @return The design block: a numeric matrix, or a two-dimensional
#'   \pkg{Matrix} object where the term built one, with `nrow` equal to the
#'   number of observations on which the term was built and `ncol` equal to
#'   [term_npar()].
#'
#' @seealso [term_predict()] for the same mapping at other rows,
#'   [term_coef_names()] for the column names, [term_npar()] for the count,
#'   and [term_refresh()] for a block that moves.
#'
#' @examples
#' d <- data.frame(x = 1:4, g = factor(c("a", "b", "a", "b")))
#'
#' term_matrix(term_build(linpar(~ x), d))
#'
#' # The column names are the term's coefficient names.
#' X <- term_matrix(term_build(linpar(~ x + g, label = "lin"), d))
#' colnames(X)
#'
#' # A grouping indicator comes back sparse, not as a base matrix.
#' R <- term_matrix(term_build(random(~ 1 | g), d))
#' c(class = class(R)[1], is.matrix = is.matrix(R))
#'
#' # A specification has nothing to return.
#' try(term_matrix(linpar(~ x)))
#'
#' @export
#' @aliases term_matrix.additive_term
term_matrix <- S7::new_generic("term_matrix", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_matrix, additive_term) <- function(term, ...) {
  .assert_built(term)
  term@X
}

#' @title Penalty of a Term
#'
#' @description
#' Returns the \pkg{penalties7} penalty attached to the whole of a term's
#' coefficients, or `NULL` where there is none. The hyperparameters, their
#' bounds and links, and every derivative in the coefficients and in the
#' hyperparameters belong to the penalty object; the term only carries it.
#'
#' @details
#' # A penalty over the whole block
#'
#' A term whose penalty reaches part of its parameters returns `NULL` here and
#' declares that penalty through [term_penalties()], which names the parameters
#' that it covers. A [seg()] term and an [nl()] term carry a penalty only where
#' a subformula on one of their parameters declares one, as in `gamma ~ 0 +
#' lasso(~ 1)` over the changes of slope of a [seg()], and their other
#' parameters stay free; [gas()] declares its penalties through
#' [term_penalties()] in the same way. Reporting one of those penalties here
#' would state that it covers the whole block, so this generic returns a penalty
#' only where that holds.
#'
#' [term_penalties()] is therefore the general function, and a fitting layer
#' enumerates penalties through it. This one is the shortcut for the common
#' case.
#'
#' # A specification carries no penalty
#'
#' The penalty is attached at [term_build()], and its width is the number of
#' columns that the data produce, so `term_penalty(ridge(~ x))` is `NULL` and
#' `term_penalty(term_build(ridge(~ x), d))` is the quadratic penalty. The same
#' holds for [term_penalties()] and, through it, for [term_smooth()].
#'
#' The one method is registered on [additive_term()] and reads the `penalty`
#' property. A structural term has no such property, so `term_penalty()` on one
#' signals S7's method-not-found error.
#'
#' @param term An object inheriting from [additive_term()], built or not.
#' @param ... Passed to methods. The shipped methods do not use it.
#'
#' @return A \pkg{penalties7} penalty object, or `NULL` when the term is
#'   unpenalized, when its penalty covers only part of its parameters, or when
#'   it has not been built.
#'
#' @seealso [term_penalties()] for the general form, [term_smooth()] for
#'   whether the result has a kink, [edf()] for what it costs, and
#'   [penalties7::penalty_value()] for what the returned object computes.
#'
#' @examples
#' d <- data.frame(x = rnorm(20), g = factor(rep(c("a", "b"), 10)))
#'
#' # Unpenalized, and unbuilt: both give NULL.
#' term_penalty(linpar(~ x))
#' term_penalty(ridge(~ x))
#'
#' # Built, it is the penalty object itself.
#' p <- term_penalty(term_build(ridge(~ x), d))
#' c(name = p@penalty_name, params = p@params, n_coef = p@n_coef)
#'
#' # A term whose penalty covers part of its parameters returns NULL here
#' # and declares it through term_penalties() instead: a lasso on the
#' # change of slope of a break-point term.
#' sb <- term_build(seg(x, npsi = 1, gamma ~ 0 + lasso(~ 1)),
#'                  data.frame(x = sort(runif(50, 0, 10))))
#' term_penalty(sb)
#' vapply(term_penalties(sb), function(z) z$name, character(1))
#'
#' @export
#' @aliases term_penalty.additive_term
term_penalty <- S7::new_generic("term_penalty", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_penalty, additive_term) <- function(term, ...) {
  term@penalty
}

#' @title Every Penalty That a Term Carries
#'
#' @description
#' Returns the penalties that a term declares: a list of entries, where each
#' entry names a subset of the term's own parameters, the penalty over them, and
#' the hyperparameters that the caller held for that penalty. A fitting layer
#' enumerates penalties through this list, and [term_smooth()] and [edf()] read
#' it.
#'
#' @details
#' # Two ways in which it generalizes [term_penalty()]
#'
#' A term may carry more than one penalty, over different parameters of its
#' own. A panel model with a population value and a departure per group leaves
#' the population value free and shrinks the departures, which is one penalty
#' over part of the parameters and no penalty over the rest. An `nl()` term that
#' develops two of its parameters by two different penalized sub-terms carries
#' two entries.
#'
#' The parameters need not be coefficients of a design block. The
#' persistence of a score-driven term, the nonlinear parameters of [nl()] and
#' the break-point of [seg()] are parameters of the term and of nothing else,
#' and a penalty needs from them only a vector of numbers and their
#' positions. For a structural term `index` gives positions in
#' [term_params()]; for an additive one, columns of the block.
#'
#' # The base method, and the name of an entry
#'
#' The method on [model_term()] derives its result from [term_penalty()], so a
#' term carrying one penalty over its whole block needs no method of its own.
#' Its single entry is named `""`, meaning the whole term.
#'
#' A name is unique within the term and is not the term's own name. Two
#' `ridge()` terms in one formula are two terms with their own
#' hyperparameters, told apart by their names in the formula;
#' [statmodels7::statmod()] records a hyperparameter under the key `term` for a
#' term with one entry and `term::entry` for a term with several.
#'
#' The list itself is unnamed: the key of an entry is its `name` field,
#' read as `e$name`, and `names(entries)` is `NULL`.
#'
#' # A specification carries no entries
#'
#' The penalty is attached at [term_build()], so `term_penalties()` on an
#' unbuilt penalized term is an empty list. Consequently [term_smooth()], which
#' reads this list, returns `TRUE` for `lasso(~ x)` and `FALSE` once the term is
#' built.
#'
#' @param term A built term. An unbuilt one is accepted, and the result lists
#'   the entries that it already declares, usually none.
#' @param ... Passed to methods. The shipped methods do not use it.
#'
#' @return An unnamed list, possibly empty, one element per declared penalty.
#'   Every element is a list with the same nine fields; an empty value stands
#'   for a field that is not set:
#'   \describe{
#'     \item{`name`}{a character label unique within the term, `""` for a
#'       penalty covering the whole of it.}
#'     \item{`index`}{integer positions among the term's own parameters:
#'       columns of the block for an additive term, positions in
#'       [term_params()] for a structural one.}
#'     \item{`penalty`}{a \pkg{penalties7} penalty over exactly those
#'       parameters, so `penalty@n_coef` equals `length(index)`.}
#'     \item{`fixed`}{the hyperparameters that the caller held, a named list,
#'       empty when all of them are estimated.}
#'     \item{`n_values`, `values`, `min_ratio`, `search`}{what the caller
#'       specified about the path over this penalty's hyperparameters, taken
#'       from the term's properties `grid`, `values`, `min_ratio` and
#'       `search`.}
#'     \item{`ids`}{the sharing labels of this penalty's hyperparameters, a
#'       named character vector, empty where the penalty shares none of its
#'       hyperparameters.}
#'   }
#'   An entry that a term carries up from one of its sub-terms keeps the values
#'   and the search that the sub-term holds.
#'
#' @seealso [term_penalty()] for the single-penalty case, [term_components()]
#'   for how a term's columns divide among its own parameters, [term_hyper()]
#'   for only the held values, and [edf()] for what the entries cost.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = sort(runif(50, 0, 10)), g = factor(rep(c("a", "b"), 25)))
#'
#' # A penalty over the whole block gives one entry, named with the empty
#' # string.
#' e <- term_penalties(term_build(ridge(~ x), d))
#' length(e)
#' str(e[[1]][c("name", "index")])
#' e[[1]]$penalty
#'
#' # An unpenalized term declares nothing.
#' term_penalties(term_build(linpar(~ x), d))
#'
#' # A penalty over part of a term's parameters: nl() with one of its two
#' # parameters developed by a lasso, so the entry covers columns 1 and 2
#' # of a block of three.
#' nb <- term_build(nl(~ a * exp(-r * x), a ~ 0 + lasso(~ g),
#'                     start = list(r = 1.3)), d)
#' ent <- term_penalties(nb)
#' vapply(ent, function(z) z$name, character(1))
#' ent[[1]]$index
#' term_npar(nb)
#'
#' # The list is unnamed: the key is the entry's own field.
#' names(ent)
#'
#' # Unbuilt, there is nothing to report yet.
#' term_penalties(lasso(~ x))
#'
#' @export
#' @aliases term_penalties.model_term
term_penalties <- S7::new_generic("term_penalties", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_penalties, model_term) <- function(term, ...) {
  pen <- tryCatch(term_penalty(term), error = function(e) NULL)
  if (is.null(pen)) return(list())
  list(.penalty_entry("", seq_len(term_npar(term)), pen,
                      fixed = term@hyper, n_values = term@grid,
                      values = term@values, min_ratio = term@min_ratio,
                      search = term@search, ids = term@ids))
}

# One entry of term_penalties(), with every field present: what is not
# given is empty.
.penalty_entry <- function(name, index, penalty, fixed = list(),
                           n_values = list(), values = list(),
                           min_ratio = numeric(0), search = character(0),
                           ids = character(0)) {
  list(name = name, index = index, penalty = penalty, fixed = fixed,
       n_values = n_values, values = values, min_ratio = min_ratio,
       search = search, ids = ids)
}

# An entry of a sub-term carried into the term that holds it, under a new
# name and index, with what the sub-term holds travelling with it.
.entry_relabel <- function(e, name, index, penalty = e$penalty) {
  nz <- function(v, empty) if (is.null(v)) empty else v
  .penalty_entry(name, index, penalty, fixed = nz(e$fixed, list()),
                 n_values = nz(e$n_values, list()),
                 values = nz(e$values, list()),
                 min_ratio = nz(e$min_ratio, numeric(0)),
                 search = nz(e$search, character(0)),
                 ids = nz(e$ids, character(0)))
}

#' Where a Held Hyperparameter Belongs
#'
#' @description
#' Splits a term's `hyper` over the entries that the term declares through
#' [term_penalties()] and checks every name against the penalty that carries
#' it. The same applies to `id`.
#'
#' @details
#' The hyperparameters of a term are not always determined when the term is
#' written. A random effect's are the effects' distribution's, and a smooth's
#' are its penalty's, which may come from a factory called at a coefficient
#' count that the data determine. This function therefore makes the check at
#' the build, which is the first point at which the penalties exist.
#'
#' A name is qualified by the entry to which it belongs where there is more
#' than one: the within-group column for a random effect, the level of the
#' factor for a smooth with one smoothing parameter per level. An unqualified
#' name is an error that lists the names available. A value is never recycled
#' over the entries, because a value applied silently to every entry would be
#' wrong for the entries to which it does not belong.
#'
#' @param entries The entries that a term declares, where each entry is a list
#'   with at least `name` and `penalty`.
#' @param hyper The term's `hyper`, already normalized.
#' @param ids The term's `ids`, which are checked against the same names.
#' @param label The term's label, for the message.
#' @param what The object to which the hyperparameters belong, named in the
#'   message.
#'
#' @return The entries, with `fixed`, `values` and `ids` filled in and
#'   checked.
#'
#' @keywords internal
.entry_hyper <- function(entries, hyper, ids, label,
                         what = "the effects' distribution") {
  qual <- function(en) {
    if (nzchar(en$name)) paste0(en$penalty@params, ".", en$name)
    else en$penalty@params
  }
  avail <- unlist(lapply(entries, qual), use.names = FALSE)
  unknown <- setdiff(names(hyper), avail)
  if (length(unknown)) {
    stop(sprintf(paste0(
      "'hyper$%s' is not a hyperparameter of %s in\n",
      "  '%s'. It carries: %s."),
      unknown[1L], what, label, paste(avail, collapse = ", ")),
      call. = FALSE)
  }
  # The labels are resolved HERE and not at the constructor for the reason
  # the values are: a random effect's hyperparameters are the effects'
  # distribution's, qualified by the entry they belong to, and neither is
  # known before the group is read.
  ids <- check_ids(ids, avail, label)
  lapply(entries, function(en) {
    keep <- hyper[intersect(names(hyper), qual(en))]
    keep_id <- ids[intersect(names(ids), qual(en))]
    if (nzchar(en$name) && length(keep)) {
      names(keep) <- sub(sprintf("\\.\\Q%s\\E$", en$name), "", names(keep))
    }
    if (nzchar(en$name) && length(keep_id)) {
      names(keep_id) <- sub(sprintf("\\.\\Q%s\\E$", en$name), "",
                            names(keep_id))
    }
    # several values are a grid for a PATH to visit, and only a penalty with a
    # kink is swept along one -- the same three checks a penalized constructor
    # runs, in the same order
    vals <- check_values(keep, en$penalty, label)
    reject_pathless_values(vals, en$penalty, label)
    en$fixed <- check_hyper(keep, en$penalty, label)
    en$values <- vals
    en$ids <- keep_id
    en
  })
}


#' @title How a Term's Columns Divide Among Its Own Parameters
#'
#' @description
#' Returns one entry per parameter in which a term is written, giving the
#' columns that belong to that parameter and the sub-terms that develop it. The
#' method returns an empty list for a term whose columns are one block with one
#' meaning.
#'
#' @details
#' # Why the term supplies the division
#'
#' A term may be written in parameters of its own: [nl()] in the parameters of
#' \eqn{f}, [seg()] in a slope, a change and a break-point. Any of them may be
#' developed over covariates, and the block then carries several groups of
#' columns meaning different things. Only the term records the group to which
#' each column belongs. A coefficient name is built for a reader, so a consumer
#' reads the division from this function and does not recover it by parsing
#' names.
#'
#' # A parameter may have several sub-terms, of different kinds
#'
#' `seg(x, psi ~ random(~ 1 | id))` develops the break-point with an
#' unpenalized intercept and a random block, so that component's `subs` has two
#' entries and only the second carries a penalty. A consumer reporting a
#' component therefore reports a sequence of sub-terms.
#'
#' `sub_index` splits `index` among those sub-terms, in the order in which the
#' block binds them, which is [component_sub_index()] applied to their
#' coefficient counts.
#'
#' # Structural terms
#'
#' A structural term contributes no design columns, and `index` gives positions
#' in [term_params()], which is the vector that indexes its state, its readable
#' quantities and its variance matrix. In both branches `index` refers to the
#' term's own coefficients.
#'
#' The base method returns an empty list, which is the answer for [linpar()],
#' [s()], [random()] and the five penalized constructors, whose columns are one
#' block with one meaning.
#'
#' @param term A built term. An unbuilt one returns an empty list.
#' @param ... Passed to methods. The shipped methods do not use it.
#'
#' @return A named list with one element per own parameter, named by that
#'   parameter. Each element is a list with
#'   \describe{
#'     \item{`name`}{the parameter's name, the same as the element's.}
#'     \item{`index`}{its columns in the term's block, or its positions in
#'       [term_params()] for a structural term.}
#'     \item{`subs`}{the sub-terms developing it, an empty list for a parameter
#'       that has no sub-terms.}
#'     \item{`sub_index`}{one integer vector per sub-term, splitting `index`
#'       among them; empty where `subs` is.}
#'   }
#'   The result is an empty list for a term whose columns do not divide.
#'
#' @seealso [component_sub_index()] for the split, [term_penalties()] for the
#'   penalties that those sub-terms bring, [term_coef_names()] for the names of
#'   the columns that are divided.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = seq(0.2, 3, length.out = 20),
#'                  g = factor(rep(c("a", "b"), 10)))
#' dd$y <- 2 * exp(-1.3 * dd$x)
#'
#' # The nl() term has two parameters, and the first is developed over a
#' # factor, which gives two columns for `a` and one for `r`.
#' b <- term_build(nl(~ a * exp(-r * x), a ~ 0 + g, start = list(r = 1.3)), dd)
#' lapply(term_components(b), function(z) z$index)
#' term_coef_names(b)
#'
#' # This lists the sub-terms of the developed parameter and splits its
#' # columns among them.
#' term_components(b)$a$sub_index
#'
#' # A term whose columns are one block returns an empty list.
#' term_components(term_build(linpar(~ x), dd))
#'
#' @export
#' @aliases term_components.model_term
term_components <- S7::new_generic("term_components", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_components, model_term) <- function(term, ...) list()
#' @title The Columns of a Component That Belong to Each of Its Sub-Terms
#'
#' @description
#' Splits a developed parameter's columns among the sub-terms developing it, in
#' the order in which the block binds them. The result fills the `sub_index`
#' field of a [term_components()] entry.
#'
#' @details
#' A developed parameter's block is its sub-terms' blocks bound side by side in
#' the order in which they were given, so the division is their coefficient
#' counts cumulated: with counts \eqn{k_1, \dots, k_m} the \eqn{i}-th sub-term
#' takes `index` at positions \eqn{k_1 + \dots + k_{i-1} + 1} to
#' \eqn{k_1 + \dots + k_i}. The counts come from [term_npar()], so every
#' sub-term must be built.
#'
#' The term computes this division, and a consumer does not, because it rests
#' on how the block was assembled.
#'
#' @param index An integer vector: the component's columns in the term's block,
#'   whose length is the sum of the sub-terms' coefficient counts.
#' @param subs A list of built sub-terms developing the parameter, in the order
#'   in which their blocks were bound.
#'
#' @return A list of integer vectors, one per sub-term, partitioning `index` in
#'   order. The result is an empty list when `subs` is empty.
#'
#' @seealso [term_components()], whose methods are its only callers;
#'   [term_npar()] for the counts by which it divides the columns.
#'
#' @keywords internal
component_sub_index <- function(index, subs) {
  if (!length(subs)) return(list())
  k <- vapply(subs, function(s) as.integer(term_npar(s)), integer(1))
  ends <- cumsum(k)
  lapply(seq_along(subs), function(i) index[(ends[[i]] - k[[i]] + 1L):ends[[i]]])
}



#' @title Number of Parameters of a Built Term
#'
#' @description
#' Returns the number of parameters of its own that a built term carries: the
#' number of columns of the design block for an additive term, and the number
#' of entries of [term_params()] for a structural one, which has no block. It
#' is the length of the vector into which a [term_penalties()] entry indexes,
#' the length against which [edf()] is measured, and the length that a fit
#' reserves for the term.
#'
#' @details
#' The two methods are the two branches. On [additive_term()] it is
#' `ncol(term@X)`, so it equals `length(term_coef_names(term))`, and a
#' specification signals an error. On [structural_term()] it is
#' `length(term_params(term))`, which is the count of the term's own
#' parameters after any subformula has expanded: `gas(p = 1, q = 1)` has three,
#' and `gas(p = 1, q = 1, omega ~ z)` has four, the level's intercept and slope
#' in place of the level.
#'
#' @param term A built term (see [term_build()]). An unbuilt additive term
#'   signals the error
#'   `"the term has not been built; call term_build(term, data) first."`.
#' @param ... Passed to methods. The shipped methods do not use it.
#'
#' @return A single whole number.
#'
#' @seealso [term_coef_names()] and [term_params()] for the names behind the
#'   count, [term_matrix()] for the block, [edf()] for the effective degrees of
#'   freedom of the term.
#'
#' @examples
#' d <- data.frame(x = rnorm(20), z = rnorm(20),
#'                 g = factor(rep(letters[1:4], 5)))
#'
#' # An additive term counts columns.
#' b <- term_build(linpar(~ x + g), d)
#' c(npar = term_npar(b), names = length(term_coef_names(b)),
#'   cols = ncol(term_matrix(b)))
#'
#' # A structural term counts its own parameters, and a subformula
#' # replaces one of them by the coefficients developing it.
#' term_npar(term_build(gas(p = 1, q = 1), d))
#' gz <- term_build(gas(p = 1, q = 1, omega ~ z), d)
#' c(npar = term_npar(gz), params = length(term_params(gz)))
#' term_params(gz)
#'
#' @export
#' @aliases term_npar.additive_term term_npar.structural_term
term_npar <- S7::new_generic("term_npar", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_npar, additive_term) <- function(term, ...) {
  .assert_built(term)
  ncol(term@X)
}

S7::method(term_npar, structural_term) <- function(term, ...) {
  length(term_params(term))
}


#' @title Where a Term's Own Coefficients Begin
#'
#' @description
#' Returns the starting coefficients of a built term, one per column of its
#' block. The base method returns zero everywhere, which is an admissible start
#' for a term whose block is a fixed design, since that block does not depend on
#' the coefficients.
#'
#' @details
#' # Why a term supplies a start of its own
#'
#' A term that recomputes its block from its coefficients ([term_refresh()])
#' needs a start of its own, because zero is degenerate for it. In [jump()] the
#' break-point is read off two coefficients as \eqn{-g_k/\delta_k}, so a vector
#' of zeros puts every break-point at the same clamped position, and with more
#' than one break-point the columns of the block coincide. In [seg()] the
#' Jacobian column is \eqn{-\gamma_k \mathbb{1}(x > \psi_k)} and vanishes
#' identically at \eqn{\gamma_k = 0}.
#'
#' Without a `target`, the break-point terms return the start that
#' [term_build()] computed: unit changes, and the break-points at the interior
#' quantiles of the covariate or at the positions that `psi` names. With a
#' `target` of the right length and no developed coefficient, the slope and the
#' changes are instead the least-squares coefficients of the target on an
#' intercept and the term's columns at the starting positions. [nl()] returns
#' the starting values of its own parameters carried through their links.
#'
#' A coefficient of zero has a meaning that depends on the block that the term
#' builds, so this start belongs to the term, as [term_start()] does for a
#' structural term.
#'
#' # The role of `target`
#'
#' `target` is the response carried onto the scale of the predictor to which the
#' term contributes. A term cannot compute it, because the term holds its
#' formula and its charts but not the distribution, the link or the equation,
#' which the fitting layer holds. [nl()] uses it to estimate its own
#' parameters from the data, over a deterministic grid on each free parameter's
#' chart, and the break-point terms use it for the least-squares coefficients
#' described above. Every other term ignores it. The argument is optional;
#' without it each term returns the start that does not depend on the response.
#'
#' @param term A built term (see [term_build()]). An unbuilt one signals an
#'   error through [term_npar()].
#' @param target Optional numeric vector, one value per observation: the
#'   response on the scale of the predictor, `NULL` by default. It is supplied
#'   only where `params_interpretation` records that the response reads the
#'   parameter directly, so a term in a scale's equation receives `NULL`.
#' @param ... Passed to methods. The shipped methods do not use it.
#'
#' @return A numeric vector of length [term_npar()], in the block's column
#'   order.
#'
#' @seealso [term_start()] for a structural term's own parameters,
#'   [term_refresh()] for the terms for which this function exists,
#'   [seg_start()] for the grid rule behind a break-point start.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(50, 0, 10)))
#'
#' # A fixed design starts at zero.
#' term_coef_start(term_build(linpar(~ x), dd))
#'
#' # A break-point term does not: with several break-points, zero would
#' # place all of them at the same position.
#' jb <- term_build(jump(x, npsi = 2), dd)
#' setNames(term_coef_start(jb), term_coef_names(jb))
#'
#' # seg() starts the change at one and the break-point inside the data.
#' sb <- term_build(seg(x, npsi = 1), dd)
#' setNames(term_coef_start(sb), term_coef_names(sb))
#' range(dd$x)
#'
#' @export
#' @aliases term_coef_start.model_term
term_coef_start <- S7::new_generic("term_coef_start", "term",
  function(term, target = NULL, ...) S7::S7_dispatch())

S7::method(term_coef_start, model_term) <- function(term, target = NULL, ...) {
  numeric(term_npar(term))
}


#' @title Coefficients That a Term Draws for Itself
#'
#' @description
#' Replaces, in a draw of the term's block coefficients, the coefficients whose
#' scale the term determines. Every other coefficient is left as it arrived, so
#' a caller drawing a truth passes what it has drawn and receives a vector that
#' it can use.
#'
#' @details
#' # Why a term draws at all
#'
#' A coefficient of a design column is drawn from a normal of the width that
#' its caller chose, and for almost every term this is sufficient: a slope is
#' measured in the response's units against the covariate's, and the term does
#' not change that. A break-point is a position on the covariate's own axis,
#' so a draw made without reading the covariate has no reason to fall inside
#' the data. A draw from a standard normal falls outside the range of a
#' covariate that is measured on another scale, and the confinement then moves
#' the position to the edge of the interval, where one of the two segments
#' holds a twentieth of the rows and the data do not reveal a break. The
#' break-point terms therefore draw their own positions on the covariate's axis.
#' [term_draw()] applies the same division to the own parameters of a
#' structural term, and each quantity is drawn by the term that defines it.
#'
#' # The width
#'
#' The width follows from the spacing of the positions from which the draw
#' spreads. Writing \eqn{w} for the width of the confinement interval (from the
#' 5th to the 95th percentile of the covariate) and \eqn{K} for the number of
#' break-points, those positions are the \eqn{j/(K+1)} quantiles of the
#' covariate. For a covariate spread evenly over the interval, consecutive
#' positions are \eqn{w/(K+1)} apart. Requiring three standard deviations to
#' stay within half of that spacing, so that two neighbors cross only if both
#' travel further than a draw usually goes, gives
#'
#' \deqn{\sigma = \frac{sd \cdot w}{6(K+1)},}
#'
#' which at one break-point is \eqn{w/12} and, for a covariate spread evenly
#' over the interval, keeps the draw inside \eqn{[0.25w, 0.75w]} of the
#' interval at three standard deviations. For another distribution of the
#' covariate the spacing is a rule of thumb, and the width is still computed
#' from \eqn{w}.
#'
#' # Where the coefficients go
#'
#' A continuous or smoothed construction holds the break-point directly, so
#' the drawn position is written where it stands. A discontinuous one does
#' not hold it: it reads it off \eqn{\psi_k = -g_k/\delta_k}, so the
#' coefficient written is \eqn{-\delta_k \psi_k} at the change of level that
#' the rest of the draw supplies. That is why the whole vector is passed in
#' and not only the term's own slice, and why the position is not the ratio of
#' two independently drawn normals, which has no mean and may fall anywhere.
#' For [jseg()] the position read from the written coefficient is near the
#' drawn one and not equal to it.
#'
#' Where the break-point carries a development, the drawn positions are
#' projected onto that development's own design by least squares, so any
#' subformula is served and the spread that the design induces is the spread
#' that a group actually receives.
#'
#' @param term A built term (see [term_build()]).
#' @param coef The coefficients drawn so far, of length [term_npar()], in the
#'   block's column order.
#' @param sd The width of the draw, `1` by default, multiplying the derived
#'   width above.
#' @param ... Passed to methods. The shipped methods do not use it.
#'
#' @return A list with `coef`, the coefficients in which the term's own have
#'   been replaced; `index`, the positions that were replaced, empty if the
#'   term replaced no coefficient; and `scale`, one width per replaced position.
#'   A caller that reports the truth of a prior over those coordinates reads
#'   `scale`, since the values there are no longer a draw from that prior.
#'
#' @seealso [term_draw()] for a structural term's own parameters,
#'   [term_coef_start()] for the point from which this draw spreads, [seg_psi()]
#'   for the positions that it produces.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(200, 0, 10)))
#'
#' # For an ordinary term the draw is returned unchanged.
#' lb <- term_build(linpar(~ x), dd)
#' identical(term_coef_draw(lb, c(1, 2))$coef, c(1, 2))
#'
#' # A break-point term places its own on the covariate's scale, where a
#' # draw from a standard normal would often fall outside the data.
#' sb <- term_build(seg(x, npsi = 1), dd)
#' d <- term_coef_draw(sb, stats::rnorm(term_npar(sb)))
#' setNames(d$coef, term_coef_names(sb))
#' range(dd$x)
#'
#' @export
#' @aliases term_coef_draw.model_term
term_coef_draw <- S7::new_generic("term_coef_draw", "term",
  function(term, coef, sd = 1, ...) S7::S7_dispatch())

S7::method(term_coef_draw, model_term) <- function(term, coef, sd = 1, ...) {
  list(coef = coef, index = integer(0), scale = numeric(0))
}

#' @title Is a Term's Block the Jacobian of Its Contribution?
#'
#' @description
#' Returns `TRUE` when the block that a term reports is the exact derivative of
#' its contribution in its own coefficients and `FALSE` when it is a working
#' linearization in which quantities are frozen at the previous iterate. The
#' base method returns `TRUE`, which a fixed design satisfies trivially and
#' which [nl()] and [seg()] satisfy by construction. A sharp [jump()] or
#' [jseg()] returns `FALSE` unless it is smoothed or held by [seg_hold()].
#'
#' @details
#' # Use of the result
#'
#' Where the block is a Jacobian, a scoring step on it is a Gauss-Newton step,
#' the fitting layer can run a line search on the model's own objective, and
#' the term is fitted inside the same system as everything else.
#'
#' Where it is a working linearization the fixed-point iteration of Fasola,
#' Muggeo and Kuchenhoff (2018) is not a descent method on the model's
#' objective. Its early steps increase the objective by design, under a scaling
#' factor that anneals, so forcing a sufficient decrease on it stalls the
#' iteration. Such a term is fitted by alternating exact working fits at the
#' frozen block with [term_refresh()], and [term_converged()] reports its
#' convergence instead of a score.
#'
#' The terms that return `FALSE` are the discontinuous constructions [jump()]
#' and [jseg()], when sharp and not held. Their weight
#' \eqn{W = 1/(2\lvert\tilde x - \psi\rvert)} is held at the previous
#' break-point, and the position is read off two coefficients and is not one
#' itself.
#'
#' # Smoothing and holding change the answer
#'
#' `jump(x, smoothed = ...)` replaces the step by a smooth surrogate, and the
#' break-point becomes an ordinary parameter with a true Jacobian, so a
#' smoothed break-point term returns `TRUE`. A term put in its held state by
#' [seg_hold()] also returns `TRUE`, because its positions are fixed and its
#' block is the exact design at them. A fitting layer therefore routes these
#' terms without a special case.
#'
#' @param term A term, built or not: the answer is a property of the
#'   construction.
#' @param ... Passed to methods. The shipped methods do not use it.
#'
#' @return A single logical.
#'
#' @references
#' Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
#' iterative algorithm for change-point detection in abrupt change models.
#' *Computational Statistics*, 33, 997--1015.
#'
#' @seealso [term_refresh()] for the block that is recomputed,
#'   [term_converged()] for the verdict on such a term, [seg()], [jump()] and
#'   [jseg()].
#'
#' @examples
#' # The continuous construction has a Jacobian block; the two sharp
#' # discontinuous ones report a working linearization.
#' vapply(list(seg = seg(x), jump = jump(x), jseg = jseg(x)),
#'        term_jacobian_block, logical(1))
#'
#' # A fixed design is a Jacobian trivially, and so is nl().
#' term_jacobian_block(linpar(~ x))
#' term_jacobian_block(nl(~ a * x, start = list(a = 1)))
#'
#' # Smoothing the step makes the break-point an ordinary parameter.
#' term_jacobian_block(jump(x, smoothed = numericals7::smooth_probit()))
#'
#' @export
#' @aliases term_jacobian_block.model_term
term_jacobian_block <- S7::new_generic("term_jacobian_block", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_jacobian_block, model_term) <- function(term, ...) TRUE

#' @title Where a Term's Objective Has a Kink
#'
#' @description
#' Returns the positions, among the term's own coefficients, of those
#' coefficients in which the term's contribution to the objective is not
#' differentiable at the given coefficients. A fitting layer reads them where a
#' line search rejects every step along a Newton or Gauss-Newton direction:
#' holding those coordinates at their values and solving for the rest is the
#' step at a kink, as holding a coordinate at a bound is the step at a
#' boundary.
#'
#' @details
#' A term whose contribution is a smooth function of its coefficients has no
#' kink, and that is the base method's answer. Two shipped terms report kinks:
#' the continuous break-point construction, [seg()], whose truncated
#' line \eqn{(x-\psi)_+} is not differentiable in \eqn{\psi} at an
#' observation, and [nl()], which reports the kinks of a [seg()] sub-term that
#' develops one of its parameters. A smoothed [seg()] and the discontinuous
#' [jump()] and [jseg()] return `integer(0)`. The kinks of a penalty are not
#' reported here; [term_smooth()] reports them.
#'
#' In a penalized fit of `seg(x, psi ~ random(~ 1 | id))`, the mode can put one
#' group's break-point exactly on an observation, where the one-sided slopes
#' differ. A scoring step on all the coefficients can then be rejected at every
#' length. Holding the coordinates that move that break-point and solving for
#' the rest gives a step that reaches the mode.
#'
#' @param term A built term.
#' @param coef The term's coefficients, `NULL` for the stored ones.
#' @param ... Passed to methods. The shipped methods do not use it.
#'
#' @return An integer vector of positions in the term's coefficients, possibly
#'   empty.
#'
#' @seealso [term_jacobian_block()], [term_refresh()], [seg()].
#'
#' @examples
#' d <- data.frame(x = seq(0, 10, length.out = 41))
#' d$y <- 1 + 0.5 * d$x + 2 * pmax(d$x - 6, 0)
#' term_kinks(term_build(linpar(~ x), d))
#' term_kinks(term_build(seg(x), d), c(0.5, 2, 6))
#'
#' @export
#' @aliases term_kinks.model_term
term_kinks <- S7::new_generic("term_kinks", "term",
  function(term, coef = NULL, ...) S7::S7_dispatch())

S7::method(term_kinks, model_term) <- function(term, coef = NULL, ...) {
  integer(0)
}

#' @title The Linear Constraints a Kink Places on a Term's Coefficients
#'
#' @description
#' Returns the rows of the linear constraints that hold the term's objective on
#' its kinks at the given coefficients: one row per constraint and one column
#' per coefficient of the term. A perturbation \eqn{d} of the coefficients keeps
#' the kinks where they are if and only if \eqn{A d = 0}.
#'
#' @details
#' Where [term_kinks()] names a coordinate, the kink is held by holding that
#' coordinate, and the row is the corresponding unit vector; that is the base
#' method's answer. A break-point whose position is developed over a design,
#' \eqn{\psi_i = z_i^\top\kappa}, sits on an observation \eqn{j} along the
#' hyperplane \eqn{z_j^\top\kappa = x_j}, so the row of that constraint is
#' \eqn{z_j} and the coefficients can still move along the hyperplane. The
#' [seg()] method returns those rows and the [nl()] method carries the rows of
#' a [seg()] sub-term to the term's coefficients. Rows that repeat are given
#' once.
#'
#' @param term A built term.
#' @param coef The term's coefficients, `NULL` for the stored ones.
#' @param ... Passed to methods. The shipped methods do not use it.
#'
#' @return A numeric matrix with one column per coefficient of the term and one
#'   row per constraint, with no rows where the term has no kink at `coef`.
#'
#' @seealso [term_kinks()], [seg()].
#'
#' @examples
#' d <- data.frame(x = seq(0, 10, length.out = 41))
#' term_kink_constraints(term_build(seg(x), d), c(0.5, 2, 6))
#' term_kink_constraints(term_build(seg(x), d), c(0.5, 2, 6.1))
#'
#' @export
#' @aliases term_kink_constraints.model_term
term_kink_constraints <- S7::new_generic("term_kink_constraints", "term",
  function(term, coef = NULL, ...) S7::S7_dispatch())

S7::method(term_kink_constraints, model_term) <- function(term, coef = NULL,
                                                         ...) {
  p <- if (is.null(coef)) term_npar(term) else length(coef)
  k <- term_kinks(term, coef = coef, ...)
  A <- matrix(0, length(k), p)
  if (length(k)) A[cbind(seq_along(k), k)] <- 1
  A
}

#' @title The Covariance Label of a Term
#'
#' @description
#' Returns the label that marks a term whose coefficients share a covariance
#' block with those of other terms, or `NA_character_` where the term carries
#' no label.
#'
#' @details
#' The label is written in the middle position of a bar formula,
#' `random(~ 1 + x | u | id)`, following \pkg{brms}: the last position is the
#' grouping variable and the middle one the label. Terms carrying the same
#' label and the same grouping describe effects that are correlated with each
#' other, and a fitting layer collects them into one block.
#'
#' The label is a name and not data. It is read as a symbol and never
#' evaluated, so a column of the data with the same name is not consulted. A
#' term that carries a label declares no penalties through
#' [term_penalties()]: the prior of the shared covariance block is built by the
#' fitting layer.
#'
#' The base method returns `NA_character_`, which is the answer for every term
#' but [random()], so a term written outside this package needs no method here.
#'
#' @param term A term, built or not.
#' @param ... Passed to methods. The shipped methods do not use it.
#'
#' @return A single string, or `NA_character_` where the term carries no label.
#'
#' @seealso [random()], the one constructor that reads a label;
#'   [term_penalties()] for the penalties that a term declares on its own.
#'
#' @examples
#' # No label: the effects of this term correlate with nothing outside it.
#' term_tag(random(~ 1 | g))
#' term_tag(linpar(~ x))
#'
#' # The label is read from the middle of the two bars.
#' term_tag(random(~ 1 + x | u | g))
#'
#' # It is a name, not data: a column called `u` is not read.
#' d <- data.frame(x = rnorm(6), u = 1:6, g = factor(rep(c("a", "b"), 3)))
#' identical(term_coef_names(term_build(random(~ 1 | u | g), d)),
#'           term_coef_names(term_build(random(~ 1 | g), d)))
#'
#' @export
#' @aliases term_tag.model_term
term_tag <- S7::new_generic("term_tag", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_tag, model_term) <- function(term, ...) NA_character_

#' @title The Grouping by Which a Term's Coefficients Are Indexed
#'
#' @description
#' For a term whose coefficients are one set per level of a grouping variable,
#' returns the expression that defines that variable, its levels, and the
#' number of columns that each level carries. The result is `NULL` for every
#' other term.
#'
#' @details
#' Two labeled terms belong to the same covariance block only if they share a
#' grouping, and deciding that requires both the expression and the levels: the
#' expression by itself is not enough, since `droplevels(id)` and `id` are
#' different expressions for the same grouping and `id` under two subsets is the
#' same expression for two different ones. The levels settle it, and the
#' expression makes a message readable.
#'
#' `dim` is the number of columns that one level carries, which together with
#' the levels gives the layout of the block: the coefficients are group-major,
#' so level \eqn{i} occupies positions \eqn{(i-1)d + 1} to \eqn{id}, where
#' \eqn{d} is `dim`. A consumer stacking two such blocks into one covariance
#' needs this layout, and `names` carries the within-group design's own column
#' names in the same order, so that each coordinate of the covariance that it
#' builds is identified as belonging to `x` or to the intercept and is not
#' merely numbered.
#'
#' The base method returns `NULL`, which is the answer for every term but
#' [random()].
#'
#' @param term A built term. An unbuilt one returns `NULL`.
#' @param ... Passed to methods. The shipped methods do not use it.
#'
#' @return A list with `expr` (a language object), `levels` (a character
#'   vector), `dim` (a single integer) and `names` (`dim` column names), or
#'   `NULL`.
#'
#' @seealso [term_tag()] for the label that, with this, identifies a
#'   covariance block; [random()] for the formula from which both are read.
#'
#' @examples
#' d <- data.frame(x = rnorm(6), g = factor(rep(c("a", "b", "c"), 2)))
#' b <- term_build(random(~ 1 + x | g), d)
#' g <- term_group(b)
#' c(dim = g$dim, levels = length(g$levels))
#' g$names
#'
#' # every other kind of term has none
#' term_group(term_build(linpar(~ x), d))
#'
#' @export
#' @aliases term_group.model_term
term_group <- S7::new_generic("term_group", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_group, model_term) <- function(term, ...) NULL

#' @title The Within-Group Design of a Term With Grouped Effects
#'
#' @description
#' Returns, at new rows, the design row that each group's effect multiplies.
#'
#' @details
#' A term whose coefficients are one effect per group, [random()], builds
#' its block by interacting a within-group design with the group indicators.
#' This returns the first factor by itself, and that factor does not depend on
#' the grouping: a caller who integrates the effect of a group that was absent
#' from the fitting data reads \eqn{z_i^\top b} from it, where \eqn{b} is drawn
#' from the prior.
#'
#' The base method returns `NULL`, which is the answer for every term but
#' [random()].
#'
#' @param term A built term.
#' @param newdata A data frame.
#' @param ... Passed to methods. The shipped methods do not use it.
#'
#' @return A numeric matrix of `nrow(newdata)` rows, or `NULL`.
#'
#' @seealso [term_group()] for the layout of the block,
#'   [term_within.RandomTerm()] for the method.
#'
#' @examples
#' d <- data.frame(x = rnorm(6), g = factor(rep(c("a", "b", "c"), 2)))
#' term_within(term_build(random(~ 1 + x | g), d), d[1:2, ])
#' term_within(term_build(linpar(~ x), d), d)
#'
#' @export
#' @aliases term_within.model_term
term_within <- S7::new_generic("term_within", "term",
  function(term, newdata, ...) S7::S7_dispatch())

S7::method(term_within, model_term) <- function(term, newdata, ...) NULL

#' @title Which of a Term's Hyperparameters Are Shared
#'
#' @description
#' Returns the sharing labels that a term carries. Each element of the answer is
#' a named character vector whose names are the hyperparameters of one penalty
#' entry and whose values are the labels. The terms carrying the same label for
#' the same hyperparameter estimate one value.
#'
#' @details
#' # What sharing is, and what it is not
#'
#' Sharing does not merge the penalties. They stay two objects with the same
#' number in them, so each keeps its own block, its own coefficients and its own
#' effective degrees of freedom, and only the value that a criterion estimates
#' is shared. Two smooths with one label are one smoothing parameter over two
#' curves, which is \pkg{mgcv}'s `id`, and two lassos with one label are one
#' threshold swept along one path.
#'
#' # Comparability of the penalties
#'
#' The comparability of two penalties is not checked, deliberately. A shared
#' value multiplies each penalty as it is, so it has the usual meaning for
#' penalties of the same shape over covariates of comparable scale, and the
#' caller is responsible for the comparison in other cases. The intended case
#' is two smooths of the same basis and dimension.
#'
#' # Where it may be written
#'
#' On any term that carries a penalty, including one written inside a
#' subformula: a `seg()` developing its slope change and its break-point over
#' two ridges may give the two one label and estimate one shrinkage for them,
#' and a `jseg()` developing its slope change and its level change may do the
#' same.
#'
#' The shape is that of [term_hyper()]: one element of the answer for each entry
#' of [term_penalties()] that shares anything, keyed by the name of the entry,
#' so the two functions read alike. A built term derives its result from its
#' entries, which allows a structural term to report what the terms of its
#' subformulas share without a method of its own, and an unbuilt term derives
#' its result from its own property.
#'
#' @param term A term, built or not.
#' @param ... Passed to methods. The shipped methods do not use it.
#'
#' @return A list, one named character vector per penalty entry that shares
#'   anything, keyed by that entry's name. The list is empty where the term
#'   shares nothing.
#'
#' @seealso [term_hyper()] for holding a hyperparameter instead of sharing it,
#'   [term_penalties()] for the entries to which the labels belong.
#'
#' @examples
#' # One hyperparameter, so the label needs no name.
#' term_ids(s(x, id = "L"))
#'
#' # With several hyperparameters the label needs a name: this shares the
#' # mixing and leaves the strength free.
#' term_ids(enet(~ x, id = c(alpha = "A")))
#'
#' # A term sharing nothing returns an empty list.
#' term_ids(s(x))
#' term_ids(linpar(~ x))
#'
#' @export
#' @aliases term_ids.model_term
term_ids <- S7::new_generic("term_ids", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_ids, model_term) <- function(term, ...) {
  ent <- tryCatch(term_penalties(term), error = function(e) list())
  if (length(ent)) {
    out <- list()
    for (e in ent) {
      v <- e$ids
      if (is.null(v) || !length(v)) next
      out[[if (is.null(e$name)) "" else e$name]] <- v
    }
    if (length(out)) return(out)
  }
  v <- term@ids
  if (!length(v)) return(list())
  stats::setNames(list(v), "")
}

#' @title Coefficient Names of a Built Term
#'
#' @description
#' Returns the names of a built additive term's coefficients, one per column of
#' its block and in the block's own order. A fit reports these names, so that a
#' reader of a coefficient table can tell which term produced which row.
#'
#' @details
#' A term with a non-empty `label` prefixes every name with it and a dot, so
#' `linpar(~ x, label = "lin")` gives `lin.(Intercept)` and `lin.x`. The five
#' penalized constructors set the label to the constructor's own name by
#' default, which is why a ridge over `x` and `g` reads `ridge.x`, `ridge.ga`
#' and so on; [linpar()] sets none, so its names are
#' [stats::model.matrix()]'s unchanged.
#'
#' The names are assigned at build time and recorded, so [term_predict()] also
#' labels its columns with them at other rows. Uniqueness is not
#' enforced here; [check_term()] checks it.
#'
#' @param term A built additive term (see [term_build()]). A specification
#'   signals the error
#'   `"the term has not been built; call term_build(term, data) first."`.
#' @param ... Passed to methods. The shipped methods do not use it.
#'
#' @return A character vector of length [term_npar()], in column order.
#'
#' @seealso [term_npar()] for the count, [term_matrix()] for the block that they
#'   name, [term_params()] for a structural term's parameter names instead.
#'
#' @examples
#' d <- data.frame(x = 1:4, g = factor(c("a", "b", "a", "b")))
#'
#' # linpar() takes model.matrix()'s names as they come.
#' term_coef_names(term_build(linpar(~ x + g), d))
#'
#' # A label prefixes every one of them.
#' term_coef_names(term_build(linpar(~ x, label = "lin"), d))
#'
#' # The penalized constructors label themselves by default.
#' term_coef_names(term_build(ridge(~ x + g), d))
#'
#' # They are the column names of the block, and of a prediction.
#' b <- term_build(ridge(~ x + g), d)
#' identical(term_coef_names(b), colnames(term_matrix(b)))
#'
#' @export
#' @aliases term_coef_names.additive_term
term_coef_names <- S7::new_generic("term_coef_names", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_coef_names, additive_term) <- function(term, ...) {
  .assert_built(term)
  term@coef_names
}

#' @title Whether a Term's Penalized Objective Is Smooth
#'
#' @description
#' Returns `TRUE` when the term's contribution to the penalized objective is
#' differentiable in the coefficients. A fitting layer reads the result to split
#' the coefficient vector into the block that a classical optimizer handles and
#' the block that needs a proximal step or a coordinate descent.
#'
#' @details
#' # The answer comes from the penalties
#'
#' The term does not declare the answer. For every entry of [term_penalties()]
#' the function evaluates the kink set through [penalties7::penalty_kinks()], at
#' a probe value inside each hyperparameter's bounds, and the result is `FALSE`
#' as soon as one entry has a point in its kink set. An unpenalized term is
#' therefore smooth, a ridge or a Gaussian prior is smooth, and lasso, SCAD, MCP
#' and the elastic net are not. The answer always agrees with the penalties
#' that the term carries.
#'
#' The probe is any admissible value, because the kink set is structural: the
#' midpoint of a bounded interval, one step inside a half-bounded one, zero
#' where the interval is the whole line.
#'
#' # The answer covers the whole term
#'
#' For a term that carries a penalty over part of its parameters and no penalty
#' over the rest, the result is determined by the part that is penalized.
#' `nl(~ a * exp(-r * x), a ~ 0 + lasso(~ g))` is not smooth: the coefficients
#' developing `a` sit at a kink, although `r` is unpenalized.
#'
#' # A specification is always smooth
#'
#' The penalty is attached at [term_build()], so `term_penalties()` on an
#' unbuilt term is empty and this function returns `TRUE` whatever penalty the
#' term will carry. `term_smooth(lasso(~ x))` is `TRUE` and
#' `term_smooth(term_build(lasso(~ x), d))` is `FALSE`. The function is
#' therefore applied to a built term.
#'
#' @param term An object inheriting from [model_term()]. Build it first, or the
#'   answer is `TRUE` by default.
#' @param ... Passed to methods. The shipped methods do not use it.
#'
#' @return A single logical, never `NA`.
#'
#' @seealso [term_penalties()] for the entries that it evaluates,
#'   [penalties7::penalty_kinks()] for the kink set, [edf()], which uses
#'   the same split to count degrees of freedom.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(30), g = factor(rep(c("a", "b"), 15)))
#'
#' # Unpenalized and quadratically penalized terms are smooth.
#' vapply(list(linpar(~ x), ridge(~ x), s(x, basis7::bspline_smooth(k = 5)), random(~ 1 | g)),
#'        function(t) term_smooth(term_build(t, d)), logical(1))
#'
#' # The four kinked penalties are not.
#' vapply(list(lasso(~ x), scad(~ x), mcp(~ x), enet(~ x)),
#'        function(t) term_smooth(term_build(t, d)), logical(1))
#'
#' # A kink on part of a term's parameters makes the term non-smooth.
#' d$y <- 2 * exp(-1.3 * d$x)
#' nb <- term_build(nl(~ a * exp(-r * x), a ~ 0 + lasso(~ g),
#'                     start = list(r = 1.3)), d)
#' term_smooth(nb)
#' c(npar = term_npar(nb), penalized = length(term_penalties(nb)[[1]]$index))
#'
#' # Unbuilt, there is no penalty to read yet.
#' c(spec = term_smooth(lasso(~ x)), built = term_smooth(term_build(lasso(~ x), d)))
#'
#' @export
#' @aliases term_smooth.model_term
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
  stop(sprintf("'%s' is not a known penalty.", as.character(kind)[1L]),
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
#' Applies a built term's recorded mapping to new rows, returning the block that
#' the term would have produced had those rows been in the data on which it was
#' built. Factor levels, contrasts, spline knots and a Demmler-Reinsch
#' reparametrization come from the `build_info` and are reused; the spreads of a
#' standardization enter the penalty only and are not applied to the block. A
#' factor level that the `build_info` does not contain is rejected, except for the
#' grouping levels of a [random()] term with `unseen = "zero"`.
#'
#' @details
#' # The identity that makes it useful
#'
#' The block returned is \eqn{\tilde{X}_t} such that
#' \eqn{\tilde{\eta} = \tilde{X}_t \beta_t} is the term's contribution at the
#' new rows, at the coefficients that the model already carries. The identity
#' holds because the mapping is reused. A rebuilt encoding gives a block of the
#' same shape that multiplies the same coefficients but means something else: a
#' factor whose new rows omit a level loses a column, and a basis
#' rebuilt on a narrower range is a different set of functions.
#'
#' For a smooth, predicting on a subset of the fitting rows reproduces the
#' corresponding rows of the original block exactly, whereas rebuilding the term
#' on the subset gives a block that multiplies the coefficients by different
#' functions. The subset check of [check_term()] makes this comparison.
#'
#' Predicting on the fitting data returns the block itself, so
#' `term_predict(b, data)` and `term_matrix(b)` agree to the last bit.
#'
#' # Which terms have a method
#'
#' The methods are registered for [linpar()], the penalized terms, [random()],
#' the smooths, [nl()] and the break-point terms. A structural term contributes
#' no block and has no method, so `term_predict()` on one signals S7's
#' method-not-found error; [term_continue()] is the corresponding operation
#' there.
#'
#' @param term A built additive term (see [term_build()]). A specification
#'   signals the error
#'   `"the term has not been built; call term_build(term, data) first."`.
#' @param newdata A data frame carrying every variable that the term names, with
#'   any number of rows. Any other object is rejected by the generic, before
#'   dispatch, with the error `"'newdata' must be a data frame."`. A factor here
#'   need carry only the levels that its own rows use; the rest come from the
#'   `build_info`.
#' @param ... Passed to methods. The method for [random()] reads `unseen`. Its
#'   default `"error"` rejects a grouping level that was absent at the build
#'   and a missing grouping value; with `unseen = "zero"` such a row receives a
#'   row of zeros. The other shipped methods do not use it.
#'
#' @return A block of `nrow(newdata)` rows and [term_npar()] columns, in the
#'   same storage that the term built: a numeric matrix, or a \pkg{Matrix}
#'   object where the block is sparse.
#'
#' @seealso [term_matrix()] for the block on the fitting data, [term_build()]
#'   for what records the `build_info`, [check_term()] for the check that takes
#'   this identity as its subject, and [term_continue()] for a structural term.
#'
#' @examples
#' d <- data.frame(x = 1:6, g = factor(rep(c("a", "b", "c"), 2)))
#' b <- term_build(linpar(~ x + g), d)
#'
#' # New rows go through the same mapping as the fitting rows.
#' term_predict(b, data.frame(x = c(0.5, 2.5), g = factor(c("a", "c"))))
#'
#' # On the fitting data it returns the block itself.
#' all.equal(term_predict(b, d), term_matrix(b))
#'
#' # A subset that drops a level keeps the build_info's columns, where a
#' # rebuild would lose one.
#' nd <- droplevels(d[d$g != "c", ])
#' dim(term_predict(b, nd))
#' dim(model.matrix(~ x + g, nd))
#'
#' # A basis is not replaced on the narrower range: reapplying agrees with
#' # the original rows exactly, rebuilding does not.
#' d2  <- data.frame(x = seq(0, 1, length.out = 60))
#' bs  <- term_build(s(x, basis7::bspline_smooth(k = 6)), d2)
#' X   <- term_matrix(bs)
#' sub <- 1:20
#' max(abs(term_predict(bs, d2[sub, , drop = FALSE]) - X[sub, ]))
#' max(abs(term_matrix(term_build(s(x, basis7::bspline_smooth(k = 6)), d2[sub, , drop = FALSE])) -
#'         X[sub, ]))
#'
#' @export
term_predict <- S7::new_generic("term_predict", "term",
  function(term, newdata, ...) {
    if (!is.data.frame(newdata)) {
      stop("'newdata' must be a data frame.", call. = FALSE)
    }
    S7::S7_dispatch()
  })

# --- printing ---------------------------------------------------------------

#' @title Print a Model Term
#' @name print.model_term
#'
#' @description
#' Prints one line naming the term's class, its label when it has one, and
#' whether it is built. A built term reports how many coefficients its block
#' carries; a specification is marked as one and names the call that would
#' build it. This is the default for every term class, and several classes
#' override it to add the information that only that class carries.
#'
#' @details
#' The two forms are
#'
#' ```
#' <SmoothTerm> 's(x)' built: 4 coefficients
#' <SmoothTerm> 's(x)' (specification; call term_build() with data)
#' ```
#'
#' The class name is the S7 class's own, the label comes from the `label`
#' property and is omitted when empty, and the count is `ncol(x@X)`.
#'
#' A structural term has no `X` to count, so the line shows the word `built` and
#' no count. Every shipped structural term registers a method of its own, so
#' this line is reached only by a structural class written outside the package.
#'
#' The classes that override this add something of their own:
#' [print.PenalizedTerm()] names the penalty and the standardization,
#' [gas()] and [regime()] report their parameters, [nl()] reports its parameters
#' and its derivative mode, and [seg()] reports its break-points.
#'
#' @param x A term, built or not.
#' @param ... Unused, and accepted so that the signature matches [print()]'s.
#'
#' @return `x`, invisibly. The function is called for the line that it writes.
#'
#' @seealso [term_is_built()] for the predicate on which it branches,
#'   [term_coef_names()] for the names behind the count.
#'
#' @examples
#' d <- data.frame(x = 1:4)
#'
#' # A specification prints first, and the same term built prints second.
#' linpar(~ x)
#' term_build(linpar(~ x), d)
#'
#' # The label is shown when there is one.
#' s(x, basis7::bspline_smooth(k = 5))
#' linpar(~ x, label = "lin")
#'
#' @keywords internal
S7::method(print, model_term) <- function(x, ...) {
  cls <- attr(S7::S7_class(x), "name")
  lab <- if (nzchar(x@label)) sprintf(" '%s'", x@label) else ""
  # A structural term has no `X` to count, so the coefficient line is asked
  # of the branch that has one. Every shipped structural term registers a
  # print method of its own; one written outside the package arrives here.
  if (S7::S7_inherits(x, structural_term)) {
    cat(sprintf("<%s>%s %s\n", cls, lab,
                if (term_is_built(x)) "built"
                else "(specification; call term_build() with data)"))
  } else if (term_is_built(x)) {
    cat(sprintf("<%s>%s built: %d coefficient%s\n", cls, lab,
                ncol(x@X), if (ncol(x@X) == 1L) "" else "s"))
  } else {
    cat(sprintf("<%s>%s (specification; call term_build() with data)\n",
                cls, lab))
  }
  invisible(x)
}
