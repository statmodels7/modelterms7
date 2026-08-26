#' @include term_classes.R
NULL

#' @title Build a Term on Data
#'
#' @description
#' Turns a term specification into a built term: an additive term computes its
#' design block from `data`, assigns the coefficient names and records the
#' blueprint that will reproduce the mapping on other rows; a structural term
#' records whatever its recursion needs, its grouping and its ordering. The
#' returned object is a copy of the specification with those properties filled,
#' and the specification is unchanged.
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
#' Building is what produces \eqn{X_t} from the data and attaches
#' \eqn{\rho_t}. [term_matrix()] then reads the block, [term_penalties()] the
#' penalties, and [term_predict()] reproduces \eqn{X_t} at other rows through
#' the blueprint.
#'
#' A structural term has no such block. [gas()] records the group each row
#' belongs to and its place in that group's series; [regime()] and the marginal
#' break-point terms record what their forward recursion reads. They report
#' themselves through [term_filter()] or [term_loglik()], and [term_matrix()]
#' has no method for them.
#'
#' # Building twice, and building on other data
#'
#' Building is not idempotent in general: a term built again on new data
#' re-derives its factor levels and its knots from those rows. That is what
#' [term_predict()] exists to avoid, and what [check_term()]'s subset check
#' tests for. Build once, predict thereafter.
#'
#' # The default
#'
#' There is one, on [model_term()], and it throws
#' `"the term class 'X' does not implement term_build()."`, naming the class,
#' so a term class that supplies nothing else says so clearly. It covers both
#' branches: every shipped structural term registers a method of its own, and
#' a structural class written elsewhere that does not is named the same way.
#'
#' @param term An object inheriting from [model_term()], built or not. A built
#'   term is rebuilt.
#' @param data A data frame carrying every variable the term names. Anything
#'   else throws `"'data' must be a data frame."` from the generic, before
#'   dispatch.
#' @param ... Passed to methods. No shipped method reads anything here.
#'
#' @return A built term of the same class as `term`. For an additive term the
#'   `X`, `coef_names` and `blueprint` properties are filled; for a structural
#'   term the `blueprint` is. [term_is_built()] is `TRUE` either way, reading
#'   whichever of the two the branch fills.
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
#' # A class that implements nothing is told which class it is.
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
#' `TRUE` for a term that [term_build()] has filled, `FALSE` for a bare
#' specification. It is the test [term_matrix()], [term_npar()],
#' [term_coef_names()], [term_predict()] and [plot()] apply before reading a
#' design block, so it decides which error a caller gets from those.
#'
#' @details
#' The two branches record being built in different places, so the predicate
#' asks each about its own. An additive term is built when it has coefficient
#' names, `length(term@coef_names) > 0L`; a structural term contributes no
#' design columns and so has none to count, and is built when its blueprint is
#' filled, `length(term@blueprint) > 0L`.
#'
#' **A built additive term with no columns would answer `FALSE`.** No shipped
#' constructor produces one: [linpar()] with an empty formula fails in the
#' class validator before it gets here.
#'
#' **`blueprint` is asked for rather than assumed.** It is declared on
#' [additive_term()] and on each of the three shipped structural classes, and
#' not on the abstract [structural_term()], so a structural class written
#' outside the package need not carry one. Where it does not, the answer is
#' `FALSE` rather than an error, this predicate promising a logical for
#' anything inheriting from [model_term()].
#'
#' @param term An object inheriting from [model_term()]. Anything else throws
#'   `"'term' must inherit from 'model_term'."`.
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
#' # It is what the accessors test, so it predicts the error.
#' try(term_matrix(linpar(~ x)))
#'
#' # A structural term carries no design block, so it is asked about the
#' # blueprint instead, and the answer is the same question either way.
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
  # in its coefficient names, a structural one in its blueprint, having no
  # design block to name coefficients of.
  #
  # The property is ASKED FOR rather than assumed. `blueprint` is declared on
  # additive_term and on each of the three shipped structural classes, and
  # not on the abstract structural_term, so a structural class written
  # outside the package need not carry one; reading it there would raise
  # S7's "Can't find property" from inside a predicate whose own guard
  # promises a logical for any model_term.
  if (S7::S7_inherits(term, structural_term)) {
    return("blueprint" %in% S7::prop_names(term) &&
             length(term@blueprint) > 0L)
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
#' Returns the \eqn{n \times k} block a built additive term contributes to the
#' linear predictor, with the term's coefficient names as column names and one
#' row per observation of the data it was built on. The block is the term's `X`
#' property, returned as it is stored.
#'
#' @details
#' The block is **not necessarily a base matrix**. [random()] builds a grouping
#' indicator as a `dgCMatrix`, since a row belongs to one group and the density
#' is \eqn{1/m}, and [linpar()], the penalized terms and a smooth with a factor
#' `by` build sparse when asked. Code that reads a block therefore tests
#' `is.matrix(x) && is.numeric(x)` or a two-dimensional S4 object;
#' `is.matrix()` alone is `FALSE` for every \pkg{Matrix} class and would reject
#' a block for being economical.
#'
#' For a term whose block moves with its own coefficients, [nl()] and [seg()],
#' this returns the block at the coefficients last committed by
#' [term_refresh()]. [term_jacobian_block()] says whether that block is a
#' Jacobian or a frozen working linearization.
#'
#' A structural term contributes no block and registers no method, so
#' `term_matrix()` on one stops with S7's method-not-found error.
#'
#' @param term A built additive term (see [term_build()]). A specification
#'   throws
#'   `"the term has not been built; call term_build(term, data) first."`.
#' @param ... Passed to methods. No shipped method reads anything here.
#'
#' @return The design block: a numeric matrix, or a two-dimensional
#'   \pkg{Matrix} object where the term built one, with `nrow` the number of
#'   observations the term was built on and `ncol` equal to [term_npar()].
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
#' # It answers only for a penalty over the whole block
#'
#' A term whose penalty reaches part of its parameters returns `NULL` here and
#' declares that penalty through [term_penalties()], which names the parameters
#' it covers. [seg()] penalizes its changes of slope and leaves the linear
#' effect and the break-points free; [nl()] and [gas()] carry the penalties of
#' whatever sub-terms develop their own parameters. Reporting one of those here
#' would say it covers the block, so this generic answers only where that is
#' true.
#'
#' [term_penalties()] is therefore the general question, and it is the one a
#' fitting layer asks. This one is the convenience for the common case.
#'
#' # A specification carries no penalty
#'
#' The penalty is attached at [term_build()], its width being the number of
#' columns the data produce, so `term_penalty(ridge(~ x))` is `NULL` and
#' `term_penalty(term_build(ridge(~ x), d))` is the quadratic penalty. The same
#' holds for [term_penalties()] and, through it, for [term_smooth()].
#'
#' The one method is registered on [additive_term()] and reads the `penalty`
#' property. A structural term has no such property, so `term_penalty()` on one
#' stops with S7's method-not-found error.
#'
#' @param term An object inheriting from [additive_term()], built or not.
#' @param ... Passed to methods. No shipped method reads anything here.
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
#' # A term whose penalty covers part of its parameters answers NULL here
#' # and declares it through term_penalties() instead.
#' sb <- term_build(seg(x, npsi = 1), data.frame(x = sort(runif(50, 0, 10))))
#' term_penalty(sb)
#'
#' @export
#' @aliases term_penalty.additive_term
term_penalty <- S7::new_generic("term_penalty", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_penalty, additive_term) <- function(term, ...) {
  term@penalty
}

#' @title Every Penalty a Term Carries
#'
#' @description
#' Returns what a term declares it wants penalized: a list of entries, each
#' naming a subset of the term's own parameters, the penalty over them, and
#' whatever the caller fixed about that penalty's hyperparameters. This is the
#' enumeration a fitting layer runs over, and it is what [term_smooth()] and
#' [edf()] read.
#'
#' @details
#' # Two ways it generalizes [term_penalty()]
#'
#' **A term may carry more than one penalty**, over different parameters of its
#' own. A panel model with a population value and a departure per group wants
#' the population value free and the departures shrunk, which is one penalty
#' over part of the parameters and none over the rest. `nl()` developing two of
#' its parameters by two different penalized sub-terms carries two entries.
#'
#' **The parameters need not be coefficients of a design block.** The
#' persistence of a score-driven term, the nonlinear parameters of [nl()] and
#' the break-point of [seg()] are parameters of the term and of nothing else,
#' and all a penalty needs from them is a vector of numbers and their
#' positions. For a structural term `index` gives positions in
#' [term_params()]; for an additive one, columns of the block.
#'
#' # The base method, and the name of an entry
#'
#' The method on [model_term()] answers from [term_penalty()], so a term
#' carrying one penalty over its whole block needs no method of its own. Its
#' single entry is named `""`, meaning the whole term.
#'
#' A name is unique **within** the term and is not the term's own name. Two
#' `ridge()` terms in one formula are two terms with their own
#' hyperparameters, and it is the caller who knows what it called each of them;
#' [statmodels7::statmod()] composes a key as `term` or `term::entry`.
#'
#' The list itself is **unnamed**: read `e$name`, not `names(entries)`.
#'
#' # A specification carries no entries
#'
#' The penalty is attached at [term_build()], so `term_penalties()` on an
#' unbuilt penalized term is an empty list. That is why [term_smooth()], which
#' reads this, answers `TRUE` for `lasso(~ x)` and `FALSE` once it is built.
#'
#' @param term A built term. An unbuilt one is accepted and reports what it has,
#'   which is usually nothing.
#' @param ... Passed to methods. No shipped method reads anything here.
#'
#' @return An unnamed list, possibly empty, one element per declared penalty.
#'   Each element is a list with
#'   \describe{
#'     \item{`name`}{a character label unique within the term, `""` for a
#'       penalty covering the whole of it.}
#'     \item{`index`}{integer positions among the term's own parameters:
#'       columns of the block for an additive term, positions in
#'       [term_params()] for a structural one.}
#'     \item{`penalty`}{a \pkg{penalties7} penalty over exactly those
#'       parameters, so `penalty@n_coef` equals `length(index)`.}
#'     \item{`fixed`}{the hyperparameters the caller held, a named list, empty
#'       when all of them are estimated.}
#'     \item{`n_values`, `values`, `min_ratio`, `search`}{what the caller said
#'       about the path over this penalty's hyperparameters, from the term's
#'       properties of the same names.}
#'   }
#'
#' @seealso [term_penalty()] for the single-penalty case, [term_components()]
#'   for how a term's columns divide among its own parameters, [term_hyper()]
#'   for the held values alone, and [edf()] for what the entries cost.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = sort(runif(50, 0, 10)), g = factor(rep(c("a", "b"), 25)))
#'
#' # One penalty over the whole block: one entry, named with the empty string.
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
  list(list(name = "", index = seq_len(term_npar(term)), penalty = pen,
            fixed = term@hyper, n_values = term@grid,
            values = term@values, min_ratio = term@min_ratio,
            search = term@search))
}

#' @title How a Term's Columns Divide Among Its Own Parameters
#'
#' @description
#' Returns one entry per parameter a term is written in, saying which of its
#' columns belong to that parameter and which sub-terms develop it. A term
#' whose columns are one block with one meaning answers with an empty list.
#'
#' @details
#' # Why the term has to say it
#'
#' A term may be written in parameters of its own: [nl()] in the parameters of
#' \eqn{f}, [seg()] in a slope, a change and a break-point. Any of them may be
#' developed over covariates, and the block then carries several groups of
#' columns meaning different things. Only the term knows which group is which. A
#' coefficient name is built for a reader, and recovering the division by
#' parsing one back is the shape of mistake this package avoids everywhere
#' else.
#'
#' # A parameter may have several sub-terms, of different kinds
#'
#' `seg(x, psi ~ random(~ 1 | id))` develops the break-point with an
#' unpenalized intercept and a random block, so that component's `subs` has two
#' entries and only the second carries a penalty. A consumer reporting a
#' component reports a sequence.
#'
#' `sub_index` splits `index` among those sub-terms, in the order the block
#' binds them, which is [component_sub_index()] applied to their coefficient
#' counts.
#'
#' # Structural terms
#'
#' A structural term contributes no design columns, and `index` gives positions
#' in [term_params()]: the vector its state, its readable quantities and its
#' variance matrix are all indexed by. In both branches the field names the
#' term's own coefficients.
#'
#' The base method returns an empty list, which is the answer for [linpar()],
#' [s()], [random()] and the five penalized constructors, whose columns are one
#' block with one meaning.
#'
#' @param term A built term. An unbuilt one returns an empty list.
#' @param ... Passed to methods. No shipped method reads anything here.
#'
#' @return A named list, one element per own parameter and named by it, each a
#'   list with
#'   \describe{
#'     \item{`name`}{the parameter's name, the same as the element's.}
#'     \item{`index`}{its columns in the term's block, or its positions in
#'       [term_params()] for a structural term.}
#'     \item{`subs`}{the sub-terms developing it, an empty list where none do.}
#'     \item{`sub_index`}{one integer vector per sub-term, splitting `index`
#'       among them; empty where `subs` is.}
#'   }
#'   An empty list for a term whose columns do not divide.
#'
#' @seealso [component_sub_index()] for the split, [term_penalties()] for the
#'   penalties those sub-terms bring, [term_coef_names()] for the names of the
#'   columns being divided.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = seq(0.2, 3, length.out = 20),
#'                  g = factor(rep(c("a", "b"), 10)))
#' dd$y <- 2 * exp(-1.3 * dd$x)
#'
#' # nl() in two parameters, the first developed over a factor: two
#' # columns for `a` and one for `r`.
#' b <- term_build(nl(~ a * exp(-r * x), a ~ 0 + g, start = list(r = 1.3)), dd)
#' lapply(term_components(b), function(z) z$index)
#' term_coef_names(b)
#'
#' # The sub-terms of the developed parameter, and its columns split
#' # among them.
#' term_components(b)$a$sub_index
#'
#' # A term whose columns are one block answers with nothing.
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
#' the order the block binds them. It is what fills the `sub_index` field of a
#' [term_components()] entry.
#'
#' @details
#' A developed parameter's block is its sub-terms' blocks bound side by side in
#' the order they were given, so the division is their coefficient counts
#' cumulated: with counts \eqn{k_1, \dots, k_m} the \eqn{i}-th sub-term takes
#' `index` at positions \eqn{k_1 + \dots + k_{i-1} + 1} to
#' \eqn{k_1 + \dots + k_i}. The counts come from [term_npar()], so every
#' sub-term must be built.
#'
#' The term computes this rather than a consumer, because it rests on how the
#' block was assembled.
#'
#' @param index An integer vector: the component's columns in the term's block,
#'   as long as the sub-terms' coefficient counts sum to.
#' @param subs A list of built sub-terms developing the parameter, in the order
#'   their blocks were bound.
#'
#' @return A list of integer vectors, one per sub-term, partitioning `index` in
#'   order. An empty list when `subs` is empty.
#'
#' @seealso [term_components()], the only caller; [term_npar()] for the counts
#'   it divides by.
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
#' How many parameters of its own a built term carries: the number of columns
#' of the design block for an additive term, and the number of entries of
#' [term_params()] for a structural one, which has no block. It is the length
#' of the vector a [term_penalties()] entry indexes into, the length
#' [edf()] measures against, and the length a fit reserves for the term.
#'
#' @details
#' The two methods are the two branches. On [additive_term()] it is
#' `ncol(term@X)`, so it equals `length(term_coef_names(term))` and a
#' specification throws. On [structural_term()] it is
#' `length(term_params(term))`, which is the count of the term's own
#' parameters after any subformula has expanded: `gas(p = 1, q = 1)` has three,
#' and `gas(p = 1, q = 1, omega ~ z)` has four, the level's intercept and slope
#' in place of the level.
#'
#' @param term A built term (see [term_build()]). An unbuilt additive term
#'   throws
#'   `"the term has not been built; call term_build(term, data) first."`.
#' @param ... Passed to methods. No shipped method reads anything here.
#'
#' @return A single whole number.
#'
#' @seealso [term_coef_names()] and [term_params()] for the names behind the
#'   count, [term_matrix()] for the block, [edf()] for what the term spends of
#'   it.
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
#' The coefficients a built term asks to be started at, one per column of its
#' block. The base method returns zero everywhere. A term whose block is a
#' fixed design wants exactly that: the objective is convex in those
#' coordinates and the fit reaches the same optimum from anywhere.
#'
#' @details
#' # Why any term needs a start of its own
#'
#' A term that recomputes its block from its coefficients ([term_refresh()]) is
#' the case this exists for, because zero there is degenerate. In [jump()] the
#' break-point is read off two coefficients as \eqn{-g_k/\delta_k}, so a vector
#' of zeros puts every break-point at the same clamped position and makes the
#' block singular. In [seg()] the Jacobian column is
#' \eqn{-\gamma_k \mathbb{1}(x > \psi_k)} and vanishes identically at
#' \eqn{\gamma_k = 0}.
#'
#' Those terms return the start [term_build()] computed: unit changes, and the
#' break-points at the interior quantiles of the covariate or at the positions
#' `psi` names. [nl()] returns the starting values of its own parameters
#' carried through their links.
#'
#' Only the term knows what a coefficient of zero means for the block it
#' builds, which is the same reason [term_start()] belongs to a structural
#' term.
#'
#' # What `target` is for
#'
#' `target` is the response carried onto the scale of the predictor the term
#' contributes to. It is the one thing a term cannot work out for itself: the
#' term knows its formula and its charts, and the fitting layer knows the
#' distribution, the link and the equation. [nl()] uses it to estimate its own
#' parameters from the data, over a deterministic grid on each free parameter's
#' chart; every other term ignores it. It is optional, so the default is what
#' every term did before it existed.
#'
#' @param term A built term (see [term_build()]). An unbuilt one throws through
#'   [term_npar()].
#' @param target Optional numeric vector, one value per observation: the
#'   response on the scale of the predictor, `NULL` by default. Supplied only
#'   where `params_interpretation` says the response reads the parameter
#'   directly, so a term in a scale's equation is handed nothing.
#' @param ... Passed to methods. No shipped method reads anything here.
#'
#' @return A numeric vector of length [term_npar()], in the block's column
#'   order.
#'
#' @seealso [term_start()] for a structural term's own parameters,
#'   [term_refresh()] for the terms this exists for, [seg_start()] for the grid
#'   rule behind a break-point start.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(50, 0, 10)))
#'
#' # A fixed design starts at zero.
#' term_coef_start(term_build(linpar(~ x), dd))
#'
#' # A break-point term does not: its block would be singular there.
#' jb <- term_build(jump(x), dd)
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

#' @title Is a Term's Block the Jacobian of Its Contribution?
#'
#' @description
#' `TRUE` when the block a term reports is the exact derivative of its
#' contribution in its own coefficients, `FALSE` when it is a working
#' linearization with quantities frozen at the previous iterate. The base
#' method returns `TRUE`, which a fixed design satisfies trivially and which
#' [nl()] and [seg()] satisfy by construction; [jump()] and [jseg()] answer
#' `FALSE`.
#'
#' @details
#' # What the answer decides
#'
#' Where the block is a Jacobian, a scoring step on it is a Gauss-Newton step,
#' a line search on the model's own objective is licensed, and the term is
#' fitted inside the same system as everything else.
#'
#' Where it is a working linearization the fixed-point iteration of Fasola,
#' Muggeo and Kuchenhoff (2018) is not a descent method on the model's
#' objective. Its early steps go uphill on purpose, under a scaling factor that
#' anneals, so forcing a sufficient decrease on it stalls the iteration.
#' Such a term is fitted by alternating exact working fits at the frozen block
#' with [term_refresh()], and its convergence is what [term_converged()]
#' answers instead of a score.
#'
#' The two discontinuous constructions are the ones that answer `FALSE`. Their
#' weight \eqn{W = 1/(2\lvert\tilde x - \psi\rvert)} is held at the previous
#' break-point, and the position is read off two coefficients rather than being
#' one.
#'
#' # Smoothing changes the answer
#'
#' `jump(x, smoothed = ...)` replaces the step by a smooth surrogate, and the
#' break-point becomes an ordinary parameter with a true Jacobian, so a
#' smoothed break-point term answers `TRUE`. That is how a fitting layer routes
#' it without a special case.
#'
#' @param term A term, built or not: the answer is a property of the
#'   construction.
#' @param ... Passed to methods. No shipped method reads anything here.
#'
#' @return A single logical.
#'
#' @references
#' Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
#' iterative algorithm for change-point detection in abrupt change models.
#' *Computational Statistics*, 33, 997--1015.
#'
#' @seealso [term_refresh()] for the block being recomputed,
#' [term_converged()] for the verdict on such a term, [seg()], [jump()] and
#' [jseg()].
#'
#' @examples
#' # The continuous construction differentiates; the two discontinuous
#' # ones report a working linearization.
#' vapply(list(seg = seg(x), jump = jump(x), jseg = jseg(x)),
#'        term_jacobian_block, logical(1))
#'
#' # A fixed design is a Jacobian trivially, and so is nl().
#' term_jacobian_block(linpar(~ x))
#' term_jacobian_block(nl(~ a * x, start = list(a = 1)))
#'
#' # Smoothing the step makes the break-point an ordinary parameter.
#' term_jacobian_block(jump(x, smoothed = penalties7::smooth_probit()))
#'
#' @export
#' @aliases term_jacobian_block.model_term
term_jacobian_block <- S7::new_generic("term_jacobian_block", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_jacobian_block, model_term) <- function(term, ...) TRUE

#' @title Coefficient Names of a Built Term
#'
#' @description
#' The names of a built additive term's coefficients, one per column of its
#' block and in the block's own order. They are the names a fit reports, so a
#' coefficient table is readable without knowing which term produced which row.
#'
#' @details
#' A term with a non-empty `label` prefixes every name with it and a dot, so
#' `linpar(~ x, label = "lin")` gives `lin.(Intercept)` and `lin.x`. The five
#' penalized constructors set the label to the constructor's own name by
#' default, which is why a ridge over `x` and `g` reads `ridge.x`, `ridge.ga`
#' and so on; [linpar()] sets none, so its names are
#' [stats::model.matrix()]'s unchanged.
#'
#' The names are assigned at build time and recorded, so they are also what
#' [term_predict()] labels its columns with at other rows. Uniqueness is not
#' enforced here; [check_term()] checks it.
#'
#' @param term A built additive term (see [term_build()]). A specification
#'   throws
#'   `"the term has not been built; call term_build(term, data) first."`.
#' @param ... Passed to methods. No shipped method reads anything here.
#'
#' @return A character vector of length [term_npar()], in column order.
#'
#' @seealso [term_npar()] for the count, [term_matrix()] for the block they
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
#' `TRUE` when the term's contribution to the penalized objective is
#' differentiable in the coefficients. A fitting layer reads it to split the
#' coefficient vector into the block a classical optimizer handles and the
#' block that needs a proximal step or a coordinate descent.
#'
#' @details
#' # The answer comes from the penalties
#'
#' The term does not declare it. Every entry of [term_penalties()] is asked for
#' its kink set through [penalties7::penalty_kinks()], at a probe value inside
#' each hyperparameter's bounds, and the answer is `FALSE` as soon as one
#' entry reports a point. So an unpenalized term is smooth, a ridge or a
#' Gaussian prior is smooth, and lasso, SCAD, MCP and the elastic net are not.
#' A term cannot disagree with its own penalties.
#'
#' The probe is any admissible value, the kink set being structural: the
#' midpoint of a bounded interval, one step inside a half-bounded one, zero
#' where the interval is the whole line.
#'
#' # It answers for the whole term
#'
#' A term carrying a penalty over part of its parameters and none over the rest
#' answers for the part. `nl(~ a * exp(-r * x), a ~ 0 + lasso(~ g))` is not
#' smooth: the coefficients developing `a` sit at a kink, although `r` is
#' unpenalized.
#'
#' # A specification is always smooth
#'
#' The penalty is attached at [term_build()], so `term_penalties()` on an
#' unbuilt term is empty and this answers `TRUE` whatever penalty the term will
#' carry. `term_smooth(lasso(~ x))` is `TRUE` and
#' `term_smooth(term_build(lasso(~ x), d))` is `FALSE`. Ask a built term.
#'
#' @param term An object inheriting from [model_term()]. Build it first, or the
#'   answer is `TRUE` by default.
#' @param ... Passed to methods. No shipped method reads anything here.
#'
#' @return A single logical, never `NA`.
#'
#' @seealso [term_penalties()] for the entries it runs over,
#'   [penalties7::penalty_kinks()] for the set it reads, [edf()], which uses
#'   the same split to count degrees of freedom.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(30), g = factor(rep(c("a", "b"), 15)))
#'
#' # Unpenalized and quadratically penalized terms are smooth.
#' vapply(list(linpar(~ x), ridge(~ x), s(x, k = 5), random(~ 1 | g)),
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
#' Applies a built term's recorded mapping to new rows, returning the block the
#' term would have produced had those rows been in the data it was built on.
#' Factor levels, contrasts, spline knots, a Demmler-Reinsch reparametrization
#' and the spreads a standardization used all come from the blueprint and are
#' reused. A factor level the blueprint does not know is rejected.
#'
#' @details
#' # The identity that makes it useful
#'
#' The block returned is \eqn{\tilde{X}_t} such that
#' \eqn{\tilde{\eta} = \tilde{X}_t \beta_t} is the term's contribution at the
#' new rows, at the coefficients the model already carries. The identity holds
#' because the mapping is reused. A rebuilt encoding gives a block
#' of the same shape multiplying the same coefficients and meaning something
#' else: a factor whose new rows omit a level loses a column, and a basis
#' rebuilt on a narrower range is a different set of functions.
#'
#' The difference is not small. On a smooth of 60 points over \eqn{[0, 1]},
#' predicting on the first 20 rows agrees with those rows of the original block
#' exactly, while rebuilding the term on them differs by 2.33 in the same
#' units. [check_term()]'s subset check is exactly this comparison.
#'
#' Predicting on the fitting data returns the block itself, so
#' `term_predict(b, data)` and `term_matrix(b)` agree to the last bit.
#'
#' # Which terms have a method
#'
#' Six do: [linpar()], the penalized terms, [random()], the smooths, [nl()] and
#' the break-point terms. A structural term contributes no block and has no
#' method, so `term_predict()` on one stops with S7's method-not-found error;
#' [term_continue()] is the corresponding operation there.
#'
#' @param term A built additive term (see [term_build()]). A specification
#'   throws
#'   `"the term has not been built; call term_build(term, data) first."`.
#' @param newdata A data frame carrying every variable the term names, with any
#'   number of rows. Anything else throws `"'newdata' must be a data frame."`
#'   from the generic, before dispatch. A factor here need carry only the
#'   levels its own rows use; the rest come from the blueprint.
#' @param ... Passed to methods. No shipped method reads anything here.
#'
#' @return A block of `nrow(newdata)` rows and [term_npar()] columns, in the
#'   same storage the term built: a numeric matrix, or a \pkg{Matrix} object
#'   where the block is sparse.
#'
#' @seealso [term_matrix()] for the block on the fitting data, [term_build()]
#'   for what records the blueprint, [check_term()] for the check this identity
#'   is the subject of, and [term_continue()] for a structural term.
#'
#' @examples
#' d <- data.frame(x = 1:6, g = factor(rep(c("a", "b", "c"), 2)))
#' b <- term_build(linpar(~ x + g), d)
#'
#' # New rows, the same mapping.
#' term_predict(b, data.frame(x = c(0.5, 2.5), g = factor(c("a", "c"))))
#'
#' # On the fitting data it returns the block itself.
#' all.equal(term_predict(b, d), term_matrix(b))
#'
#' # A subset that drops a level keeps the blueprint's columns, where a
#' # rebuild would lose one.
#' nd <- droplevels(d[d$g != "c", ])
#' dim(term_predict(b, nd))
#' dim(model.matrix(~ x + g, nd))
#'
#' # A basis is not replaced on the narrower range: reapplying agrees with
#' # the original rows exactly, rebuilding does not.
#' d2  <- data.frame(x = seq(0, 1, length.out = 60))
#' bs  <- term_build(s(x, k = 6), d2)
#' X   <- term_matrix(bs)
#' sub <- 1:20
#' max(abs(term_predict(bs, d2[sub, , drop = FALSE]) - X[sub, ]))
#' max(abs(term_matrix(term_build(s(x, k = 6), d2[sub, , drop = FALSE])) -
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
#' carries; a specification says so and names the call that would build it.
#' This is the default for every term class, and several classes override it to
#' add what they alone carry.
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
#' A structural term has no `X` to count, so it gets the word `built` and no
#' count. Every shipped structural term registers a method of its own, so this
#' line is reached only by a structural class written outside the package.
#'
#' The classes that override this add something of their own:
#' [print.PenalizedTerm()] names the penalty and the standardization,
#' [gas()]'s and [regime()]'s report their parameters, [nl()]'s its formula and
#' [seg()]'s its break-points.
#'
#' @param x A term, built or not.
#' @param ... Unused, and accepted so that the signature matches [print()]'s.
#'
#' @return `x`, invisibly. Called for the line it writes.
#'
#' @seealso [term_is_built()] for the predicate it branches on,
#'   [term_coef_names()] for the names behind the count.
#'
#' @examples
#' d <- data.frame(x = 1:4)
#'
#' # A specification, and the same term built.
#' linpar(~ x)
#' term_build(linpar(~ x), d)
#'
#' # The label is shown when there is one.
#' s(x, k = 5)
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
