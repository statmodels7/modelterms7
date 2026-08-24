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
#' # The two defaults
#'
#' `term_build.model_term` throws
#' `"the term class 'X' does not implement term_build()."`, naming the class,
#' so a term class that supplies nothing else says so clearly.
#'
#' `term_build.structural_term` throws
#' `"structural terms are reserved for a later release; none is implemented yet."`
#' Every structural term that ships overrides it, so the message is reachable
#' only from a structural class written elsewhere.
#'
#' @param term An object inheriting from [model_term()], built or not. A built
#'   term is rebuilt.
#' @param data A data frame carrying every variable the term names. Anything
#'   else throws `"'data' must be a data frame."` from the generic, before
#'   dispatch.
#' @param ... Passed to methods. No shipped method reads anything here.
#'
#' @return A built term of the same class as `term`. For an additive term the
#'   `X`, `coef_names` and `blueprint` properties are filled and
#'   [term_is_built()] is `TRUE`; for a structural term the `blueprint` is
#'   filled and [term_is_built()] stays `FALSE`, that predicate testing for a
#'   design block.
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

S7::method(term_build, structural_term) <- function(term, data, ...) {
  stop("structural terms are reserved for a later release; none is implemented yet.",
       call. = FALSE)
}

#' @title Whether a Term Carries a Design Block
#'
#' @description
#' `TRUE` for an additive term that [term_build()] has filled, `FALSE` for a
#' bare specification. It is the test [term_matrix()], [term_npar()],
#' [term_coef_names()], [term_predict()] and [plot()] apply before reading a
#' block, so it decides which error a caller gets from those.
#'
#' @details
#' The predicate is `S7::S7_inherits(term, additive_term)` together with
#' `length(term@coef_names) > 0L`, so it asks whether the term has coefficients
#' to name. Two consequences follow, and both matter to a caller writing
#' against the class rather than against one term.
#'
#' **A structural term answers `FALSE` whether or not it is built.** [gas()]
#' and [regime()] contribute no design columns, so they have no coefficient
#' names to count; a built one has its blueprint filled and answers every
#' generic on its own branch. The test covering both branches is
#' `length(term@blueprint) > 0L`.
#'
#' **A built additive term with no columns would also answer `FALSE`.** No
#' shipped constructor produces one: [linpar()] with an empty formula fails in
#' the class validator before it gets here.
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
#' # A structural term has no block, so it answers FALSE once built.
#' gb <- term_build(gas(p = 1, q = 1), data.frame(y = rnorm(20)))
#' c(predicate = term_is_built(gb), blueprint_filled = length(gb@blueprint) > 0L)
#'
#' @seealso [term_build()], [term_matrix()], [term_coef_names()], [term_npar()]
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
#' The penalty attached to the whole of the term's coefficients, or
#' `NULL` when there is none. The hyperparameters, their bounds and
#' links, and every derivative in the coefficients and the hyperparameters
#' are the penalty object's, not the term's.
#'
#' @details
#' A term whose penalty reaches only part of its parameters returns
#' `NULL` here and declares that penalty through
#' [term_penalties()], which names the parameters it covers:
#' [seg()] penalizes the changes and not the linear effect or the
#' break-points, and the developments of [nl()] and
#' [gas()] carry their sub-terms' penalties. Reading
#' a partial penalty here would say that it covers the block, so the
#' question this generic asks is answered only where the answer is the whole
#' of it.
#'
#' @param term An object inheriting from class [additive_term()].
#' @param ... Passed to methods.
#'
#' @return A penalty object, or `NULL`.
#'
#' @examples
#' term_penalty(linpar(~x))
#'
#' @seealso [term_penalties()], [term_smooth()], [edf()]
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
#' [term_penalty()] answers for the common case, one penalty over the
#' whole of a term's design block, and this generalizes it in two directions a
#' model layer needs.
#'
#' A term may carry **more than one** penalty, over different parameters of
#' its own. A panel model with a population value and a departure per group
#' wants the population value free and the departures shrunk, which is one
#' penalty over part of the parameters and none over the rest.
#'
#' The parameters need **not be coefficients of a design block**. The
#' persistence of a score-driven term, the nonlinear parameters of
#' [nl()], the break-point of [seg()] are parameters of the
#' term and nothing else, and everything a penalty needs from them is a vector
#' of numbers and their positions.
#'
#' The base method answers from `term_penalty()`, so a term that carries
#' one penalty over its whole block -- every term shipped here -- needs no
#' method of its own and behaves exactly as before. Its single entry is named
#' with the empty string, meaning the whole term, so a caller that keys the
#' hyperparameters by term name keys them exactly as it did.
#'
#' @param term A built term.
#' @param ... Passed to methods.
#'
#' @return A list, possibly empty. Each entry has `name` (a label unique
#'   WITHIN the term, empty for a penalty over the whole of it), `index`
#'   (positions among the term's parameters) and `penalty` (a
#'   \pkg{penalties7} object). The name is not the term's: two `ridge()`
#'   terms in one formula are two terms with their own hyperparameters, and it
#'   is the caller that knows what it called each one.
#'
#' @examples
#' term_penalties(term_build(ridge(~x), data.frame(x = rnorm(20))))
#' term_penalties(term_build(linpar(~x), data.frame(x = rnorm(20))))
#'
#' @seealso [term_penalty()], [term_npar()]
#' @export
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
#' One entry per parameter the term is written in, saying which columns of its
#' block belong to that parameter and which sub-terms develop it, or an empty
#' list for a term whose columns answer to nothing above them.
#'
#' @details
#' A term written in parameters of its own -- [nl()] in the
#' parameters of \eqn{f}, [seg()] in a slope, a change and a
#' break-point -- may develop any of them over covariates, and then its block
#' carries several groups of columns that mean different things. What divides
#' them is the TERM's answer and cannot be recovered from the coefficient
#' names: a name is built for a reader and parsing one back is the shape of
#' mistake this package avoids everywhere else.
#'
#' A parameter may be developed by SEVERAL sub-terms at once, and they need
#' not be of one kind: `seg(x, psi ~ random(~1 | id))` develops the
#' break-point with an unpenalized intercept AND a random block, so the
#' component's `subs` has two entries and only the second carries a
#' penalty. A consumer that reports a component therefore reports a sequence
#' and not a single kind.
#'
#' A STRUCTURAL term contributes no design columns, and there `index`
#' gives positions in [term_params()] instead: the vector its
#' state, its readable quantities and its variance matrix are all indexed by.
#' In both cases the field names the term's own coefficients.
#'
#' The base method returns an empty list, which says that the term's columns
#' are its own and divide no further. That is the honest answer for
#' [linpar()], [s()], [random()] and the
#' penalized constructors, whose columns are one block with one meaning.
#'
#' @param term A built term.
#' @param ... Passed to methods.
#'
#' @return A list, one entry per own parameter, each with `name` (the
#'   parameter), `index` (its columns in the term's block),
#'   `subs` (the sub-terms developing it, empty where there are none)
#'   and `sub_index` (the columns belonging to each of those
#'   sub-terms). Empty for a term whose columns do not divide.
#'
#' @examples
#' dd <- data.frame(x = seq(0.2, 3, length.out = 20),
#'                  g = factor(rep(c("a", "b"), 10)))
#' dd$y <- 2 * exp(-1.3 * dd$x)
#' b <- term_build(nl(~ a * exp(-r * x), a ~ 0 + g, start = list(r = 1.3)), dd)
#' lapply(term_components(b), function(z) z$index)
#'
#' # a term whose columns are one block answers with nothing
#' term_components(term_build(linpar(~ x), dd))
#'
#' @seealso [term_penalties()], [term_coef_names()]
#' @export
term_components <- S7::new_generic("term_components", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_components, model_term) <- function(term, ...) list()
#' The Columns of a Component That Belong to Each of Its Sub-Terms
#'
#' @description
#' Splits a developed parameter's columns among the sub-terms developing it,
#' in the order the block binds them.
#'
#' @details
#' A developed parameter's block is its sub-terms' blocks bound side by side
#' in the order they were given, so the division is their coefficient counts
#' cumulated. It is computed by the term rather than left to a consumer
#' because it rests on how the block is assembled, which is the term's
#' business and not something a name can be parsed for.
#'
#' @param index The component's columns in the term's block.
#' @param subs The sub-terms developing the parameter.
#'
#' @return A list of integer vectors, one per sub-term, empty where there are
#'   no sub-terms.
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
#' The coefficients a built term asks to be started at, one per column of
#' its block. The base method returns zero everywhere, which is what a term
#' whose block is a fixed design wants: the fit reaches the same optimum
#' from anywhere, the objective being convex in those coordinates.
#'
#' @details
#' A term that recomputes its block from its coefficients
#' ([term_refresh()]) is the case this exists for, because zero
#' is not a neutral point there but a degenerate one. In
#' [jump()] the break-point is read off two coefficients as
#' \eqn{-g_k/\delta_k}, so a vector of zeros puts every break-point at the
#' same clamped position and makes the block singular; in
#' [seg()] the Jacobian column is \eqn{-\gamma_k\,\mathbb{1}(x >
#' \psi_k)} and vanishes identically. Those terms return the start
#' [term_build()] computed -- unit changes and the break-points
#' at the interior quantiles of the covariate, or the positions
#' `psi` names -- and [nl()] returns the starting values of
#' its own parameters carried through their links.
#'
#' The value belongs to the term for the reason [term_start()]
#' records for a structural one: only the term knows what a coefficient of
#' zero means for the block it builds.
#'
#' `target` is the response carried onto the scale of the predictor the
#' term contributes to, which is what a term needs to estimate parameters of
#' its own from the data and is the one thing it cannot work out for itself:
#' the term knows its formula and its charts, the fitting layer knows the
#' distribution, the link and the equation. It is optional, and a term that
#' has no use for it ignores it, so the default is the behaviour every term
#' had before it existed.
#'
#' @param term A built term (see [term_build()]).
#' @param target Optional numeric vector, one value per observation: the
#'   response on the scale of the predictor. [nl()] uses it to
#'   estimate its own parameters; every other term ignores it.
#' @param ... Passed to methods.
#'
#' @return A numeric vector of length [term_npar()].
#'
#' @examples
#' dd <- data.frame(x = sort(runif(50, 0, 10)))
#' term_coef_start(term_build(linpar(~x), dd))
#' term_coef_start(term_build(jump(x), dd))
#'
#' @seealso [term_start()], [term_refresh()]
#' @export
term_coef_start <- S7::new_generic("term_coef_start", "term",
  function(term, target = NULL, ...) S7::S7_dispatch())

S7::method(term_coef_start, model_term) <- function(term, target = NULL, ...) {
  numeric(term_npar(term))
}

#' @title Is a Term's Block the Jacobian of Its Contribution?
#'
#' @description
#' `TRUE` when the design block a term reports is the exact derivative
#' of its contribution in its own coefficients, `FALSE` when it is a
#' working linearization with quantities frozen at the previous iterate. The
#' base method returns `TRUE`, which is what a fixed design satisfies
#' trivially and what [nl()] and [seg()] satisfy by
#' construction.
#'
#' @details
#' The distinction decides how a fitting layer may treat the block. Where
#' the block is a Jacobian, a scoring step on it is a Gauss--Newton step and
#' a line search on the model's own objective is licensed, so the term can
#' be fitted inside the same system as everything else. Where it is a
#' working linearization -- [jump()] and [jseg()], whose
#' weight \eqn{W = 1/(2\lvert \tilde x - \psi\rvert)} is held at the
#' previous break-point and whose position is read off two coefficients --
#' the fixed-point iteration of \cite{fasola2018} is not a descent method on
#' the model's objective, and forcing a sufficient decrease on it stalls the
#' iteration. Such a term is fitted by alternating exact working fits at the
#' frozen block with [term_refresh()], and its convergence is what
#' [term_converged()] answers rather than a score.
#'
#' @param term A term (built or not; the answer is a property of the
#'   construction).
#' @param ... Passed to methods.
#'
#' @return A single logical.
#'
#' @references
#' Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
#' iterative algorithm for change-point detection in abrupt change
#' models. *Computational Statistics*, 33, 997--1015.
#'
#' @examples
#' term_jacobian_block(seg(x))
#' term_jacobian_block(jump(x))
#'
#' @seealso [term_refresh()], [term_converged()]
#' @export
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
#' differentiable in the coefficients. The answer is read from the penalties
#' rather than declared by the term: an unpenalized term is smooth, and a
#' penalized one is smooth exactly when no penalty it carries declares a
#' kink, so a term cannot disagree with its own penalties. The model layer
#' uses this flag to split the coefficient vector into the block the
#' classical optimizers handle and the block that needs non-smooth
#' strategies.
#'
#' @details
#' The enumeration is [term_penalties()], so a term carrying one
#' penalty over part of its parameters and none over the rest answers for
#' the part: `seg(x, penalty = penalties7::lasso_penalty)` is not
#' smooth, its slope
#' changes sitting at a kink, although its linear effect and its
#' break-points are unpenalized.
#'
#' @param term An object inheriting from class [model_term()].
#' @param ... Passed to methods.
#'
#' @return A logical scalar.
#'
#' @examples
#' term_smooth(linpar(~x))
#'
#' @seealso [term_penalties()], [term_penalty()], [edf()]
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
#' The classes that override this add something of their own:
#' [print.PenalizedTerm()] names the penalty and the standardization,
#' [gas()]'s and [regime()]'s report their parameters, [nl()]'s its formula and
#' [seg()]'s its break-points. A structural term is never reported as built
#' here, [term_is_built()] testing for a design block.
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
  if (term_is_built(x)) {
    cat(sprintf("<%s>%s built: %d coefficient%s\n", cls, lab,
                ncol(x@X), if (ncol(x@X) == 1L) "" else "s"))
  } else {
    cat(sprintf("<%s>%s (specification; call term_build() with data)\n",
                cls, lab))
  }
  invisible(x)
}
