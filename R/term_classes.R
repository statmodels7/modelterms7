#' @title S7 Base Class for Model Terms
#' @name model_term
#'
#' @description
#' The abstract root of the term hierarchy. A term records what a formula
#' names: the recipe turning a data frame into a contribution to the model,
#' together with the metadata a fit reads. A term written in a formula is a
#' **specification**, carrying only what the call said; [term_build()] turns it
#' into a **built** term carrying the design block or the state the recipe
#' produces on given data.
#'
#' [interpret_formula()] recognizes a term by this class. Any call on the right
#' of a formula whose value inherits from `model_term` becomes a term, so a
#' term class defined outside the package works in a formula the day it is
#' written, with nothing to register.
#'
#' @details
#' # The two branches
#'
#' [additive_term()] contributes a block of design columns \eqn{X_j\beta_j} to
#' a linear predictor. [structural_term()] rewrites the likelihood instead:
#' [gas()]'s predictor is a recursion, [regime()]'s contribution is a
#' likelihood mixed over latent states. The branch decides which generics a
#' consumer may call, so [statmodels7::statmod()] routes on it.
#'
#' # What the root carries, and why
#'
#' Beyond `label` the properties are all about a hyperparameter path.
#' **Which hyperparameters are estimated is said by the term and by nothing
#' else**: the term is where the penalty is named, so it is where a held value,
#' a grid size, a written-out set of values, the depth of the path and the way
#' several hyperparameters are combined all belong. An outer criterion is put
#' to every hyperparameter of a model, smooth ones included, and carries none
#' of this.
#'
#' Each of `hyper`, `grid` and `values` is a named list keyed by the penalty's
#' own hyperparameter names, and each is empty by default. `min_ratio` and
#' `search` are single values rather than one per hyperparameter, because only
#' the path over the size of the kink uses a ratio: a bounded hyperparameter is
#' swept over its own interval and a shape over a geometric grid above its
#' lower bound.
#'
#' # It cannot be constructed
#'
#' `model_term` is abstract, as are `additive_term` and `structural_term`.
#' `model_term()` throws
#' `"Can't construct an object from abstract class <model_term>"`. Use it as a
#' parent when writing a term class, and as the test
#' `S7::S7_inherits(x, model_term)` when asking whether an object is a term.
#'
#' A class inheriting from it supplies [term_build()] at the very least: the
#' method registered here throws, naming the class that did not implement it.
#' Most other generics carry a usable default, so a new term class starts from
#' a working object and overrides what it has reason to.
#'
#' @param label A character string prefixed to the term's coefficient names
#'   when non-empty, and used as the title of [plot()] and the tag of
#'   [print()]. `character(0)` and `""` both mean no label.
#' @param hyper The hyperparameters of the term's penalty the caller **held**,
#'   as a named list keyed by the penalty's names. Empty, the default, means
#'   every one of them is estimated. See [term_hyper()].
#' @param grid How many values a path visits for each of the term's
#'   hyperparameters, as a named list of single whole numbers. Empty, the
#'   default, leaves the number to the fitting layer. See [term_grid()].
#' @param values The values a path visits, for each hyperparameter the caller
#'   wrote out, as a named list of numeric vectors. Empty, the default, has the
#'   path build them. See [term_values()].
#' @param min_ratio How far down the path over the size of the kink reaches, as
#'   a fraction of the value that empties the block: one number in \eqn{(0, 1)},
#'   or `numeric(0)` for the fitting layer's own. See [term_path_min()].
#' @param search How the term's own hyperparameters are covered when it has
#'   more than one carrying a kink: `"grid"` for every combination of them,
#'   `"cyclic"` for one at a time, or `character(0)` for the default. See
#'   [term_search()].
#'
#' @return Nothing: the class is abstract and cannot be instantiated. As a type
#'   it is the parent of every term, with the six properties above.
#'
#' @seealso [additive_term()] and [structural_term()] for the two branches,
#'   [term_build()] for turning a specification into a built term,
#'   [interpret_formula()] for how a formula is read into terms, and
#'   [check_term()] for validating one.
#'
#' @examples
#' # Every term inherits from this class, whichever branch it is on.
#' vapply(list(linpar(~ 1), ridge(~ x), s(x, k = 5), gas(p = 1, q = 1)),
#'        function(t) S7::S7_inherits(t, model_term), logical(1))
#'
#' # The branch is what a consumer routes on.
#' c(additive = S7::S7_inherits(s(x, k = 5), additive_term),
#'   structural = S7::S7_inherits(s(x, k = 5), structural_term))
#' c(additive = S7::S7_inherits(gas(p = 1, q = 1), additive_term),
#'   structural = S7::S7_inherits(gas(p = 1, q = 1), structural_term))
#'
#' # The path properties are empty until a caller sets one.
#' r <- ridge(~ x)
#' lengths(list(hyper = r@hyper, grid = r@grid, values = r@values,
#'              min_ratio = r@min_ratio, search = r@search))
#'
#' # Holding lambda puts it in `hyper`, where term_hyper() reads it.
#' term_hyper(ridge(~ x, lambda = 0.5))
#'
#' # Abstract: there is nothing to construct.
#' try(model_term())
#'
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
    # The values themselves, where the caller wrote them out instead of
    # having them built. This is the third state of one argument: NULL has
    # the path construct the grid, one number holds the hyperparameter, and
    # several are the grid, used as given. A term may be in a different
    # state for each of its hyperparameters. See term_values().
    values = S7::new_property(S7::class_list, default = quote(list())),
    # How far DOWN the path over the size of the kink reaches, as a
    # fraction of the value that empties the block. It belongs beside the
    # grid size and for the same reason, and it is one number rather than
    # one per hyperparameter because only that path uses it: a bounded
    # hyperparameter is swept over its own interval and a shape over a
    # geometric grid above its lower bound, where a ratio means nothing.
    min_ratio = S7::new_property(S7::class_numeric,
                                 default = quote(numeric(0))),
    # HOW its own hyperparameters are covered when it has several with a
    # kink: every combination of them, or one at a time. It belongs beside
    # the grid size and the depth for the same reason -- a term whose
    # penalty has a kink is fitted by a scheme of its own, and how that
    # scheme sweeps the term's own hyperparameters is a property of the
    # term. A criterion is asked of every hyperparameter of the model,
    # smooth ones included, and has no business carrying it.
    search = S7::new_property(S7::class_character,
                              default = quote(character(0)))
  )
)

#' @title S7 Class for Additive Terms
#' @name additive_term
#'
#' @description
#' The branch of [model_term()] whose terms contribute a block of design
#' columns \eqn{X_j \beta_j} to the linear predictor of one distribution
#' parameter. Everything a formula can write that is not a filter or a latent
#' mixture is on this branch: [linpar()], the five penalized terms, [random()],
#' [s()] and [te()], [nl()] and the break-point terms.
#'
#' A built term on this branch carries four things beyond the root's
#' properties: the block itself, the coefficient names, the blueprint that
#' reproduces the block on new rows, and the penalty on its coefficients.
#'
#' @details
#' # What a build fills, and what stays empty
#'
#' A specification has `X` empty and `blueprint` an empty list; [term_build()]
#' fills all four. The blueprint is the reason [term_predict()] can **reapply**
#' the mapping instead of rebuilding it: the factor levels, the contrasts, the
#' knots a basis was placed on and the spreads a standardization used are
#' recorded there at build time. [check_term()]'s subset check is what tests
#' that they really are.
#'
#' `penalty` is `NULL` for an unpenalized term. Reading it directly answers for
#' a term carrying one penalty over its whole block; [term_penalties()] is the
#' general form and covers a term whose penalty reaches part of its parameters
#' or which carries several.
#'
#' # The contract on this branch
#'
#' [term_matrix()] returns `X`, [term_coef_names()] returns `coef_names`,
#' [term_npar()] counts the columns and [term_penalty()] returns `penalty`, all
#' from the properties, so a subclass gets them without writing anything. What
#' a subclass owes is [term_build()] and [term_predict()].
#'
#' A term whose block moves with its own coefficients, as [nl()] and [seg()]
#' do, also implements [term_refresh()] and [term_value()], and answers
#' [term_jacobian_block()] to say whether the block it returns is a Jacobian or
#' a frozen working linearization.
#'
#' # The block need not be dense
#'
#' `X` is `class_any`, so a block may be a base matrix or any \pkg{Matrix}
#' class. A grouping indicator has one non-zero per row, and [random()] builds
#' it sparse; [linpar()] and the penalized terms build sparse when asked. Code
#' reading a block tests `is.matrix(x) && is.numeric(x)` or a two-dimensional
#' S4 object, since `is.matrix()` alone is `FALSE` for every \pkg{Matrix}
#' class.
#'
#' @inheritParams model_term
#' @param X The design block, one row per observation: a numeric matrix or a
#'   two-dimensional \pkg{Matrix}. Empty until [term_build()] fills it.
#' @param coef_names The block's coefficient names, one per column, prefixed by
#'   `label` when there is one. Filled by [term_build()].
#' @param blueprint A named list of everything needed to reproduce the mapping
#'   on new rows. Its contents are the subclass's business; nothing outside the
#'   term reads them. Empty until [term_build()] fills it.
#' @param penalty A \pkg{penalties7} penalty on the block's coefficients, or
#'   `NULL` when the term is unpenalized.
#'
#' @return Nothing: the class is abstract and cannot be instantiated.
#'   `additive_term()` throws. As a type it is the parent of every term that
#'   contributes design columns, carrying the four properties above beside
#'   [model_term()]'s six.
#'
#' @seealso [structural_term()] for the other branch, [term_build()],
#'   [term_predict()] and [term_matrix()] for the contract, and [linpar()] for
#'   the simplest term on this branch.
#'
#' @examples
#' d <- data.frame(x = rnorm(20), g = factor(rep(letters[1:4], 5)))
#'
#' # A specification carries no block; a build fills all four properties.
#' spec <- ridge(~ x + g)
#' c(X = length(spec@X), blueprint = length(spec@blueprint))
#' b <- term_build(spec, d)
#' dim(term_matrix(b))
#' term_coef_names(b)
#' b@penalty
#'
#' # An unpenalized term has a NULL penalty.
#' is.null(term_build(linpar(~ x), d)@penalty)
#'
#' # A grouping indicator is built sparse, so the block is not a base matrix.
#' class(term_matrix(term_build(random(~ 1 | g), d)))
#'
#' # Abstract: there is nothing to construct.
#' try(additive_term())
#'
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
#' The branch of [model_term()] whose terms rewrite the likelihood instead of
#' adding design columns. [gas()] drives one distribution parameter by a
#' score-driven recursion, so its contribution to the predictor is a state
#' rather than a product \eqn{X\beta}; [regime()] and the marginal break-point
#' terms replace the log-likelihood outright with one mixed over latent states.
#' A term on this branch has no block, no coefficients and no
#' [term_matrix()] method.
#'
#' @details
#' # Parameters, not coefficients
#'
#' What such a term estimates are its **own parameters**, named by
#' [term_params()] and each riding a chart named by [term_links()], so a
#' persistence stays inside \eqn{(-1, 1)} and a loading stays positive whatever
#' an optimizer proposes. [term_npar()] counts those parameters,
#' [term_start()] says where they start, and a [term_penalties()] entry indexes
#' into them, where an additive term's entry indexes into columns.
#'
#' The class adds no properties of its own: everything a built structural term
#' records goes in the root's `label` and in the subclass's own slots.
#'
#' # The two shapes
#'
#' A **filter** reports a predictor and its exact derivative: [term_filter()]
#' returns the level at each observation together with the Jacobian in the
#' term's parameters, propagated beside the state because the recursion is the
#' only place it can be computed. [term_adjoint()], [term_curvature()] and
#' [term_third()] carry the first three orders through the same recursion.
#' [gas()] has this shape.
#'
#' A **likelihood** term reports no predictor at all. [term_loglik()] returns
#' its own contribution, [term_posterior()] the smoothed latent states and
#' [term_hessian()] the observed information of the mixture. [regime()] and the
#' marginal break-point terms have this shape.
#'
#' A fitting layer tells them apart by which of the two generics answers, and
#' at most one structural term is allowed per model formula.
#'
#' @inheritParams model_term
#' @param blueprint A named list of everything needed to reproduce the term
#'   on new rows: the design's state, the levels a factor had, whatever the
#'   subclass's own filter or recursion needs. Its contents are that
#'   subclass's business and nothing outside the term reads them. Empty in a
#'   specification and filled by [term_build()], which is what
#'   [term_is_built()] reads on this branch, a structural term having no
#'   coefficient names to record being built in.
#'
#' @return Nothing: the class is abstract and cannot be instantiated.
#'   `structural_term()` throws. As a type it is the parent of [gas()],
#'   [regime()] and the marginal break-point terms, and carries
#'   [model_term()]'s six properties together with `blueprint`, the branch's
#'   record of having been built, empty in a specification.
#'
#' @seealso [additive_term()] for the other branch, [gas()] and [regime()] for
#'   the two shapes, [term_params()] and [term_links()] for what such a term
#'   estimates.
#'
#' @examples
#' # gas() and regime() are on this branch; every additive term is not.
#' vapply(list(gas(p = 1, q = 1), regime(k = 2), s(x, k = 5), ridge(~ x)),
#'        function(t) S7::S7_inherits(t, structural_term), logical(1))
#'
#' # A structural term names its own parameters instead of coefficients.
#' g <- gas(p = 1, q = 2)
#' term_params(g)
#' vapply(term_links(g), function(l) l@link_name, character(1))
#'
#' # And has no design block at all.
#' try(term_matrix(g))
#'
#' # Abstract: there is nothing to construct.
#' try(structural_term())
#'
#' @export
structural_term <- S7::new_class(
  name = "structural_term",
  parent = model_term,
  abstract = TRUE,
  properties = list(
    # The branch's own record of having been built, and the reason it is
    # declared here rather than on each subclass: every consumer reads it --
    # term_build() fills one, term_is_built() and the three print methods
    # test its length -- so a structural class written outside the package
    # inherits the contract instead of having to know about it. The default
    # is an empty list, which is what a specification carries.
    blueprint = S7::class_list
  )
)

#' @title S7 Class for the Unpenalized Parametric Term
#' @name LinparTerm
#'
#' @description
#' The subclass of [additive_term()] holding an unpenalized parametric block
#' built from a one-sided formula through [stats::model.matrix()]. It is what
#' [linpar()] constructs, and what [interpret_formula()] collects a formula's
#' bare covariates into, so most models carry one whether or not the caller
#' wrote it.
#'
#' Beyond the additive branch's four properties it records the formula, whether
#' the block is sparse, and the contrasts used for the formula's factors.
#'
#' @details
#' # The three properties of its own
#'
#' `formula` is the one-sided formula, kept with its environment, so a symbol
#' the data do not carry is still found where the formula was written.
#'
#' `sparse` holds what the caller asked for and `NULL` where nothing was asked.
#' The build leaves it alone and records the storage it settled on in
#' `blueprint$sparse`, which is the value [term_predict()] reads, so the block
#' and a prediction from it never differ in storage. A formula naming a factor
#' of many levels has one non-zero per row, and the dense model matrix of it is
#' the memory the choice exists to avoid.
#'
#' `contrasts` is a named list, one entry per factor, or empty for the
#' session's own `options("contrasts")`. Whatever is used is recorded in the
#' blueprint and reapplied, so a fit and a prediction never disagree about the
#' coding.
#'
#' The class carries no penalty: `penalty` is `NULL` on every `LinparTerm`, and
#' [term_penalties()] returns an empty list, so [edf()] counts its columns
#' exactly.
#'
#' @inheritParams additive_term
#' @param formula The one-sided formula defining the block, such as
#'   `~ x + log(z) + f`. Its environment is kept and used when a symbol is
#'   absent from the data.
#' @param sparse `TRUE` to build the block as a `dgCMatrix` through
#'   [Matrix::sparse.model.matrix()], `FALSE` for a base matrix, `NULL` to let
#'   the build choose. It is kept as given; the storage actually used is in
#'   `blueprint$sparse`.
#' @param contrasts A named list of contrasts for the formula's factors, in
#'   [stats::model.matrix()]'s own form, or an empty list for the session's
#'   defaults.
#'
#' @return An S7 object of class `LinparTerm`, inheriting from
#'   [additive_term()] and [model_term()], with the three properties above
#'   beside the ten they supply.
#'
#' @seealso [linpar()], the constructor to use; [interpret_formula()], which
#'   builds one implicitly; [term_build()] and [term_predict()].
#'
#' @examples
#' d <- data.frame(x = rnorm(20), g = factor(rep(letters[1:4], 5)))
#'
#' # linpar() is the constructor; the class is what it returns.
#' tm <- linpar(~ x + g)
#' S7::S7_inherits(tm, LinparTerm)
#' tm@formula
#'
#' # The property keeps what was asked for; the blueprint records what the
#' # build settled on, and that is what a prediction reads.
#' bt <- term_build(tm, d)
#' c(asked = is.null(bt@sparse), settled = bt@blueprint$sparse)
#' term_build(linpar(~ g, sparse = TRUE), d)@blueprint$sparse
#'
#' # The contrasts are recorded and reapplied, so the coding cannot drift.
#' b <- term_build(linpar(~ g, contrasts = list(g = "contr.sum")), d)
#' term_coef_names(b)
#'
#' # It is never penalized.
#' c(penalty = is.null(b@penalty), entries = length(term_penalties(b)),
#'   edf = edf(b))
#'
#' @export
LinparTerm <- S7::new_class(
  name = "LinparTerm",
  parent = additive_term,
  properties = list(
    formula = S7::class_any,
    # WHERE the block is built sparse rather than compressed afterwards: a
    # formula carrying a factor of many levels has one non-zero per row, and
    # the dense model matrix of it is the memory the choice exists to avoid.
    # NULL until the build settles it from the design; the settled value is
    # what the blueprint carries.
    sparse = S7::class_any,
    contrasts = S7::new_property(S7::class_list, default = quote(list()))
  )
)
