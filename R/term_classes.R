#' @title S7 Base Class for Model Terms
#' @name model_term
#'
#' @description
#' `model_term` is the abstract root of the term hierarchy. A term records what
#' a formula names: the recipe turning a data frame into a contribution to the
#' model, together with the metadata that a fit reads. A term written in a
#' formula is a specification that carries only what is written in the call;
#' [term_build()] turns it into a built term that carries the design block or
#' the state that the recipe produces on given data.
#'
#' [interpret_formula()] recognizes a term by this class. Any call on the right
#' of a formula whose value inherits from `model_term` becomes a term, so a
#' term class defined outside the package works in a formula as soon as it is
#' written, with nothing to register.
#'
#' @details
#' # The two branches
#'
#' [additive_term()] contributes a block of design columns \eqn{X_j\beta_j} to
#' a linear predictor. [structural_term()] rewrites the likelihood instead:
#' [gas()]'s predictor is a recursion, [regime()]'s contribution is a
#' likelihood mixed over latent states. The branch determines which generics a
#' consumer may call, so [statmodels7::statmod()] routes on it.
#'
#' # The properties of the root
#'
#' Beyond `label` the properties all describe a hyperparameter path. Only the
#' term determines which hyperparameters are estimated, because the penalty is
#' named in the term, so a held value, a grid size, a written-out set of
#' values, the depth of the path and the way in which several hyperparameters
#' are combined all belong to it. An outer criterion applies to every
#' hyperparameter of a model, including the smooth ones, and records none of
#' this.
#'
#' Each of `hyper`, `grid` and `values` is a named list keyed by the penalty's
#' own hyperparameter names, and each is empty by default. `min_ratio` and
#' `search` are single values rather than one per hyperparameter, because only
#' the path over the size of the kink uses a ratio: a bounded hyperparameter is
#' swept over its own interval and a shape over a geometric grid above its
#' lower bound.
#'
#' # Abstract classes
#'
#' `model_term` is abstract, as are `additive_term` and `structural_term`.
#' Calling `model_term()` signals the error
#' `"Can't construct an object from abstract class <model_term>"`. The class
#' serves as a parent when writing a term class, and as the test
#' `S7::S7_inherits(x, model_term)` when testing whether an object is a term.
#'
#' A class inheriting from it supplies [term_build()] at the very least: the
#' method registered here signals an error that names the class that did not
#' implement it. Most other generics have a usable default, so a new term class
#' starts from a working object and overrides the defaults that it has reason
#' to change.
#'
#' @param label A character string prefixed to the term's coefficient names
#'   when non-empty, and used as the title of [plot()] and the tag of
#'   [print()]. The empty string `""` means no label.
#' @param hyper The hyperparameters of the term's penalty that the caller
#'   held, as a named list keyed by the penalty's names. Empty, the
#'   default, means every one of them is estimated. See [term_hyper()].
#' @param grid How many values a path visits for each of the term's
#'   hyperparameters, as a named list of single whole numbers. Empty, the
#'   default, leaves the number to the fitting layer. See [term_grid()].
#' @param values The values that a path visits, for each hyperparameter that
#'   the caller wrote out, as a named list of numeric vectors. Empty, the
#'   default, has the path build them. See [term_values()].
#' @param min_ratio How far down the path over the size of the kink reaches, as
#'   a fraction of the value that empties the block: one number in \eqn{(0, 1)},
#'   or `numeric(0)` for the fitting layer's own. See [term_path_min()].
#' @param search How the term's own hyperparameters are covered when it has
#'   more than one carrying a kink: `"grid"` for every combination of them,
#'   `"cyclic"` for one at a time, or `character(0)` for the default. See
#'   [term_search()].
#' @param ids Which of the term's hyperparameters are shared with those of
#'   other terms, and under what label: a character vector named by the
#'   term's own hyperparameters, or `character(0)` for none. The terms
#'   carrying the same label for the same hyperparameter estimate one value.
#'   See [term_ids()].
#'
#' @return Nothing: the class is abstract and cannot be instantiated. As a type
#'   it is the parent of every term, with the properties above.
#'
#' @seealso [additive_term()] and [structural_term()] for the two branches,
#'   [term_build()] for turning a specification into a built term,
#'   [interpret_formula()] for how a formula is read into terms, and
#'   [check_term()] for validating one.
#'
#' @examples
#' # Every term inherits from this class, on either branch.
#' vapply(list(linpar(~ 1), ridge(~ x), s(x, basis7::bspline_smooth(k = 5)), gas(p = 1, q = 1)),
#'        function(t) S7::S7_inherits(t, model_term), logical(1))
#'
#' # A consumer routes on the branch of a term.
#' c(additive = S7::S7_inherits(s(x, basis7::bspline_smooth(k = 5)), additive_term),
#'   structural = S7::S7_inherits(s(x, basis7::bspline_smooth(k = 5)), structural_term))
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
                              default = quote(character(0))),
    # WHICH of its hyperparameters are shared with those of other terms, and
    # under what name: a character vector whose names are the term's own
    # hyperparameters and whose values are the labels. The terms carrying the
    # same label for the same hyperparameter estimate ONE value.
    #
    # It belongs here beside the other four for the same reason they do: what
    # a term's hyperparameters are and how they are estimated is the term's
    # own answer. Sharing does not merge the penalties -- they stay two
    # objects with the same number in them -- so this is a constraint on the
    # search and nothing else, which is what keeps every reader of a
    # hyperparameter working unchanged.
    ids = S7::new_property(S7::class_character,
                           default = quote(character(0)))
  )
)

#' @title S7 Class for Additive Terms
#' @name additive_term
#'
#' @description
#' `additive_term` is the branch of [model_term()] whose terms contribute a
#' block of design columns \eqn{X_j \beta_j} to the linear predictor of one
#' distribution parameter. Every term that a formula can write and that is not
#' a filter or a latent mixture is on this branch: [linpar()], the five
#' penalized terms, [random()], [s()] and [te()], [nl()] and the break-point
#' terms [seg()], [jump()] and [jseg()]. With `marginal = TRUE` the break-point
#' terms are structural terms instead.
#'
#' A built term on this branch carries four things beyond the root's
#' properties: the block itself, the coefficient names, the blueprint that
#' reproduces the block on new rows, and the penalty on its coefficients.
#'
#' @details
#' # The build
#'
#' A specification has `X` empty and `blueprint` an empty list; [term_build()]
#' fills `X`, `coef_names` and `blueprint`, and `penalty` as well for a term
#' that carries a penalty. The blueprint allows [term_predict()] to reapply
#' the mapping instead of rebuilding it: the factor levels, the
#' contrasts, the knots on which a basis was placed and the spreads that a
#' standardization used are recorded there at build time. [check_term()]'s
#' subset check tests that they are.
#'
#' `penalty` is `NULL` for an unpenalized term. Reading it directly gives the
#' penalty of a term that carries one penalty over its whole block;
#' [term_penalties()] is the general form and covers a term whose penalty
#' reaches part of its parameters or which carries several.
#'
#' # The methods of this branch
#'
#' [term_matrix()] returns `X`, [term_coef_names()] returns `coef_names`,
#' [term_npar()] counts the columns and [term_penalty()] returns `penalty`, all
#' from the properties, so a subclass inherits them without writing anything. A
#' subclass must supply [term_build()] and [term_predict()].
#'
#' A term whose block moves with its own coefficients, as [nl()] and [seg()]
#' do, also implements [term_refresh()] and [term_value()].
#' [term_jacobian_block()] reports whether the block is a Jacobian or a frozen
#' working linearization; it is `TRUE` by default, and a term whose block is
#' the latter overrides it.
#'
#' # The block need not be dense
#'
#' `X` is `class_any`, so a block may be a base matrix or any \pkg{Matrix}
#' class. A grouping indicator has one non-zero per row, and [random()] builds
#' it sparse; [linpar()] and the penalized terms build sparse where their
#' `sparse` argument is `TRUE`, or, when it is `NULL`, where the factor part of
#' the design is large. Code reading a block tests
#' `is.matrix(x) && is.numeric(x)` or a two-dimensional S4 object, since
#' `is.matrix()` alone is `FALSE` for every \pkg{Matrix} class.
#'
#' @inheritParams model_term
#' @param X The design block, one row per observation: a numeric matrix or a
#'   two-dimensional \pkg{Matrix}. It is empty until [term_build()] fills it.
#' @param coef_names The block's coefficient names, one per column, prefixed by
#'   `label` when there is one. [term_build()] fills it.
#' @param blueprint A named list of everything needed to reproduce the mapping
#'   on new rows. Its contents belong to the subclass, and only the term reads
#'   them. It is empty until [term_build()] fills it.
#' @param penalty A \pkg{penalties7} penalty on the block's coefficients, or
#'   `NULL` when the term is unpenalized.
#'
#' @return Nothing: the class is abstract and cannot be instantiated.
#'   Calling `additive_term()` signals an error. As a type it is the parent of
#'   every term that contributes design columns, carrying the four properties
#'   above beside those of [model_term()].
#'
#' @seealso [structural_term()] for the other branch, [term_build()],
#'   [term_predict()] and [term_matrix()] for the methods of the branch, and
#'   [linpar()] for the simplest term on this branch.
#'
#' @examples
#' d <- data.frame(x = rnorm(20), g = factor(rep(letters[1:4], 5)))
#'
#' # A specification carries no block; a build fills the block, the names and
#' # the blueprint, and the penalty where the term has one.
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
#' `structural_term` is the branch of [model_term()] whose terms rewrite the
#' likelihood instead of adding design columns. [gas()] drives one
#' distribution parameter by a score-driven recursion, so its contribution to
#' the predictor is a state rather than a product \eqn{X\beta}; [regime()] and
#' the marginal break-point terms replace the log-likelihood outright with one
#' mixed over latent states. A term on this branch has no block, no coefficients
#' and no [term_matrix()] method.
#'
#' @details
#' # The parameters of a structural term
#'
#' Such a term estimates its own parameters, which [term_params()] names and
#' each of which rides a chart named by [term_links()], so a persistence stays
#' inside \eqn{(-1, 1)} and a loading stays positive whatever an optimizer
#' proposes. [term_npar()] counts those parameters,
#' [term_start()] gives their starting values, and a [term_penalties()] entry
#' indexes into them, where an additive term's entry indexes into columns.
#'
#' The class adds one property of its own, `blueprint`, which records that the
#' term has been built. Everything else that a built structural term records is
#' held in the root's properties and in the subclass's own slots.
#'
#' # The two shapes
#'
#' A filter reports a predictor and its exact derivative: [term_filter()]
#' returns the level at each observation together with the Jacobian in the
#' term's parameters, propagated beside the state because the recursion is the
#' only place in which that Jacobian can be computed. [term_adjoint()],
#' [term_curvature()], [term_third()] and [term_fourth()] carry the derivatives
#' of the first to the fourth order through the same recursion. [gas()] has this
#' shape.
#'
#' A likelihood term reports no predictor at all. [term_loglik()] returns
#' its own contribution, [term_posterior()] the smoothed latent states and
#' [term_hessian()] the observed information of the mixture. [regime()] and the
#' marginal break-point terms have this shape.
#'
#' The two shapes are distinguished by which of the two generics
#' [term_filter()] and [term_loglik()] has a method for the term's class.
#' [statmodels7::statmod()] admits at most one structural term per model
#' formula.
#'
#' @inheritParams model_term
#' @param blueprint A named list of everything needed to reproduce the term
#'   on new rows: the design's state, the levels that a factor had, whatever
#'   the subclass's own filter or recursion needs. Its contents belong to that
#'   subclass and only the term reads them. It is empty in a specification and
#'   [term_build()] fills it. [term_is_built()] reads it on this branch,
#'   because a structural term has no coefficient names that would show
#'   whether it is built.
#'
#' @return Nothing: the class is abstract and cannot be instantiated.
#'   Calling `structural_term()` signals an error. As a type it is the parent
#'   of [gas()], [regime()] and the marginal break-point terms, and carries
#'   the properties of [model_term()] together with `blueprint`, the branch's
#'   record of having been built, empty in a specification.
#'
#' @seealso [additive_term()] for the other branch, [gas()] and [regime()] for
#'   the two shapes, [term_params()] and [term_links()] for what such a term
#'   estimates.
#'
#' @examples
#' # gas() and regime() are on this branch; the additive terms are on the other.
#' vapply(list(gas(p = 1, q = 1), regime(k = 2), s(x, basis7::bspline_smooth(k = 5)), ridge(~ x)),
#'        function(t) S7::S7_inherits(t, structural_term), logical(1))
#'
#' # A structural term names its own parameters instead of coefficients.
#' g <- gas(p = 1, q = 2)
#' term_params(g)
#' vapply(term_links(g), function(l) l@link_name, character(1))
#'
#' # A structural term has no design block at all.
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
#' `LinparTerm` is the subclass of [additive_term()] that holds an unpenalized
#' parametric block built from a one-sided formula through
#' [stats::model.matrix()]. [linpar()] constructs it, and [interpret_formula()]
#' collects a formula's bare covariates into it, so most models carry one
#' whether or not the caller wrote it.
#'
#' Beyond the additive branch's four properties it records the formula, whether
#' the block is sparse, the contrasts used for the formula's factors, and
#' whether the column `(Intercept)` is dropped after the coding.
#'
#' @details
#' # The properties of its own
#'
#' `formula` is the one-sided formula, kept with its environment, so a symbol
#' that the data do not carry is still found where the formula was written.
#'
#' `sparse` holds the value that the caller gave, and `NULL` where the caller
#' gave no value. The build leaves it unchanged and records in
#' `blueprint$sparse` the storage that it resolved, which is the value that
#' [term_predict()] reads, so the block and a prediction from it never differ
#' in storage. A formula naming a factor of many levels has one non-zero per
#' row, and building the dense model matrix of it is the memory use that
#' sparse storage avoids.
#'
#' `contrasts` is a named list, one entry per factor, or empty for the
#' session's own `options("contrasts")`. Whatever is used is recorded in the
#' blueprint and reapplied, so a fit and a prediction always use the same
#' coding.
#'
#' `drop_intercept` is `FALSE` by default. When it is `TRUE`, the build removes
#' the column `(Intercept)` after the coding, and the blueprint records the
#' setting so that a prediction applies it again.
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
#' @param drop_intercept `TRUE` to remove the column `(Intercept)` after the
#'   model matrix is built, so that the factors keep the coding that they have
#'   in a formula with an intercept. [interpret_formula()] sets it on a
#'   `linpar()` written in a formula that already carries an intercept. The
#'   default is `FALSE`.
#'
#' @return An S7 object of class `LinparTerm`, inheriting from
#'   [additive_term()] and [model_term()], with the four properties above
#'   beside the properties that it inherits.
#'
#' @seealso [linpar()], the constructor to use; [interpret_formula()], which
#'   builds one implicitly; [term_build()] and [term_predict()].
#'
#' @examples
#' d <- data.frame(x = rnorm(20), g = factor(rep(letters[1:4], 5)))
#'
#' # linpar() is the constructor and returns an object of this class.
#' tm <- linpar(~ x + g)
#' S7::S7_inherits(tm, LinparTerm)
#' tm@formula
#'
#' # The property keeps the value given; the blueprint records the storage
#' # that the build used, and a prediction reads the blueprint.
#' bt <- term_build(tm, d)
#' c(sparse_is_null = is.null(bt@sparse), settled = bt@blueprint$sparse)
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
    contrasts = S7::new_property(S7::class_list, default = quote(list())),
    # the formula's intercept lives in the implicit block, so a linpar()
    # written beside it drops its own intercept column AFTER the coding: a
    # formula rewritten as ~ . - 1 would code its factors with full
    # indicators and alias one of them against that intercept
    drop_intercept = S7::new_property(S7::class_logical, default = FALSE)
  )
)
