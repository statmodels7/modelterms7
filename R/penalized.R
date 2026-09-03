#' @include term_classes.R generics.R
NULL

#' @title S7 Class for Penalized Parametric Terms
#' @name PenalizedTerm
#'
#' @description
#' The subclass of [additive_term()] for a parametric block whose coefficients
#' carry a \pkg{penalties7} penalty. [ridge()], [lasso()], [enet()], [scad()]
#' and [mcp()] all construct it; the five differ only in the penalty their
#' factory attaches at build time, so every derivative, hyperparameter, bound,
#' link and kink belongs to the penalty object and none of it is restated here.
#'
#' @details
#' # The five properties of its own
#'
#' `input` is the block as it was given, a one-sided formula or a matrix.
#' `input_expr` is the expression that produced a matrix input, kept so that
#' [term_predict()] can re-evaluate it in new data.
#'
#' `factory` maps a coefficient count to the penalty. It is called at build,
#' when the count is finally known, and with the diagonal map as a second
#' argument when `standardize` asks for one, so a factory that will never be
#' standardized may take the count alone.
#'
#' `standardize` is a single logical on the specification; after a build the
#' spreads themselves are in `blueprint$standardize`, one per column.
#'
#' `sparse` governs the formula route only, a matrix input being kept in
#' whatever storage it arrives in. `NULL` until the build settles it.
#'
#' # What the class does not decide
#'
#' Which hyperparameters are estimated, over what grid and how far down a path
#' are recorded in [model_term()]'s own properties, read through
#' [term_hyper()], [term_grid()], [term_values()], [term_path_min()] and
#' [term_search()]. [term_smooth()] is `TRUE` for a built [ridge()] and `FALSE`
#' for the other four, read from each penalty's kink set.
#'
#' @inheritParams additive_term
#' @param input The block as given: a one-sided formula or a numeric matrix.
#' @param input_expr The expression that produced a matrix input, kept so that
#'   [term_predict()] can re-evaluate it in new data. `NULL` on the formula
#'   route.
#' @param factory A function of the coefficient count, returning the
#'   \pkg{penalties7} penalty over that many coefficients, and taking the
#'   diagonal map as a second argument where `standardize` asks for one.
#' @param sparse `TRUE`, `FALSE` or `NULL` for the formula route's storage;
#'   see [linpar()] for the rule `NULL` is settled by. A matrix input needs no
#'   such argument.
#' @param standardize A single logical: whether the block's columns are put on
#'   a common scale by the penalty's diagonal map. `FALSE` by default.
#'
#' @return An S7 object of class `PenalizedTerm`, inheriting from
#'   [additive_term()] and [model_term()], with the five properties above
#'   beside the ten they supply.
#'
#' @seealso [penalized_terms()] for what the five constructors share,
#'   [ridge()] and its four siblings for their own hyperparameters,
#'   [print.PenalizedTerm()] for how one displays.
#'
#' @examples
#' set.seed(3)
#' dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
#'
#' # Every one of the five returns this class.
#' vapply(list(ridge(~ x1), lasso(~ x1), enet(~ x1), scad(~ x1), mcp(~ x1)),
#'        function(t) S7::S7_inherits(t, PenalizedTerm), logical(1))
#'
#' # The penalty is attached at build, when the column count is known.
#' spec <- lasso(~ x1 + x2)
#' c(spec = is.null(spec@penalty))
#' b <- term_build(spec, dd)
#' b@penalty
#' b@penalty@n_coef
#'
#' # standardize is a flag before the build and the spreads after it.
#' dd$x3 <- 1000 * dd$x2
#' bs <- term_build(lasso(~ x1 + x3, standardize = TRUE), dd)
#' bs@standardize
#' bs@blueprint$standardize
#'
#' @export
PenalizedTerm <- S7::new_class(
  name = "PenalizedTerm",
  parent = additive_term,
  properties = list(
    input = S7::class_any,
    input_expr = S7::class_any,
    factory = S7::class_function,
    standardize = S7::new_property(S7::class_logical, default = FALSE),
    # the FORMULA route's storage. A matrix input is kept in whatever
    # storage it arrives in and needs no argument to say so. NULL until the
    # build settles it from the design.
    sparse = S7::class_any
  )
)

.penalized_spec <- function(x, expr, label, standardize, factory,
                            hyper = list(), extra = list(), grid = list(),
                            min_ratio = NULL, search = NULL, sparse = NULL,
                            ids = NULL) {
  # An argument named after ANOTHER penalty's hyperparameter is the mistake
  # this catches: `mcp(x, a = 3)` writes SCAD's shape on an MCP, whose own
  # is gamma, and reaches nothing. R would report it as an unused argument,
  # which says that it was not read and not what to write instead.
  #
  # It is also where an argument the term does not have arrives -- `by` on a
  # separable penalty, a grid size on a ridge -- so the reply names what was
  # written and then the hyperparameters, which is what the caller was most
  # likely reaching for, rather than claiming to list every argument.
  if (length(extra)) {
    pen <- tryCatch(factory(1L), error = function(e) NULL)
    stop(sprintf(paste0("'%s' has no argument '%s'.\n  Its hyperparameters",
                        " are: %s, each held by naming it and estimated when",
                        " left NULL."),
                 label, names(extra)[[1L]],
                 if (is.null(pen)) "none" else paste(pen@params,
                                                     collapse = ", ")),
         call. = FALSE)
  }
  if (!is.logical(standardize) || length(standardize) != 1L ||
      is.na(standardize)) {
    stop("'standardize' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !nzchar(label)) {
    stop("'label' must be a single non-empty character string.",
         call. = FALSE)
  }
  is_formula <- inherits(x, "formula")
  if (is_formula) {
    if (length(x) != 2L) {
      stop("a formula input must be one-sided, e.g. ~ x1 + x2.",
           call. = FALSE)
    }
  } else if (isS4(x) && methods::is(x, "Matrix")) {
    # A Matrix is KEPT as it is. Densifying a caller's sparse design here
    # would undo exactly the saving it was passed for, and the factor is
    # 1/density: measured on a 4000 x 60 indicator design at density 0.017,
    # 0.050 MB becomes 1.920 MB. A logical Matrix is carried to double
    # rather than rejected, an indicator being the commonest sparse input.
    if (!methods::is(x, "dMatrix")) {
      x <- tryCatch(methods::as(x, "dMatrix"), error = function(e) {
        stop("a Matrix input must be numeric or logical.", call. = FALSE)
      })
    }
  } else {
    x <- as.matrix(x)
    if (!is.numeric(x)) {
      stop("a matrix input must be numeric.", call. = FALSE)
    }
  }
  vals <- check_values(hyper, factory, label)
  reject_pathless_values(vals, tryCatch(factory(1L), error = function(e) NULL),
                         label)
  # a MATRIX input is kept in whatever storage it arrives in, so `sparse`
  # governs the FORMULA route, where the model matrix would otherwise be
  # built dense whatever the caller passed
  PenalizedTerm(label = label, input = x, input_expr = expr,
                factory = factory, standardize = standardize,
                sparse = .check_design_opts(sparse, NULL, label)$sparse,
                hyper = check_hyper(hyper, factory, label),
                values = vals,
                grid = check_grid(grid, factory, label),
                min_ratio = check_min_ratio(min_ratio, label),
                search = check_search(search, label),
                ids = check_ids(ids, .factory_params(factory), label),
                X = NULL, coef_names = character(0),
                blueprint = list(), penalty = NULL)
}

# The spread of each column, read without densifying: for a sparse block
# the two-pass form would have to subtract the mean from every entry, so
# the sum of squares is used there and the two-pass elsewhere, where the
# cancellation it avoids is real for a column with a large mean.
#' @title Column Standard Deviations of a Block, in Either Storage
#'
#' @description
#' The standard deviation of each column of a design block, computed without
#' densifying a \pkg{Matrix}. It is the spread `standardize = TRUE` divides
#' each coefficient by.
#'
#' @details
#' On a base matrix it is [stats::sd()] column by column. On a \pkg{Matrix} it
#' is assembled from the column sums and sums of squares,
#' \eqn{s_j^2 = (\sum_i x_{ij}^2 - n\bar{x}_j^2)/(n-1)}, so a sparse block is
#' never expanded. The sum of squares is floored at zero before the square
#' root, that expression being able to go slightly negative by rounding on a
#' nearly constant column.
#'
#' A column with no spread takes \eqn{s_j = 1}, which penalizes its
#' coefficient on its own scale instead of dividing by zero. That covers a
#' constant column, a single-row block and any column whose spread comes back
#' non-finite.
#'
#' @param X A design block: a numeric matrix or a two-dimensional \pkg{Matrix}.
#'
#' @return A numeric vector of length `ncol(X)`, every entry finite and
#'   strictly positive.
#'
#' @seealso [penalized_terms()] for the standardization this feeds,
#'   [penalties7::map_diagonal()] for the map it becomes.
#'
#' @keywords internal
.block_sd <- function(X) {
  if (isS4(X) && methods::is(X, "Matrix")) {
    n <- nrow(X)
    m <- Matrix::colMeans(X)
    s <- sqrt(pmax(Matrix::colSums(X^2) - n * m^2, 0) / max(n - 1L, 1L))
  } else {
    s <- apply(X, 2L, stats::sd)
  }
  s <- as.numeric(s)
  # a constant column has nothing to standardize by, and its coefficient is
  # penalized on its own scale rather than divided by zero
  s[!is.finite(s) | s <= 0] <- 1
  s
}

#' What the Penalized Terms Share
#'
#' @description
#' The input handling, the storage, the standardization and the prediction of
#' the five penalized terms, documented once. [ridge()], [lasso()], [enet()],
#' [scad()] and [mcp()] each have a page of their own carrying the penalty's
#' formula, its hyperparameters and where those may lie.
#'
#' Each takes its block as a one-sided formula or as a numeric matrix, and
#' attaches the corresponding \pkg{penalties7} object to the block's
#' coefficients at build time. The hyperparameters, their bounds and links, the
#' derivatives and the kink set are the penalty's and are never restated by the
#' term.
#'
#' @details
#' # The two inputs
#'
#' A **formula** goes through the [stats::model.matrix()] machinery with the
#' intercept removed: a penalized block does not penalize an intercept, and the
#' model's intercept lives in the parametric block. Its blueprint records the
#' terms, the factor levels and the contrasts, exactly as [linpar()] does, so
#' `~ g` and `~ 0 + g` give the same four columns for a four-level factor.
#'
#' The exception is a formula whose intercept is all it has. `ridge(~ 1)` is a
#' block of that one column under the penalty, removing it leaving no block at
#' all. That is the form a subformula on another term's parameter uses:
#' `gamma ~ 0 + lasso(~ 1)` says the parameter itself carries a lasso, there
#' being no parametric block of its own to hold an unpenalized intercept.
#'
#' A **matrix** is used as given, and its columns are named after its own
#' column names, or numbered `1`, `2`, ... when it has none.
#'
#' # Predicting a matrix input
#'
#' Prediction re-evaluates the expression that produced the matrix **in the new
#' data, and only there**. So the intended use is a matrix column of the model
#' data frame:
#'
#' ```r
#' dd$R <- R
#' ridge(R)      # in the formula
#' ```
#'
#' A subset of `dd` then carries the matching rows of `R`. A free-standing
#' matrix from the calling environment builds, its value having been captured,
#' and prediction is refused with a message naming the expression: resolving it
#' outside the new data would silently reuse the build-time rows whenever the
#' counts happened to agree.
#'
#' @section Sparse storage:
#' A **matrix** input is kept in whatever storage it arrives in, so a
#' `dgCMatrix` passed to any of the five stays one and nothing needs to be
#' said.
#'
#' `sparse` governs the **formula** route, where the model matrix would
#' otherwise be built dense whatever the columns look like. `TRUE` goes through
#' [Matrix::sparse.model.matrix()], which builds the block sparse; building a
#' dense one and compressing it would cost the memory the choice exists to
#' avoid.
#'
#' It pays where the formula carries a factor of many levels, whose indicator
#' columns hold one non-zero per row, and `lasso(~ 0 + g)` over hundreds of
#' groups is the case. On numeric covariates the block is dense whatever is
#' asked for, and the sparse storage then costs more than it saves. Left
#' `NULL`, the default, the storage is settled at build by [.resolve_sparse()]:
#' the dense indicator part holds `n` times its column count in cells against
#' one non-zero per row, and the two routes cross at about \eqn{10^5} of those
#' cells.
#'
#' Standardization does not interfere. It is a diagonal map on the **penalty**
#' and never an operation on the design, so a sparse block stays sparse under
#' it.
#'
#' @section Standardization:
#' A hyperparameter is comparable across coordinates only where the coordinates
#' share a scale. Without `standardize` a lasso penalizes a column measured in
#' meters more than the same column measured in kilometers, and a reader of
#' \eqn{\lambda} has no way to tell.
#'
#' `standardize = TRUE` divides each coefficient by the standard deviation of
#' its own column, through the penalty's diagonal map. With \eqn{z_j = x_j/s_j}
#' the coefficients satisfy \eqn{\beta_{z,j} = s_j\beta_{x,j}}, so
#'
#' \deqn{\lambda\sum_j \lvert\beta_{z,j}\rvert
#'   = \lambda\sum_j s_j\lvert\beta_{x,j}\rvert
#'   = \rho(S\beta_x), \qquad S = \mathrm{diag}(s),}
#'
#' which is the standardized penalty read on the original scale. Three things
#' follow. The design is never rescaled, so a sparse block stays sparse.
#' \eqn{\lambda} stays one number, and the coefficients are already on the
#' scale the data came in, with nothing to map back. And centering, which is
#' what would destroy sparsity, is not needed: the fit is invariant to a
#' translation of a penalized column wherever an intercept is free.
#'
#' The spread is computed from the built block by [.block_sd()] and frozen in
#' `blueprint$standardize`, so the same term standardizes identically in every
#' equation of a distributional model and does not move with the working
#' weights of a fit. A constant column takes \eqn{s_j = 1}. [print()] shows the
#' values, a number that changes the meaning of \eqn{\lambda} having to be
#' legible.
#'
#' For SCAD and MCP the map is not a rescaling of \eqn{\lambda} alone.
#' Substituting \eqn{s_j\beta_j} gives \eqn{\lambda_j = \lambda s_j} **and**
#' \eqn{a_j = a/s_j} (or \eqn{\gamma_j = \gamma/s_j}), a composition of both
#' hyperparameters per coordinate, which the map expresses exactly.
#'
#' [random()] does not standardize and takes no such argument. Its columns are
#' grouping indicators and its penalty is a variance component with a meaning
#' of its own; weighting it by the size of the groups would change the model.
#'
#' @param x A one-sided formula, such as `~ x1 + x2` or `~ 0 + g`, or a numeric
#'   matrix, ideally a matrix column of the model data frame. A two-sided
#'   formula throws.
#' @param label A single non-empty character string prefixed to the coefficient
#'   names as `label.name`. Each constructor defaults it to its own name, so a
#'   ridge over `x1` reads `ridge.x1`, and two penalized terms in one formula
#'   stay apart by their labels.
#' @param standardize A single logical, `FALSE` by default: whether to penalize
#'   each coefficient on the scale of its own column. See the section above.
#' @param sparse Governs the formula route: `TRUE` builds a `dgCMatrix` through
#'   [Matrix::sparse.model.matrix()], `FALSE` a dense model matrix, and `NULL`,
#'   the default, settles it at build from the size of the design. A matrix
#'   input needs no such argument. See the section above.
#' @param ... Not used, and reported. An argument named after another penalty's
#'   hyperparameter is the mistake this catches: `scad(~ x, gamma = 1)` throws
#'   `"'scad' has no argument 'gamma'."` and lists the ones it does have.
#'
#' @return An unbuilt [PenalizedTerm()]: a specification, with `X`,
#'   `coef_names`, `blueprint` and `penalty` all empty until [term_build()]
#'   fills them.
#'
#' @seealso [ridge()], [lasso()], [enet()], [scad()] and [mcp()] for the five
#'   penalties; [linpar()] for the unpenalized block; [random()] and [s()] for
#'   the penalized structures; [term_penalty()] and [edf()] for what a built
#'   one reports.
#'
#' @examples
#' set.seed(3)
#' dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20),
#'                  g = factor(rep(letters[1:4], 5)))
#'
#' # The intercept is removed, so a four-level factor gives four columns.
#' built <- term_build(lasso(~ x1 + x2), dd)
#' term_coef_names(built)
#' term_coef_names(term_build(lasso(~ g), dd))
#'
#' # Unless the intercept is all the formula has.
#' term_coef_names(term_build(ridge(~ 1), dd))
#'
#' # The hyperparameters and the kink are the penalty's.
#' term_penalty(built)@params
#' term_smooth(built)
#'
#' # A matrix column of the data is the input that predicts. The name must
#' # resolve where the term is written as well as in the data.
#' R <- matrix(rnorm(60), 20, 3, dimnames = list(NULL, c("a", "b", "c")))
#' dd$R <- R
#' bm <- term_build(ridge(R), dd)
#' term_coef_names(bm)
#' dim(term_predict(bm, dd[1:5, ]))
#'
#' # Standardizing puts the spreads on the penalty's map, not on the design.
#' dd$x3 <- 1000 * dd$x2
#' bs <- term_build(lasso(~ x1 + x3, standardize = TRUE), dd)
#' term_penalty(bs)@map
#' apply(term_matrix(bs), 2, sd)
#'
#' # An argument named after another penalty's hyperparameter is reported.
#' try(scad(~ x1, gamma = 1))
#'
#' @name penalized_terms
#' @keywords internal
NULL

#' @title Build a Penalized Block and Attach Its Penalty
#' @name term_build.PenalizedTerm
#'
#' @description
#' Builds the block of a [ridge()], [lasso()], [enet()], [scad()] or [mcp()]
#' term, names its coefficients, records the blueprint, and calls the term's
#' factory to produce the \pkg{penalties7} penalty over exactly that many
#' coefficients. The penalty exists only after this: a specification's
#' `penalty` property is `NULL`, its width being unknown until the data are
#' seen.
#'
#' @details
#' # Two input routes, two blueprints
#'
#' A **formula** is built through [stats::model.frame()] and
#' [.design_matrix()] with the intercept removed by
#' `update(formula, ~ . - 1)`, unless the intercept is all the formula has, in
#' which case it is kept and is the block. The blueprint records
#' `kind = "formula"`, the terms object, the levels, the contrasts and the
#' settled storage.
#'
#' A **matrix** is used as it stands, and its row count must equal
#' `nrow(data)`; anything else throws with both numbers. Its columns take the
#' matrix's own names, or `1`, `2`, ... where it has none. The blueprint
#' records `kind = "matrix"`, the expression that produced it and those base
#' names.
#'
#' Either way every column name is prefixed with the term's label, and the
#' block's row names are dropped.
#'
#' # Standardization is applied to the penalty
#'
#' With `standardize = TRUE` the column spreads come from [.block_sd()], go
#' into `blueprint$standardize` and become a [Matrix::Diagonal()] map. The
#' design is not touched, so a sparse block stays sparse.
#'
#' The map is passed to the factory **as a second argument, and only when there
#' is one**. A factory written before standardization existed takes the count
#' alone and goes on working, provided it is never standardized.
#'
#' @param term An unbuilt or built [PenalizedTerm()].
#' @param data A data frame with as many rows as a matrix input has, carrying
#'   every variable a formula input names.
#' @param ... Unused.
#'
#' @return The term with `X`, `coef_names`, `blueprint` and `penalty` filled.
#'   `penalty@n_coef` equals `ncol(X)`.
#'
#' @seealso [penalized_terms()] for what the five constructors share,
#'   [term_predict.PenalizedTerm()] for the block at new rows,
#'   [.block_sd()] for the spreads.
#'
#' @examples
#' set.seed(3)
#' dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
#'
#' # The penalty appears at build, over exactly the columns built.
#' b <- term_build(lasso(~ x1 + x2), dd)
#' c(cols = ncol(term_matrix(b)), n_coef = b@penalty@n_coef)
#'
#' # A matrix input must have as many rows as the data.
#' M <- matrix(rnorm(30), 10, 3)
#' try(term_build(ridge(M), dd))
#'
#' @keywords internal
S7::method(term_build, PenalizedTerm) <- function(term, data, ...) {
  if (inherits(term@input, "formula")) {
    # The intercept is dropped because the model's intercept lives in the
    # parametric block and is not penalized -- unless it is all the formula
    # has, where dropping it leaves no block at all. `ridge(~1)` then means
    # the single coefficient under the penalty, which is what a subformula
    # on a term's own parameter needs: `gamma ~ lasso(~1)` says the change
    # itself carries a lasso, there being no parametric block of its own to
    # hold an unpenalized intercept.
    f <- if (.intercept_only(term@input)) term@input else {
      stats::update(term@input, ~ . - 1)
    }
    environment(f) <- environment(term@input)
    mf <- stats::model.frame(f, data,
                             na.action = stats::na.pass,
                             drop.unused.levels = FALSE)
    tt <- attr(mf, "terms")
    sp <- .resolve_sparse(term@sparse, nrow(mf), .indicator_cols(tt, mf))
    b <- .design_matrix(tt, mf, NULL, sp)
    X <- b$X
    base_names <- colnames(X)
    term@blueprint <- list(
      kind = "formula",
      terms = stats::delete.response(tt),
      xlev = stats::.getXlevels(tt, mf),
      contrasts = b$contrasts,
      # the SETTLED storage, so a prediction reproduces the build's kind
      sparse = sp
    )
  } else {
    X <- term@input
    if (nrow(X) != nrow(data)) {
      stop(sprintf("the matrix input has %d rows and 'data' has %d.",
                   nrow(X), nrow(data)), call. = FALSE)
    }
    base_names <- colnames(X)
    if (is.null(base_names)) base_names <- as.character(seq_len(ncol(X)))
    term@blueprint <- list(
      kind = "matrix",
      expr = term@input_expr,
      base_names = base_names
    )
  }
  cn <- paste(term@label, base_names, sep = ".")
  colnames(X) <- cn
  rownames(X) <- NULL
  term@X <- X
  term@coef_names <- cn
  # standardization is a weight per coordinate on the penalty and never an
  # operation on X: with z_j = x_j/s_j the coefficient is beta_z,j = s_j
  # beta_x,j, so rho(s beta) is the standardized penalty read on the
  # original scale. The design keeps its scale and its sparsity, lambda
  # stays one number, and the coefficients need no mapping back.
  map <- NULL
  if (isTRUE(term@standardize)) {
    s <- stats::setNames(.block_sd(X), cn)
    term@blueprint$standardize <- s
    map <- Matrix::Diagonal(x = as.numeric(s))
  }
  # the factory's second argument is the map, and it is passed only when
  # there is one: a subclass written before standardization existed carries a
  # one-argument factory, and it is never standardized, so the old contract
  # goes on working unchanged
  term@penalty <- if (is.null(map)) term@factory(ncol(X)) else
    term@factory(ncol(X), map)
  term
}

# Whether a formula carries an intercept and nothing else, which is the one
# case where removing the intercept would leave no block. `~1` and `~+1`
# reach here; `~0` and `~-1` do not, having neither an intercept nor a term
# and being an empty block whichever rule is applied.
.intercept_only <- function(f) {
  tt <- stats::terms(f)
  attr(tt, "intercept") == 1L && length(attr(tt, "term.labels")) == 0L
}

#' @title A Penalized Block at New Rows
#' @name term_predict.PenalizedTerm
#'
#' @description
#' Reproduces a built penalized block at `newdata`. A formula input is
#' reapplied through the recorded terms, levels, contrasts and storage, exactly
#' as [linpar()]'s is. A matrix input has its expression re-evaluated **in
#' `newdata` alone**, and the result is checked to have the right number of
#' rows and columns before it is returned.
#'
#' @details
#' The matrix branch evaluates `blueprint$expr` with `newdata` as the data and
#' [baseenv()] as the enclosure, so the calling environment is not on the
#' search path. A free-standing matrix therefore cannot be found, and the error
#' names the expression and says to supply it as a column of the data. That is
#' deliberate: resolving the name outside `newdata` would silently return the
#' build-time rows whenever the row counts happened to agree, which is a wrong
#' answer rather than an error.
#'
#' Two shape checks follow the evaluation, because an expression that resolves
#' in `newdata` may still give the wrong block: a row count differing from
#' `nrow(newdata)` and a column count differing from the build's each throw
#' with both numbers.
#'
#' @param term A built [PenalizedTerm()]. An unbuilt one throws
#'   `"the term has not been built; call term_build(term, data) first."`.
#' @param newdata A data frame. For a matrix input it must carry the matrix as
#'   a column of the same name the term was written with.
#' @param ... Unused.
#'
#' @return A block of `nrow(newdata)` rows and [term_npar()] columns, in the
#'   storage the build used, with the term's coefficient names as column names
#'   and no row names.
#'
#' @seealso [term_predict()] for the generic and the identity it satisfies,
#'   [penalized_terms()] for why a matrix input belongs in the data frame.
#'
#' @examples
#' set.seed(3)
#' dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
#'
#' # A formula input reapplies its blueprint.
#' b <- term_build(lasso(~ x1 + x2), dd)
#' all.equal(term_predict(b, dd), term_matrix(b))
#' dim(term_predict(b, dd[1:5, ]))
#'
#' # A matrix column of the data predicts on a subset of those rows.
#' R <- matrix(rnorm(60), 20, 3)
#' dd$R <- R
#' bm <- term_build(ridge(R), dd)
#' dim(term_predict(bm, dd[1:5, ]))
#'
#' # A free-standing matrix builds and cannot be predicted from.
#' Rfree <- matrix(rnorm(60), 20, 3)
#' bf <- term_build(ridge(Rfree), dd)
#' try(term_predict(bf, dd[1:5, ]))
#'
#' @keywords internal
S7::method(term_predict, PenalizedTerm) <- function(term, newdata, ...) {
  .assert_built(term)
  bp <- term@blueprint
  if (identical(bp$kind, "formula")) {
    mf <- stats::model.frame(bp$terms, newdata,
                             na.action = stats::na.pass,
                             xlev = bp$xlev)
    X <- .design_matrix(bp$terms, mf, bp$contrasts, isTRUE(bp$sparse))$X
  } else {
    # kept sparse for the reason the constructor keeps it: a prediction that
    # densified would spend at new data what the build was careful not to
    X <- tryCatch(.as_block(eval(bp$expr, newdata, baseenv())),
                  error = function(e) {
                    stop(sprintf(paste("the matrix expression `%s` could not be",
                                       "evaluated in 'newdata' (%s); supply the",
                                       "matrix as a column of the data."),
                                 deparse(bp$expr), conditionMessage(e)),
                         call. = FALSE)
                  })
    if (nrow(X) != nrow(newdata)) {
      stop(sprintf("the matrix expression yields %d rows against %d in 'newdata'.",
                   nrow(X), nrow(newdata)), call. = FALSE)
    }
    if (ncol(X) != length(term@coef_names)) {
      stop(sprintf("the matrix expression yields %d columns against %d at build time.",
                   ncol(X), length(term@coef_names)), call. = FALSE)
    }
  }
  colnames(X) <- term@coef_names
  rownames(X) <- NULL
  X
}
