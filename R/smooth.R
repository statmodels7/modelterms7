#' @include term_classes.R generics.R
NULL

#' @title S7 Class for Smooth Terms
#' @name SmoothTerm
#'
#' @description
#' The subclass of [additive_term()] for a penalized smooth: a \pkg{basis7}
#' expansion of one or more covariates with a roughness penalty on its
#' coefficients. [s()] constructs it for one covariate and [te()] for several,
#' and the two differ in the basis they build and the penalty they attach, not
#' in the class.
#'
#' @details
#' # The four properties of its own
#'
#' `vars` holds the covariate expressions, one for [s()] and two or more for
#' [te()], kept unevaluated so that a build reads them in whatever data it is
#' given. `by` is the expression the smooth varies with, or `NULL`.
#'
#' `spec` carries the construction settings the build reads: the basis or
#' bases, the dimension `k`, the degree, whether the linear part is separated
#' out, and for [te()] whether the penalty is anisotropic. What a build then
#' computes from the data goes into the blueprint instead: the
#' Demmler-Reinsch transform, the centering constraint, the `by` levels.
#'
#' `sparse` is `NULL` until the build settles it. A smooth's block is sparse
#' only under a **factor** `by`, where each row sits in the block of its own
#' level: the basis itself is dense by construction, the Demmler-Reinsch
#' rotation making it so, and a numeric `by` merely multiplies it.
#'
#' @inheritParams additive_term
#' @param vars A list of the covariate expressions being smoothed.
#' @param by The expression the smooth varies with, or `NULL`.
#' @param spec A named list of construction settings: the basis, its dimension
#'   and degree, whether the linear part is carried separately, and for [te()]
#'   the `anisotropic` flag.
#' @param sparse `TRUE`, `FALSE` or `NULL` for the block's storage; only a
#'   factor `by` admits `TRUE`. See [s()].
#'
#' @return An S7 object of class `SmoothTerm`, inheriting from
#'   [additive_term()] and [model_term()], with the four properties above
#'   beside the ten they supply.
#'
#' @seealso [s()] and [te()], the two constructors; [term_penalty()] for the
#'   roughness penalty; [edf()] for what a fitted smooth spends.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(80)), z = runif(80))
#'
#' # Both constructors return this class.
#' c(s = S7::S7_inherits(s(x), SmoothTerm),
#'   te = S7::S7_inherits(te(x, z), SmoothTerm))
#'
#' # The settings are on `spec`; what the data decide is in the blueprint.
#' tm <- s(x, k = 8)
#' names(tm@spec)
#' names(term_build(tm, dd)@blueprint)
#'
#' @export
SmoothTerm <- S7::new_class(
  name = "SmoothTerm",
  parent = additive_term,
  properties = list(
    vars = S7::class_list,
    by = S7::class_any,
    spec = S7::class_list,
    # WHERE the block is sparse: a factor `by` puts each row in the block of
    # its own level, a density of 1/m. The basis itself is dense by
    # construction -- the Demmler-Reinsch rotation makes it so -- and so is
    # a numeric `by`, which merely multiplies it. NULL until the build
    # settles it; the settled value is what the blueprint carries.
    sparse = S7::class_any
  )
)

#' Penalized Smooth of One Covariate
#'
#' @description
#' A smooth function of one covariate, expanded in a \pkg{basis7} basis and
#' penalized for roughness. The default is a cubic B-spline basis under the
#' Demmler-Reinsch reparametrization, which separates the linear effect from
#' the nonlinear deviation and turns the roughness penalty into the identity on
#' the deviation.
#'
#' As the smoothing parameter grows the fit approaches a straight line, and
#' [edf()] falls to exactly one.
#'
#' @details
#' # The block and its penalty
#'
#' The block has one column for the linear effect, centered and scaled,
#' followed by the reparametrized basis, so `s(x, k = 8)` gives seven columns
#' named `s(x).lin`, `s(x).z1` ... `s(x).z6`. That ordering is what the penalty
#' reads: it is the quadratic penalty of \eqn{\mathrm{diag}(0, 1, \dots, 1)},
#' rank deficient by exactly one, so the linear effect is unpenalized and the
#' deviation is shrunk toward zero.
#'
#' Two consequences a reader of a fit needs. At a large smoothing parameter the
#' fit tends to a straight line, so `edf()` runs from `k - 1` down to 1 and
#' never to 0. And the linear column is orthogonal to the rest over the
#' observed covariate, so the linear and the nonlinear parts of a fitted smooth
#' are separately readable.
#'
#' `linear = FALSE` drops that first column, and the penalty is then the
#' identity over the whole block, of full rank.
#'
#' # The construction is empirical
#'
#' The Demmler-Reinsch transform (Demmler and Reinsch, 1975; used for effect
#' selection by Bach and Klein, 2024) takes the inner product **at the observed
#' covariate values**, so it is computed when the term is built and stored in
#' the blueprint. Prediction is the parent basis evaluated at the new points
#' and multiplied by that same transform, so the separation of the linear from
#' the nonlinear part holds at new rows as it does at old.
#'
#' Rebuilding on other rows instead would place the knots on their range and
#' compute another transform. Measured on 80 points, predicting on the first
#' ten agrees with those rows of the original block exactly and rebuilding
#' differs by 2.85.
#'
#' @section Varying the smooth by another variable:
#' A **factor** `by` gives one smooth per level: the block is the smooth
#' multiplied by each level's indicator, and the penalty is the same matrix
#' repeated blockwise, so one smoothing parameter governs every level.
#' `s(x, k = 5, by = g)` over a four-level factor has 16 columns.
#'
#' A **numeric** `by` gives a varying-coefficient term: the smooth multiplies
#' that variable, and the fitted function is the coefficient of `by` as it
#' changes with the covariate.
#'
#' @section Sparse storage:
#' A factor `by` is the one place a smooth's block can be sparse, each row
#' sitting in the block of its own level and nowhere else, a density of
#' \eqn{1/m}. `sparse = TRUE` builds it that way instead of building the dense
#' matrix and compressing it: measured at 2000 rows, \eqn{k = 10} and 200
#' levels, 0.35 MB against 28.93 MB with the numbers identical.
#'
#' `sparse = NULL`, the default, settles it at build from the size of the
#' block through [.resolve_sparse()], the dense form holding \eqn{n m k} cells
#' against \eqn{n k} non-zeros.
#'
#' An explicit `TRUE` is **refused** without a factor `by`, and the message
#' says why: the basis is dense by construction and a numeric `by` merely
#' multiplies it, so there would be nothing to build on.
#'
#' The block alone is sparse. The penalty of a factor `by` is the same matrix
#' repeated blockwise and \pkg{penalties7} returns it dense, 25.92 MB at those
#' sizes; that is a property of that package's contract.
#'
#' @param x The covariate, an expression evaluated in the data.
#' @param by An optional factor or numeric variable, given as a bare
#'   expression; `NULL` by default. See the section above.
#' @param k The basis dimension before reparametrization, `10` by default. It
#'   must exceed `degree`: a cubic spline needs at least four basis functions,
#'   and anything smaller throws. The block has `k - 1` columns with `linear =
#'   TRUE` and `k - 2` without it.
#' @param degree The spline degree, `3` by default, a cubic spline.
#' @param basis An optional \pkg{basis7} basis used in place of the default
#'   B-spline. Its range is taken as given, so a basis built on one interval is
#'   not re-placed on the data's.
#' @param linear Whether the linear effect is carried in the block and left
#'   unpenalized there, `TRUE` by default.
#' @param label A single non-empty string prefixed to the coefficient names.
#'   `NULL`, the default, builds one from the covariate: `s(x)`.
#' @param lambda The smoothing parameter, held at the value given and
#'   **estimated** when left `NULL`, which is the default.
#' @param sparse `TRUE`, `FALSE`, or `NULL` to settle it at build. Only a
#'   factor `by` admits `TRUE`; without one it is refused rather than ignored.
#'   See the section above.
#'
#' @return An unbuilt [SmoothTerm()]: a specification, with `X`, `coef_names`,
#'   `blueprint` and `penalty` empty until [term_build()] fills them.
#'
#' @references
#' Demmler, A. and Reinsch, C. (1975). Oscillation matrices with spline
#' smoothing. *Numerische Mathematik*, 24, 375--382.
#'
#' Bach, P. and Klein, N. (2024). Bayesian effect selection in additive models
#' with an application to time-to-event data.
#'
#' @seealso [te()] for several covariates, [random()] for a grouped effect,
#'   [nl()] for a parametric nonlinear shape, [edf()] for what a fitted smooth
#'   spends.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(80)), g = factor(rep(letters[1:4], 20)))
#' dd$y <- sin(2 * pi * dd$x) + rnorm(80, sd = 0.2)
#'
#' # k = 8 gives seven columns: the linear effect and six deviations.
#' b <- term_build(s(x, k = 8), dd)
#' term_coef_names(b)
#'
#' # The penalty is diag(0, 1, ..., 1): the linear column is free.
#' penalties7::penalty_matrix(term_penalty(b), list(lambda = 1))
#'
#' # The linear column really is the linear effect, and is orthogonal to
#' # the rest over the observed covariate.
#' X <- term_matrix(b)
#' cor(X[, 1], dd$x)
#' max(abs(crossprod(X[, 1], X[, -1])))
#'
#' # So edf runs from k - 1 down to one, not to zero.
#' H <- crossprod(X)
#' cf <- rnorm(ncol(X))
#' vapply(c(1e-8, 1, 1e12),
#'        function(l) edf(b, coef = cf, hessian = H, theta = list(lambda = l)),
#'        numeric(1))
#'
#' # A factor `by` is one smooth per level under one smoothing parameter.
#' bf <- term_build(s(x, k = 5, by = g), dd)
#' c(npar = term_npar(bf), levels = nlevels(dd$g))
#'
#' # The transform is computed on the data and reapplied, never rebuilt.
#' max(abs(term_predict(b, dd[1:10, ]) - X[1:10, ]))
#' max(abs(term_matrix(term_build(s(x, k = 8), dd[1:10, ])) - X[1:10, ]))
#'
#' # Sparsity needs a factor `by`, and says so when there is none.
#' try(term_build(s(x, k = 5, sparse = TRUE), dd))
#'
#'
#' # Fitted. The data are simulated from a known truth, so the
#' # estimates below can be read against it.
#' if (requireNamespace("statmodels7", quietly = TRUE)) {
#'   set.seed(4)
#'   fd <- data.frame(z = sort(runif(200, -3, 3)))
#'   fd$y <- sin(fd$z) + rnorm(200, sd = 0.3)
#'   ft <- statmodels7::statmod(y ~ s(z),
#'                              distributions7::gaussian1_distrib(), fd)
#'   # the smoothing parameter is chosen by REML, and the fit follows sin()
#'   round(c(edf = sum(ft@edf$edf),
#'           rmse = sqrt(mean((fitted(ft) - sin(fd$z))^2))), 3)
#' }
#' @export
s <- function(x, by = NULL, k = 10, degree = 3, basis = NULL,
              linear = TRUE, label = NULL, lambda = NULL,
              sparse = NULL) {
  xe <- substitute(x)
  .smooth_spec(list(xe), substitute(by), k, degree,
               if (is.null(basis)) NULL else list(basis), linear, label,
               sprintf("s(%s)", deparse(xe)), lambda, sparse = sparse)
}

#' Penalized Smooth of Several Covariates
#'
#' @description
#' A tensor-product smooth: a \pkg{basis7} basis in each covariate, combined by
#' [basis7::tensor_basis()], with a roughness penalty on the product
#' coefficients. By default each margin keeps a smoothing parameter of its own,
#' so the surface may be rough in one direction and smooth in another.
#'
#' @details
#' # The block and its penalty
#'
#' The block is the tensor basis evaluated at the covariates, and the penalty
#' is built from the marginal roughness penalties carried into the product,
#' \eqn{P_v = I \otimes \cdots \otimes P_v \otimes \cdots \otimes I}, each
#' penalizing curvature in one direction.
#'
#' With `anisotropic = TRUE`, the default, those components go to
#' [penalties7::additive_penalty()] and keep one smoothing parameter each,
#' named `lambda1`, `lambda2`, and so on. That is the usual reason for fitting
#' a tensor smooth. With `anisotropic = FALSE` they are summed first and one
#' `lambda` governs the total, which costs one hyperparameter instead of one
#' per margin.
#'
#' The marginal bases are **not** reparametrized, so the marginal linear
#' effects are not separated out as [s()] separates its one. They lie in the
#' null space of the tensor penalty, and a strongly penalized fit is shrunk
#' toward no surface at all rather than toward a plane.
#'
#' @section Centering:
#' The tensor product of the marginal bases contains the constant, which the
#' penalty's null space contains as well, so beside an intercept the block
#' would be rank deficient by exactly one and the penalty would not cover the
#' deficiency.
#'
#' The block therefore carries the sum-to-zero constraint over the observed
#' covariates ([basis7::constrain_basis()]). The term has **one column fewer**
#' than the product of its marginal dimensions, so `te(x, z, k = 4)` gives 15
#' and not 16; every column sums to zero over the data it was built on, to
#' machine precision;
#' and the penalty follows by congruence with its rank unchanged, the direction
#' removed having been one of its null directions.
#'
#' The transform is stored in the blueprint and reapplied by [term_predict()],
#' as the Demmler-Reinsch transform of [s()] is.
#'
#' The level of the surface is then the model's intercept, so a formula that
#' removes it, `y ~ te(x, z) - 1`, fits a surface constrained to average zero
#' over the data.
#'
#' @param ... The covariates, bare expressions evaluated in the data, at least
#'   two of them. One throws `"'te' needs at least two covariates; use s() for
#'   one."`.
#' @param by An optional factor or numeric variable, as in [s()], with the same
#'   two readings and the same sparsity rule.
#' @param k The basis dimension per margin, `5` by default, recycled to the
#'   number of covariates. As in [s()] it must exceed `degree`.
#' @param degree The spline degree per margin, `3` by default, recycled.
#' @param bases An optional list of \pkg{basis7} bases, one per covariate, used
#'   in place of the default B-splines.
#' @param anisotropic `TRUE`, the default, for one smoothing parameter per
#'   margin; `FALSE` for one over their sum. Anything that is not a single
#'   logical throws.
#' @param label A single non-empty string prefixed to the coefficient names.
#'   `NULL`, the default, builds one from the covariates: `te(x,z)`.
#' @param lambda The smoothing parameters, held at the values given and
#'   **estimated** when left `NULL`, which is the default. An anisotropic
#'   product carries one per margin, so a vector of that length, or a named one
#'   holding some of them.
#' @param sparse `TRUE`, `FALSE`, or `NULL` to settle it at build. Only a
#'   factor `by` admits `TRUE`. See [s()].
#'
#' @return An unbuilt [SmoothTerm()]: a specification, with `X`, `coef_names`,
#'   `blueprint` and `penalty` empty until [term_build()] fills them.
#'
#' @references
#' Wood, S. N. (2006). Low-rank scale-invariant tensor product smooths for
#' generalized additive mixed models. *Biometrics* 62, 1025--1036.
#'
#' Wood, S. N. (2017). *Generalized Additive Models: An Introduction with R*,
#' 2nd edition. Chapman and Hall/CRC.
#'
#' @seealso [s()] for one covariate, [basis7::tensor_basis()] for the product,
#'   [penalties7::additive_penalty()] for the anisotropic penalty.
#'
#' @examples
#' set.seed(2)
#' dd <- data.frame(x = runif(120), z = runif(120))
#' dd$y <- dd$x * dd$z + rnorm(120, sd = 0.1)
#'
#' # Four by four margins give fifteen columns: the centering removes one.
#' b <- term_build(te(x, z, k = 4), dd)
#' c(npar = term_npar(b), product = 4 * 4)
#'
#' # Every column sums to zero over the data, to machine precision.
#' max(abs(colSums(term_matrix(b))))
#'
#' # Anisotropic by default: one smoothing parameter per margin.
#' term_penalty(b)@penalty_name
#' term_penalty(b)@params
#' term_penalty(term_build(te(x, z, k = 4, anisotropic = FALSE), dd))@params
#'
#' # Holding both of them.
#' term_hyper(te(x, z, k = 4, lambda = c(1, 5)))
#'
#' # The centering transform is reapplied, not recomputed.
#' max(abs(term_predict(b, dd[1:10, ]) - term_matrix(b)[1:10, ]))
#'
#' # One covariate is s(), not te().
#' try(te(x, k = 4))
#'
#'
#' # Fitted. The data are simulated from a known truth, so the
#' # estimates below can be read against it.
#' if (requireNamespace("statmodels7", quietly = TRUE)) {
#'   set.seed(5)
#'   fd <- data.frame(a = runif(300, -2, 2), b = runif(300, -2, 2))
#'   fd$y <- sin(fd$a) * cos(fd$b) + rnorm(300, sd = 0.3)
#'   ft <- statmodels7::statmod(y ~ te(a, b, k = 5),
#'                              distributions7::gaussian1_distrib(), fd)
#'   # one smoothing parameter per margin, against a truth of sin(a) cos(b)
#'   round(sqrt(mean((fitted(ft) - sin(fd$a) * cos(fd$b))^2)), 3)
#' }
#' @export
te <- function(..., by = NULL, k = 5, degree = 3, bases = NULL,
               anisotropic = TRUE, label = NULL, lambda = NULL,
               sparse = NULL) {
  vars <- as.list(substitute(list(...)))[-1L]
  if (length(vars) < 2L) {
    stop("'te' needs at least two covariates; use s() for one.",
         call. = FALSE)
  }
  if (!is.logical(anisotropic) || length(anisotropic) != 1L ||
      is.na(anisotropic)) {
    stop("'anisotropic' must be TRUE or FALSE.", call. = FALSE)
  }
  # ANISOTROPIC is one smoothing parameter per margin, so that is how many
  # names there are to hold; isotropic is the single one of a quadratic
  # penalty. Either way the names are the penalty's own, which is what the
  # summary prints and what the fit keys its hyperparameters by.
  nms <- if (anisotropic) paste0("lambda", seq_along(vars)) else "lambda"
  sp <- .smooth_spec(vars, substitute(by), k, degree, bases, FALSE, label,
                     sprintf("te(%s)",
                             paste(vapply(vars, deparse, character(1)),
                                   collapse = ",")),
                     lambda, nms, sparse = sparse)
  sp@spec$anisotropic <- anisotropic
  sp
}

.smooth_spec <- function(vars, by, k, degree, bases, linear, label,
                         default_label, lambda = NULL, names = "lambda",
                         sparse = NULL) {
  nv <- length(vars)
  chk <- function(v, nm, lo) {
    if (!is.numeric(v) || anyNA(v) || any(v < lo) || any(v != round(v))) {
      stop(sprintf("'%s' must be a whole number of at least %d.", nm, lo),
           call. = FALSE)
    }
    as.integer(rep_len(v, nv))
  }
  k <- chk(k, "k", 3L)
  degree <- chk(degree, "degree", 1L)
  # a B-spline of degree m spans at least m + 1 functions, and finding that
  # out inside term_build() would report it several frames from the call
  if (is.null(bases) && any(k <= degree)) {
    stop(sprintf(paste("'k' must exceed 'degree': a spline of degree %d needs",
                       "at least %d basis functions."),
                 max(degree), max(degree) + 1L), call. = FALSE)
  }
  if (!is.null(bases)) {
    if (!is.list(bases) || length(bases) != nv ||
        !all(vapply(bases, function(b) S7::S7_inherits(b, basis7::basis),
                    logical(1)))) {
      stop("'bases' must be a list of one basis7 basis per covariate.",
           call. = FALSE)
    }
  }
  if (is.null(label)) label <- default_label
  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !nzchar(label)) {
    stop("'label' must be a single non-empty character string.",
         call. = FALSE)
  }
  if (!is.null(sparse) &&
      (!is.logical(sparse) || length(sparse) != 1L || is.na(sparse))) {
    stop(paste0("'sparse' must be TRUE, FALSE, or NULL to settle it from",
                " the design."), call. = FALSE)
  }
  # ASKED AT CONSTRUCTION, where the caller can see it: the sparsity of a
  # smooth comes from a factor `by`, whose indicators put each row in one
  # block. Without a `by` there is only the basis, which is dense by
  # construction, so the argument would have nothing to build on.
  if (isTRUE(sparse) && is.null(by)) {
    stop(paste0("'sparse' has nothing to build on here: a smooth's basis is",
                " dense by\n  construction. Sparsity comes from a FACTOR",
                " 'by', whose indicators put each\n  row in the block of its",
                " own level."), call. = FALSE)
  }
  SmoothTerm(label = label, vars = vars, by = by, sparse = sparse,
             spec = list(k = k, degree = degree, bases = bases,
                         linear = isTRUE(linear)),
             hyper = smooth_hyper(lambda, names, label),
             X = NULL, coef_names = character(0),
             blueprint = list(), penalty = NULL)
}

# the covariate values, one column per variable
.smooth_x <- function(vars, data) {
  out <- lapply(vars, function(e) as.numeric(eval(e, data, baseenv())))
  n <- unique(lengths(out))
  if (length(n) != 1L || n != nrow(data)) {
    stop("every covariate of a smooth must give one value per row.",
         call. = FALSE)
  }
  if (any(vapply(out, anyNA, logical(1)))) {
    stop("a smooth's covariates must not contain missing values.",
         call. = FALSE)
  }
  out
}

# One copy of the basis per level of a factor `by`: row i is non-zero only
# in the block of its own level, which is the shape .random_block() has and
# a density of 1/m. Built sparse where the term asks for it, rather than
# built dense and compressed.
.smooth_by_block <- function(Z, g, m, sparse = FALSE) {
  n <- nrow(Z)
  cols <- ncol(Z)
  gi <- as.integer(g)
  if (!isTRUE(sparse)) {
    G <- outer(g, levels(g), `==`) * 1
    Zb <- matrix(0, n, m * cols)
    for (l in seq_len(m)) {
      Zb[, (l - 1L) * cols + seq_len(cols)] <- G[, l] * Z
    }
    return(Zb)
  }
  if (anyNA(gi)) {
    stop(paste0("a 'by' with missing values has no block to place those rows",
                " in;\n  sparse = TRUE cannot build them."), call. = FALSE)
  }
  Matrix::sparseMatrix(
    i = rep(seq_len(n), each = cols),
    j = as.vector(t(outer((gi - 1L) * cols, seq_len(cols), `+`))),
    x = as.vector(t(as.matrix(Z))),
    dims = c(n, m * cols))
}

# the by variable: a factor gives one copy of the block per level, a numeric
# multiplies it
.smooth_by <- function(by, data, levels = NULL) {
  if (is.null(by)) return(NULL)
  v <- eval(by, data, baseenv())
  if (is.factor(v) || is.character(v) || is.logical(v)) {
    f <- if (is.null(levels)) factor(v) else factor(v, levels = levels)
    if (any(is.na(f) & !is.na(v))) {
      stop("a 'by' level absent at build time cannot be predicted.",
           call. = FALSE)
    }
    list(kind = "factor", value = f, levels = levels(f))
  } else {
    list(kind = "numeric", value = as.numeric(v))
  }
}

#' @title Build a Smooth Term
#' @name term_build.SmoothTerm
#'
#' @description
#' Builds the basis of an [s()] or [te()] term at the observed covariates,
#' applies the reparametrization that construction calls for, multiplies by a
#' `by` variable where there is one, and attaches the roughness penalty. Two
#' quantities are computed **from the data** here and recorded in the
#' blueprint, so that [term_predict()] reapplies them instead of deriving them
#' again: the Demmler-Reinsch transform for [s()], and the centering constraint
#' for [te()].
#'
#' @details
#' # The default basis
#'
#' Where no basis was supplied, each margin gets a
#' [basis7::bspline_basis()] over the observed range of its covariate, padded
#' by a thousandth of that range at each end so that the extreme observations
#' are strictly inside. A basis given through `basis` or `bases` is used with
#' its own range, untouched.
#'
#' # One covariate
#'
#' [basis7::dr_basis()] takes the inner product at the observed values, and the
#' block is that basis evaluated there. With `linear = TRUE` a column
#' \eqn{(x - \bar{x})/s_x} is prepended and the penalty is
#' \eqn{\mathrm{diag}(0, 1, \dots, 1)}; the center and the scale go into the
#' blueprint with the transform. Without it the penalty is the identity.
#'
#' # Several covariates
#'
#' [basis7::tensor_basis()] combines the margins, and one penalty component per
#' margin is carried into the product as \eqn{I \otimes \cdots \otimes P_v
#' \otimes \cdots \otimes I}, each \eqn{P_v} the margin's second-derivative
#' Gram normalized to a maximum entry of one.
#'
#' The Kronecker product is taken over the **reversed** blocks, because
#' [basis7::tensor_basis()] varies the first margin fastest. The block is then
#' centered by [basis7::constrain_basis()] over the observed covariates, so it
#' has one column fewer than the product of the marginal dimensions.
#'
#' # The `by` variable
#'
#' A factor `by` interacts the basis with the level indicators, giving `m`
#' copies of the block and a penalty repeated blockwise; the levels are
#' recorded so a prediction uses the same set. A numeric `by` multiplies the
#' basis. The storage is settled here and recorded, an explicit `sparse = TRUE`
#' being refused where there is no factor `by`.
#'
#' @param term An unbuilt or built [SmoothTerm()].
#' @param data A data frame carrying the covariates and the `by` variable.
#' @param ... Unused.
#'
#' @return The term with `X`, `coef_names`, `blueprint` and `penalty` filled.
#'
#' @seealso [s()] and [te()] for the two constructions,
#'   [term_predict.SmoothTerm()] for the block at new rows.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(80)), z = runif(80))
#'
#' # What the build records for a one-covariate smooth.
#' b <- term_build(s(x, k = 8), dd)
#' names(b@blueprint)
#' b@blueprint$core$kind
#'
#' # And for a tensor product: one column fewer than 4 x 4.
#' bt <- term_build(te(x, z, k = 4), dd)
#' term_npar(bt)
#'
#' @keywords internal
S7::method(term_build, SmoothTerm) <- function(term, data, ...) {
  xs <- .smooth_x(term@vars, data)
  sp <- term@spec
  nv <- length(xs)

  marg <- lapply(seq_len(nv), function(j) {
    if (!is.null(sp$bases)) return(sp$bases[[j]])
    r <- range(xs[[j]])
    pad <- diff(r) * 0.001 + .Machine$double.eps
    basis7::bspline_basis(lower = r[1] - pad, upper = r[2] + pad,
                          dimension = sp$k[j], degree = sp$degree[j])
  })

  if (nv == 1L) {
    b <- marg[[1L]]
    d <- basis7::dr_basis(b, xs[[1L]])
    Z <- basis7::basis_eval(d, xs[[1L]])
    nm <- paste0("z", seq_len(ncol(Z)))
    if (sp$linear) {
      ctr <- mean(xs[[1L]])
      scl <- stats::sd(xs[[1L]])
      if (!is.finite(scl) || scl == 0) scl <- 1
      Z <- cbind((xs[[1L]] - ctr) / scl, Z)
      nm <- c("lin", nm)
      pen_diag <- c(0, rep(1, ncol(Z) - 1L))
      lin <- list(center = ctr, scale = scl)
    } else {
      pen_diag <- rep(1, ncol(Z))
      lin <- NULL
    }
    P <- diag(pen_diag, length(pen_diag))
    core <- list(kind = "dr", dr = d, lin = lin)
  } else {
    tb <- basis7::tensor_basis(marg)
    xm <- do.call(cbind, xs)
    # the marginal roughness penalties carried into the product: the second
    # derivative Gram of each margin, the identity in the others
    dims <- vapply(marg, function(b) b@dimension, integer(1))
    comps <- lapply(seq_len(nv), function(j) {
      Pj <- basis7::basis_gram(marg[[j]], order = 2L)
      Pj <- Pj / max(1, max(abs(Pj)))
      blocks <- lapply(seq_len(nv), function(i)
        if (i == j) Pj else diag(dims[i]))
      # tensor_basis varies the FIRST margin fastest, so the Kronecker
      # product is taken in reverse order
      Reduce(kronecker, rev(blocks))
    })
    # The tensor product contains the constant, which the penalty's null space
    # also contains, so a model carrying an intercept would be rank deficient
    # by exactly one and nothing in the penalty would cover it. The sum-to-zero
    # constraint over the observed covariates removes that direction: every
    # column of the constrained block sums to zero, hence is orthogonal to an
    # intercept column. The penalty follows by congruence, and since the
    # direction removed lies in its null space the rank is unchanged.
    tb <- basis7::constrain_basis(tb, colSums(basis7::basis_eval(tb, xm)))
    tmat <- tb@transform
    comps <- lapply(comps, function(Pk) {
      M <- crossprod(tmat, Pk %*% tmat)
      (M + t(M)) / 2
    })
    Z <- basis7::basis_eval(tb, xm)
    nm <- paste0("z", seq_len(ncol(Z)))
    P <- if (isTRUE(sp$anisotropic)) comps else Reduce(`+`, comps)
    core <- list(kind = "tensor", basis = tb)
  }

  by_blocks <- 1L
  # settled here and carried to the blueprint: without a factor `by` there is
  # nothing sparse to build, so the answer is FALSE whatever was asked
  sp <- FALSE
  by <- .smooth_by(term@by, data)
  if (!is.null(by)) {
    if (identical(by$kind, "factor")) {
      m <- nlevels(by$value)
      nm <- as.character(t(outer(levels(by$value), nm,
                                 function(a, b) paste(a, b, sep = "."))))
      # the block is m copies of a basis of ncol(Z), one row per level, so
      # the dense form holds n * m * ncol(Z) cells against n * ncol(Z)
      # non-zeros -- the same product .resolve_sparse() reads for a formula
      sp <- .resolve_sparse(term@sparse, nrow(Z), m * ncol(Z))
      Z <- .smooth_by_block(Z, by$value, m, sp)
      # ONE COPY PER LEVEL, and the penalty says so rather than being handed
      # the assembled product. What a quadratic penalty needs from its matrix
      # -- the rank, the log pseudo-determinant, a basis of the null space --
      # all follow from one block, the eigenvalues of I (x) P being P's
      # repeated m times, so nothing of size (mk)^2 is ever decomposed.
      # Measured at m = 200 over a basis of ten, 565 times faster to build
      # and 0.03 MB stored against 32.
      #
      # The ANISOTROPIC branch still assembles: additive_penalty() reads one
      # eigendecomposition of the SUM of its components, which is not a
      # blockwise quantity, so the same shortcut does not apply to it.
      by_blocks <- m
      if (is.list(P)) {
        P <- lapply(P, function(Pk) kronecker(diag(m), Pk))
        by_blocks <- 1L
      }
    } else {
      # a numeric `by` MULTIPLIES the basis, so the block is as dense as the
      # basis is, and there is nothing for a sparse storage to hold back
      if (isTRUE(term@sparse)) {
        stop(paste0("'sparse' has nothing to build on here: a numeric 'by'",
                    " multiplies the basis,\n  so the block is as dense as",
                    " the basis. Sparsity comes from a FACTOR 'by', whose",
                    "\n  indicators put each row in one block."),
             call. = FALSE)
      }
      Z <- by$value * Z
    }
  }

  cn <- paste(term@label, nm, sep = ".")
  colnames(Z) <- cn
  rownames(Z) <- NULL
  term@X <- Z
  term@coef_names <- cn
  term@blueprint <- list(core = core, marg = marg, spec = sp,
                         vars = term@vars, by = term@by,
                         by_levels = if (!is.null(by) &&
                                         identical(by$kind, "factor"))
                           by$levels else NULL,
                         sparse = sp,
                         nblock = ncol(Z))
  term@penalty <- if (is.list(P)) penalties7::additive_penalty(P)
                  else penalties7::quadratic_penalty(P, blocks = by_blocks)
  term
}

#' @title A Smooth Term's Block at New Rows
#' @name term_predict.SmoothTerm
#'
#' @description
#' Evaluates a built smooth's recorded basis at the covariates in `newdata` and
#' applies the transforms the build computed: the Demmler-Reinsch basis and its
#' centering and scaling for [s()], the centered tensor basis for [te()], then
#' the `by` variable against the recorded levels. Nothing is derived from the
#' new rows.
#'
#' @details
#' The basis object in `blueprint$core` carries the knots and the transform, so
#' the columns at new rows are the same functions of the covariate as the
#' fitted ones, so \eqn{\tilde{X}\beta} is the fitted smooth evaluated
#' there. Rebuilding instead would place the knots on the new range
#' and compute another transform: measured on 80 points, predicting on the
#' first ten agrees with those rows of the block exactly and rebuilding differs
#' by 2.85.
#'
#' A factor `by` is expanded against `blueprint$by_levels`, so `newdata` need
#' carry only the levels its own rows use and still gets every column, and the
#' block is built in the storage the build settled on.
#'
#' New covariate values outside the range the basis was placed on are
#' evaluated, not refused. A B-spline is zero beyond its knots, so the fitted
#' function flattens rather than extrapolating a trend; read a prediction far
#' outside the fitting range with that in mind.
#'
#' @param term A built [SmoothTerm()]. An unbuilt one throws
#'   `"the term has not been built; call term_build(term, data) first."`.
#' @param newdata A data frame carrying the covariates and the `by` variable.
#' @param ... Unused.
#'
#' @return A block of `nrow(newdata)` rows and [term_npar()] columns, in the
#'   storage the build settled on, with the term's coefficient names as column
#'   names and no row names.
#'
#' @seealso [term_predict()] for the generic and the identity it satisfies,
#'   [term_build.SmoothTerm()] for what recorded the transform.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(80)), g = factor(rep(letters[1:4], 20)))
#' b <- term_build(s(x, k = 8), dd)
#' X <- term_matrix(b)
#'
#' # Reapplying is exact; rebuilding on the same rows is a different basis.
#' max(abs(term_predict(b, dd[1:10, ]) - X[1:10, ]))
#' max(abs(term_matrix(term_build(s(x, k = 8), dd[1:10, ])) - X[1:10, ]))
#'
#' # A factor `by` keeps every level's columns at a subset that has two.
#' bf <- term_build(s(x, k = 5, by = g), dd)
#' nd <- droplevels(dd[dd$g %in% c("a", "b"), ])
#' c(levels_here = nlevels(nd$g), cols = ncol(term_predict(bf, nd)))
#'
#' @keywords internal
S7::method(term_predict, SmoothTerm) <- function(term, newdata, ...) {
  .assert_built(term)
  bp <- term@blueprint
  xs <- .smooth_x(bp$vars, newdata)
  core <- bp$core

  if (identical(core$kind, "dr")) {
    Z <- basis7::basis_eval(core$dr, xs[[1L]])
    if (!is.null(core$lin)) {
      Z <- cbind((xs[[1L]] - core$lin$center) / core$lin$scale, Z)
    }
  } else {
    Z <- basis7::basis_eval(core$basis, do.call(cbind, xs))
  }

  by <- .smooth_by(bp$by, newdata, levels = bp$by_levels)
  if (!is.null(by)) {
    if (identical(by$kind, "factor")) {
      # the STORAGE is part of the blueprint: a prediction that densified
      # would spend at new data what the build was careful not to
      Z <- .smooth_by_block(Z, by$value, length(bp$by_levels),
                            isTRUE(bp$sparse))
    } else {
      Z <- by$value * Z
    }
  }
  colnames(Z) <- term@coef_names
  rownames(Z) <- NULL
  Z
}
