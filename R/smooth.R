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
#' `spec` carries what the build reads: the \pkg{basis7} smoothers, one per
#' covariate, which hold the whole construction; `by_hyper`, saying whether
#' the levels of a factor `by` share a smoothing parameter or carry one each;
#' and for [te()] whether the penalty is anisotropic. What a build then
#' computes from the data goes into the blueprint instead: the
#' Demmler-Reinsch transform, the centering constraint, the `by` levels and
#' the penalties the term declares.
#'
#' `sparse` is `NULL` until the build settles it. A smooth's block is sparse
#' only under a **factor** `by`, where each row sits in the block of its own
#' level: the basis itself is dense by construction, the Demmler-Reinsch
#' rotation making it so, and a numeric `by` merely multiplies it.
#'
#' @inheritParams additive_term
#' @inheritParams model_term
#' @param vars A list of the covariate expressions being smoothed.
#' @param by The expression the smooth varies with, or `NULL`.
#' @param spec A named list of what the build reads: `smoothers`, one
#'   \pkg{basis7} smoother per covariate; `by_hyper`; and for [te()] the
#'   `anisotropic` flag.
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
#' tm <- s(x, basis7::bspline_smooth(k = 8))
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
#' followed by the reparametrized basis, so `s(x, bspline_smooth(k = 8))` gives seven columns
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
#' `bspline_smooth(null_space = "drop")` drops that first column, and the
#' penalty is then the identity over the whole block, of full rank.
#'
#' # A penalty of your own
#'
#' The smoother's `penalty` argument replaces the roughness matrix with a
#' \pkg{penalties7} penalty of its own, built at the coefficient count only
#' the data settle: `s(x, bspline_smooth(k = 20, penalty =
#' penalties7::lasso_penalty))`. A \pkg{penalties7} constructor passes bare,
#' and anything else is a function of the count.
#'
#' What it buys is that the Demmler-Reinsch coordinates are ordered from the
#' smoothest to the most wiggly and the roughness penalty on them is the
#' identity, so an \eqn{\ell_1} penalty takes whole directions of wiggliness
#' to exactly zero and chooses the smooth's effective dimension. Measured at
#' \eqn{n = 300} with \eqn{k = 20}, the coordinates surviving as the
#' smoothing parameter grows are 15, 10, 4, 2, 1 and 0 of 18, and the fit at
#' FOUR of them is the closest to the truth: the root mean square errors at
#' those six counts are 0.2453, 0.2427, 0.2412, 0.2425, 0.2709 and 0.4967. So
#' the \eqn{\ell_1} buys a smaller smooth and a slightly better one at once.
#' A heavy-tailed penalty is the robust reading of the same block, and a
#' structured one estimates the correlation of the coefficients rather than
#' fixing it.
#'
#' Two things follow, and neither is a detail. The penalty covers **the
#' penalized coordinates alone**: a roughness matrix carries a zero row for
#' the linear column and leaves it free by its own arithmetic, where a
#' separable penalty has no such row and would shrink it, so the free columns
#' are outside the entry it declares. And a penalty **with a kink** has no
#' second derivative at zero, which a marginal criterion needs, so its
#' smoothing parameter is chosen by a path over `sparse_criterion` -- BIC by
#' default -- and not by REML, whatever `outer_criterion` says. That is a
#' change in what the number means and not only in how it is computed.
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
#' repeated blockwise.
#' `s(x, bspline_smooth(k = 5), by = g)` over a four-level factor has 16 columns.
#'
#' `by_hyper` says whether the levels are smoothed **together or apart**.
#' `"shared"`, the default, estimates one smoothing parameter for all of them,
#' which is what this package has always done. `"level"` declares one penalty
#' per level through [term_penalties()], so each is estimated on its own; it
#' is what \pkg{mgcv} does by default, and what `id` there undoes. Which to
#' want is a question about the data rather than about the code: levels that
#' differ in how wiggly they are, or in how much of them there is, ask for
#' different amounts of smoothing. Measured on three levels of one curve at
#' amplitudes 1, 0.15 and 2.5, the separately estimated smoothing parameters
#' are 2.03, 111.9 and 0.468 -- a factor of 239 apart -- and the effective
#' degrees of freedom fall from 25.5 to 21.7 as the flattest level is smoothed
#' away.
#'
#' Under `"level"` a held hyperparameter is qualified by the level it belongs
#' to, `hyper = c(lambda.a = 2)`, and so is an `id`.
#'
#' A **numeric** `by` gives a varying-coefficient term: the smooth multiplies
#' that variable, and the fitted function is the coefficient of `by` as it
#' changes with the covariate. It has no levels, so `by_hyper = "level"` is
#' rejected rather than read as `"shared"`.
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
#' @param smoother How the smooth is built: a \pkg{basis7} smoother, which
#'   carries the basis, the roughness penalty, the null space and the
#'   reparametrization together. [basis7::bspline_smooth()] is the default
#'   and is the construction this term has always used;
#'   [basis7::fourier_smooth()] is the periodic one and
#'   [basis7::legendre_smooth()] the global polynomial one.
#' @param by An optional factor or numeric variable, given as a bare
#'   expression; `NULL` by default. See the section above.
#' @param by_hyper With a factor `by`, whether the levels share one smoothing
#'   parameter (`"shared"`, the default) or carry one each (`"level"`). See
#'   the section above. Rejected without a factor `by`, which has no levels
#'   to give one each.
#' @param hyper The hyperparameters of the smoother's penalty to hold, as a
#'   named numeric vector such as `c(lambda = 2)`. What names there are
#'   depends on the penalty; anything left out is **estimated**, which is
#'   the default. Under `by_hyper = "level"` a name is qualified by the
#'   level, `c(lambda.a = 2)`.
#' @param label A single non-empty string prefixed to the coefficient names.
#'   `NULL`, the default, builds one from the covariate: `s(x)`.
#' @param id A label sharing this smooth's smoothing parameter with those
#'   of other terms carrying the same one: they are then estimated at a
#'   single value, so several curves are smoothed together. It is what
#'   `id` does in \pkg{mgcv}, and what it means best between smooths of
#'   the same basis and dimension. `NULL`, the default, shares nothing.
#'   See [term_ids()].
#' @param sparse `TRUE`, `FALSE`, or `NULL` to settle it at build. Only a
#'   factor `by` admits `TRUE`; without one it is refused rather than ignored.
#'   See the section above.
#' @param ... Not used, and accepted only so that an argument `s()` no
#'   longer takes is reported with its replacement rather than as R's own
#'   "unused argument", which names the argument and not what to write.
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
#' b <- term_build(s(x, basis7::bspline_smooth(k = 8)), dd)
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
#' bf <- term_build(s(x, basis7::bspline_smooth(k = 5), by = g), dd)
#' c(npar = term_npar(bf), levels = nlevels(dd$g))
#' length(term_penalties(bf))
#'
#' # by_hyper = "level" gives them one each: one entry per level, over that
#' # level's own columns.
#' bl <- term_build(s(x, basis7::bspline_smooth(k = 5), by = g,
#'                    by_hyper = "level"), dd)
#' vapply(term_penalties(bl), function(e) e$name, character(1))
#' vapply(term_penalties(bl), function(e) range(e$index), integer(2))
#'
#' # A penalty factory covers the PENALIZED coordinates: the free linear
#' # column is left out, a separable penalty having no zero row with which
#' # to leave it alone.
#' bp <- term_build(s(x, basis7::bspline_smooth(
#'   k = 5, penalty = penalties7::lasso_penalty)), dd)
#' term_penalties(bp)[[1]]$index
#' term_penalties(bp)[[1]]$penalty@params
#'
#' # The transform is computed on the data and reapplied, never rebuilt.
#' max(abs(term_predict(b, dd[1:10, ]) - X[1:10, ]))
#' max(abs(term_matrix(term_build(s(x, basis7::bspline_smooth(k = 8)), dd[1:10, ])) - X[1:10, ]))
#'
#' # Sparsity needs a factor `by`, and says so when there is none.
#' try(term_build(s(x, basis7::bspline_smooth(k = 5), sparse = TRUE), dd))
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
s <- function(x, smoother = basis7::bspline_smooth(), by = NULL,
              by_hyper = c("shared", "level"), hyper = NULL, id = NULL,
              label = NULL, sparse = NULL, ...) {
  xe <- substitute(x)
  retired_smooth_args(list(...), "s")
  .smooth_spec(list(xe), substitute(by), list(smoother), label,
               sprintf("s(%s)", deparse(xe)), hyper, sparse = sparse,
               ids = id, by_hyper = by_hyper)
}


#' Report an Argument That s() and te() No Longer Take
#'
#' @description
#' Names the replacement for each of the arguments the smooth constructors
#' carried before the construction became an object. They are reported here
#' rather than left to R's own "unused argument", which names the argument
#' and not what to write instead.
#'
#' @details
#' The four decisions a smooth is made of -- the basis, the penalty, the
#' null space and the reparametrization -- now travel together on a
#' \pkg{basis7} smoother, so the settings that described the first two
#' belong there. `k` and `degree` are a B-spline's and were never general;
#' `linear` named the null space as "linear", which it is only when the
#' penalty is of order 2.
#'
#' @param dots The `...` of the constructor.
#' @param fn `"s"` or `"te"`, for the message.
#'
#' @return `NULL`, invisibly. Called for the error it signals.
#'
#' @keywords internal
retired_smooth_args <- function(dots, fn) {
  nm <- names(dots)
  if (is.null(nm) || !length(nm)) {
    if (length(dots)) {
      stop(sprintf("'%s' takes no unnamed argument beyond its covariate.", fn),
           call. = FALSE)
    }
    return(invisible(NULL))
  }
  arg <- if (identical(fn, "s")) "smoother" else "smooths"
  moved <- c(
    k = sprintf("%s = bspline_smooth(k = ...)", arg),
    degree = sprintf("%s = bspline_smooth(degree = ...)", arg),
    basis = sprintf("%s = a smoother, or bspline_smooth(lower =, upper =)",
                    arg),
    bases = sprintf("%s = a list of smoothers, one per covariate", arg),
    linear = sprintf("%s = bspline_smooth(null_space = \"keep\"/\"drop\")",
                     arg),
    lambda = "hyper = c(lambda = ...)"
  )
  hit <- intersect(nm, names(moved))
  if (length(hit)) {
    stop(sprintf(paste0(
      "'%s' no longer takes '%s'. The construction of a smooth is now an",
      " object:\n  write %s.\n  See ?basis7::bspline_smooth."
    ), fn, hit[[1L]], moved[[hit[[1L]]]]), call. = FALSE)
  }
  stop(sprintf("'%s' has no argument '%s'.", fn, nm[nzchar(nm)][[1L]]),
       call. = FALSE)
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
#' than the product of its marginal dimensions, so
#' `te(x, z, smooths = bspline_smooth(k = 4))` gives 15
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
#' @param by_hyper With a factor `by`, whether the levels share the term's
#'   smoothing parameters (`"shared"`, the default) or carry a set each
#'   (`"level"`), as in [s()]. Under `anisotropic = TRUE` a set is one per
#'   margin, so `"level"` gives one per margin per level.
#' @param smooths How each margin is built: one \pkg{basis7} smoother used
#'   for every covariate, or a list of one per covariate.
#'   `bspline_smooth(k = 5)` is the default. A margin carrying a `penalty`
#'   factory is rejected: a tensor product's penalty is the sum of the
#'   marginal roughnesses, and a factory does not say how it composes with
#'   that sum.
#'
#'   What the product reads from a margin is its **basis** and its
#'   **roughness matrix**, so `k`, `degree`, `order`, `measure`, `lower` and
#'   `upper` all act. The constraint, the null space and the
#'   reparametrization are the product's own and not a margin's: the tensor
#'   contains the constant whatever its margins do, and it is the product
#'   that is centered. A margin asking for one of those is rejected rather
#'   than ignored.
#' @param anisotropic `TRUE`, the default, for one smoothing parameter per
#'   margin; `FALSE` for one over their sum. Anything that is not a single
#'   logical throws.
#' @param label A single non-empty string prefixed to the coefficient names.
#'   `NULL`, the default, builds one from the covariates: `te(x,z)`.
#' @param hyper The hyperparameters to hold, as a named numeric vector, with
#'   anything left out **estimated**. An anisotropic product carries one per
#'   margin, so a vector of that length, or a named one holding some of
#'   them.
#' @param id Labels sharing this term's smoothing parameters with those of
#'   other terms carrying the same ones. An anisotropic product carries
#'   one per margin, so the labels are named after them,
#'   `c(lambda1 = "A")`; an isotropic one carries `lambda` alone and a
#'   single unnamed string will do. `NULL`, the default, shares nothing.
#'   See [term_ids()].
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
#' b <- term_build(te(x, z, smooths = basis7::bspline_smooth(k = 4)), dd)
#' c(npar = term_npar(b), product = 4 * 4)
#'
#' # Every column sums to zero over the data, to machine precision.
#' max(abs(colSums(term_matrix(b))))
#'
#' # Anisotropic by default: one smoothing parameter per margin.
#' term_penalty(b)@penalty_name
#' term_penalty(b)@params
#' ti <- te(x, z, smooths = basis7::bspline_smooth(k = 4),
#'           anisotropic = FALSE)
#' term_penalty(term_build(ti, dd))@params
#'
#' # Holding both of them.
#' term_hyper(te(x, z, smooths = basis7::bspline_smooth(k = 4), hyper = c(1, 5)))
#'
#' # The centering transform is reapplied, not recomputed.
#' max(abs(term_predict(b, dd[1:10, ]) - term_matrix(b)[1:10, ]))
#'
#' # One covariate is s(), not te().
#' try(te(x, smooths = basis7::bspline_smooth(k = 4)))
#'
#'
#' # Fitted. The data are simulated from a known truth, so the
#' # estimates below can be read against it.
#' if (requireNamespace("statmodels7", quietly = TRUE)) {
#'   set.seed(5)
#'   fd <- data.frame(a = runif(300, -2, 2), b = runif(300, -2, 2))
#'   fd$y <- sin(fd$a) * cos(fd$b) + rnorm(300, sd = 0.3)
#'   ft <- statmodels7::statmod(y ~ te(a, b, smooths = basis7::bspline_smooth(k = 5)),
#'                              distributions7::gaussian1_distrib(), fd)
#'   # one smoothing parameter per margin, against a truth of sin(a) cos(b)
#'   round(sqrt(mean((fitted(ft) - sin(fd$a) * cos(fd$b))^2)), 3)
#' }
#' @export
te <- function(..., smooths = basis7::bspline_smooth(k = 5), by = NULL,
               by_hyper = c("shared", "level"), anisotropic = TRUE,
               hyper = NULL, id = NULL, label = NULL, sparse = NULL) {
  vars <- as.list(substitute(list(...)))[-1L]
  # THE COVARIATES ARE THE UNNAMED ARGUMENTS, so a retired one such as
  # `k = 4` would otherwise be taken for a covariate called "k" -- silently,
  # since te() has no dots left to catch it in.
  vn <- names(vars)
  if (!is.null(vn) && any(nzchar(vn))) {
    retired_smooth_args(as.list(vars)[nzchar(vn)], "te")
  }
  if (length(vars) < 2L) {
    stop("'te' needs at least two covariates; use s() for one.",
         call. = FALSE)
  }
  if (!is.logical(anisotropic) || length(anisotropic) != 1L ||
      is.na(anisotropic)) {
    stop("'anisotropic' must be TRUE or FALSE.", call. = FALSE)
  }
  smooths <- te_smoothers(smooths, length(vars))
  # A TENSOR PRODUCT'S PENALTY IS A SUM OVER THE MARGINS, and a factory
  # would have to say how the penalty it builds composes with that sum --
  # whether it replaces one margin's roughness, all of them, or the sum. The
  # question has an answer only once somebody asks it of a model, so the
  # argument is rejected here rather than answered by a default.
  if (any(vapply(smooths, function(sm) !is.null(sm@penalty), logical(1)))) {
    stop(paste0(
      "a margin's 'penalty' factory is not read by te(): a tensor",
      " product's penalty\n  is the sum of the marginal roughnesses, and a",
      " factory does not say how it\n  composes with that sum. Use s() for",
      " a penalty of your own over one covariate."), call. = FALSE)
  }
  # ANISOTROPIC is one smoothing parameter per margin, so that is how many
  # names there are to hold; isotropic is the single one of a quadratic
  # penalty. Either way the names are the penalty's own, which is what the
  # summary prints and what the fit keys its hyperparameters by.
  nms <- if (anisotropic) paste0("lambda", seq_along(vars)) else "lambda"
  sp <- .smooth_spec(vars, substitute(by), smooths, label,
                     sprintf("te(%s)",
                             paste(vapply(vars, deparse, character(1)),
                                   collapse = ",")),
                     hyper, nms, sparse = sparse, ids = id,
                     by_hyper = by_hyper)
  sp@spec$anisotropic <- anisotropic
  sp
}


#' One Smoother per Margin of a Tensor Product
#'
#' @description
#' Resolves the `smooths` argument of [te()] into one smoother per covariate,
#' and rejects a margin asking for what the product does not read from it.
#'
#' @details
#' The product reads a margin's basis and its roughness matrix, so `k`,
#' `degree`, `order`, `measure` and the interval all act. The constraint, the
#' null space and the reparametrization are the **product's**: the tensor
#' basis contains the constant whatever its margins do, so the block carries
#' the sum-to-zero constraint over the observed covariates rather than a
#' marginal one, and the marginal linear effects are not separated out as
#' [s()] separates its one.
#'
#' Those three are therefore rejected at a non-default value rather than
#' accepted and ignored, which would report a fit of a model the caller did
#' not ask for.
#'
#' @param smooths One smoother, or a list of one per covariate.
#' @param nv How many covariates.
#'
#' @return A list of `nv` smoothers.
#'
#' @keywords internal
te_smoothers <- function(smooths, nv) {
  if (S7::S7_inherits(smooths, basis7::smoother)) {
    smooths <- rep(list(smooths), nv)
  }
  if (!is.list(smooths) || length(smooths) != nv ||
      !all(vapply(smooths, function(s) S7::S7_inherits(s, basis7::smoother),
                  logical(1)))) {
    stop(sprintf(paste0(
      "'smooths' must be one basis7 smoother, or a list of %d, one per",
      " covariate."), nv), call. = FALSE)
  }
  for (j in seq_len(nv)) {
    s <- smooths[[j]]
    bad <- c(
      if (!is.null(s@constrain)) "constrain",
      if (!identical(s@null_space, "keep")) "null_space",
      if (!identical(s@reparam, "dr")) "reparam"
    )
    if (length(bad)) {
      stop(sprintf(paste0(
        "the smoother of margin %d sets '%s', which a tensor product does",
        " not read\n  from a margin: the constraint, the null space and the",
        " coordinates belong to\n  the product, which is centered over the",
        " observed covariates and is not\n  reparametrized. What a margin",
        " supplies is its basis and its roughness\n  matrix: k, degree,",
        " order, measure, lower and upper."
      ), j, bad[[1L]]), call. = FALSE)
    }
  }
  smooths
}

.smooth_spec <- function(vars, by, smoothers, label,
                         default_label, hyper = NULL, names = "lambda",
                         sparse = NULL, ids = NULL,
                         by_hyper = c("shared", "level")) {
  nv <- length(vars)
  by_hyper <- match.arg(by_hyper)
  if (!is.list(smoothers) || length(smoothers) != nv ||
      !all(vapply(smoothers, function(s) S7::S7_inherits(s, basis7::smoother),
                  logical(1)))) {
    stop(sprintf(paste0(
      "the smoother must be a basis7 smoother%s. See ?basis7::bspline_smooth."),
      if (nv > 1L) sprintf(", or a list of %d", nv) else ""),
      call. = FALSE)
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
  # WHICH HYPERPARAMETERS THERE ARE is not always known here. With the
  # roughness matrix and one smoothing parameter shared across the levels the
  # names are `names`, and a wrong one is reported where it is written. A
  # penalty factory names its own, at a coefficient count only the data
  # settle, and one smoothing parameter per level needs the levels; both are
  # checked at the build by .entry_hyper(), the first point at which the
  # penalties exist. The strict check is kept where it can be made, so the
  # common call is unaffected.
  deferred <- identical(by_hyper, "level") ||
    any(vapply(smoothers, function(sm) !is.null(sm@penalty), logical(1)))
  SmoothTerm(label = label, vars = vars, by = by, sparse = sparse,
             spec = list(smoothers = smoothers, by_hyper = by_hyper),
             hyper = if (deferred) as_hyper(hyper, label)
                     else smooth_hyper(hyper, names, label),
             ids = check_ids(ids, if (deferred) NULL else names, label),
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
#' b <- term_build(s(x, basis7::bspline_smooth(k = 8)), dd)
#' names(b@blueprint)
#' b@blueprint$core$kind
#'
#' # And for a tensor product: one column fewer than 4 x 4.
#' bt <- term_build(te(x, z, smooths = basis7::bspline_smooth(k = 4)), dd)
#' term_npar(bt)
#'
#' @keywords internal
S7::method(term_build, SmoothTerm) <- function(term, data, ...) {
  xs <- .smooth_x(term@vars, data)
  sp <- term@spec
  nv <- length(xs)

  sms <- sp$smoothers
  marg <- lapply(seq_len(nv), function(j) {
    basis7::smoother_basis(sms[[j]], xs[[j]])
  })

  if (nv == 1L) {
    # THE WHOLE CONSTRUCTION IS THE SMOOTHER'S: the basis, the roughness
    # penalty, the null space and the coordinates. What is left here is
    # what belongs to a term in a formula -- the `by`, the label, the
    # storage and the penalty object.
    out <- basis7::smoother_build(sms[[1L]], xs[[1L]])
    Z <- out$X
    nm <- out$names
    P <- out$S
    # HOW MANY LEADING COLUMNS THE ROUGHNESS LEAVES ALONE, read off the
    # construction rather than derived from the arguments: it is one at
    # order 2, two at order 3, and none for a periodic basis, whose null
    # space is the constant the model's intercept already carries.
    unpen <- out$unpenalized
    core <- list(kind = "smoother", smoother = sms[[1L]],
                 blueprint = out$blueprint)
  } else {
    tb <- basis7::tensor_basis(marg)
    xm <- do.call(cbind, xs)
    # the marginal roughness penalties carried into the product: the second
    # derivative Gram of each margin, the identity in the others
    dims <- vapply(marg, function(b) b@dimension, integer(1))
    comps <- lapply(seq_len(nv), function(j) {
      # the margin's own roughness matrix, at the order and the measure its
      # smoother carries; at the defaults this is basis_gram(order = 2)
      Pj <- basis7::smoother_gram(sms[[j]], marg[[j]], xs[[j]])
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
    # the tensor block is centered over the observed covariates, so nothing
    # of it is left free: the constant, which is what the marginal null
    # spaces have in common, is the direction the constraint removed
    unpen <- 0L
    core <- list(kind = "tensor", basis = tb)
  }

  # THE PENALTY OF ONE LEVEL, kept before the `by` expansion widens it. A
  # per-level reading needs it as it stands here, one copy rather than m.
  P_block <- P
  kcol <- ncol(Z)

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

  lev <- if (!is.null(by) && identical(by$kind, "factor")) by$levels else NULL
  # THE FACTORY IS THE SMOOTHER'S, and only a one-covariate smooth has one
  # smoother to read it from. A tensor product's penalty is a sum over the
  # margins, and a factory would have to say how it composes with that sum,
  # so te() rejects the argument at construction rather than choosing.
  fac <- if (nv == 1L) sms[[1L]]@penalty else NULL
  parts <- .smooth_parts(P_block, kcol, unpen,
                         length(term@spec$smoothers) == 1L, lev, fac,
                         identical(term@spec$by_hyper, "level"), term@label)
  if (!is.null(parts)) {
    # the names are checked HERE, the first point at which these penalties
    # exist: a factory names its own hyperparameters at a coefficient count
    # the data settle, and one per level needs the levels
    parts <- .entry_hyper(parts, term@hyper, term@ids, term@label,
                          "this smooth's penalty")
    term@hyper <- do.call(c, c(list(list()), lapply(parts, function(en) {
      if (!length(en$fixed) || !nzchar(en$name)) return(en$fixed)
      stats::setNames(en$fixed, paste0(names(en$fixed), ".", en$name))
    })))
  }

  term@blueprint <- list(core = core, marg = marg, spec = sp,
                         vars = term@vars, by = term@by,
                         by_levels = if (!is.null(by) &&
                                         identical(by$kind, "factor"))
                           by$levels else NULL,
                         sparse = sp,
                         penalties = parts,
                         nblock = ncol(Z))
  # ONE PENALTY OVER THE WHOLE BLOCK is the reading every default fit runs,
  # and it is left exactly as it shipped. Where the term declares entries of
  # its own there is no such penalty, and the property says so rather than
  # carrying one of them: term_penalties() is the general question, and two
  # answers that could disagree are worse than one.
  term@penalty <- if (!is.null(parts)) NULL
                  else if (is.list(P)) penalties7::additive_penalty(P)
                  else penalties7::quadratic_penalty(P, blocks = by_blocks)
  term
}


# THE PENALTIES A BUILT SMOOTH DECLARES.
#
# NULL means the one the term carries over its whole block, which is the base
# reading of term_penalties() and what every fit written before this ran; a
# list is one entry per level, or one entry over a subset of the columns.
#
# The subset is what a factory needs and a matrix does not. A roughness matrix
# carries a zero row for a direction it does not see, so it may be handed the
# whole block and leaves the free columns alone by its own arithmetic. A
# penalty from a factory has no such row -- a lasso shrinks every coordinate
# it is given -- so it is told which coordinates are its own, and the
# unpenalized directions are not among them.
.smooth_parts <- function(P_block, kcol, unpen, one_var, lev, factory,
                          per_level, label) {
  m <- if (is.null(lev)) 1L else length(lev)
  if (per_level && is.null(lev)) {
    stop(sprintf(paste0(
      "'%s' asks for by_hyper = \"level\" and has no factor 'by'. One",
      " smoothing\n  parameter per level needs levels to have one each;",
      " give a factor 'by', or\n  leave by_hyper at \"shared\"."), label),
      call. = FALSE)
  }
  if (is.null(factory) && !per_level) return(NULL)

  whole <- seq_len(kcol)
  pen <- if (unpen > 0L) unpen + seq_len(kcol - unpen) else whole
  at <- function(j, within) (j - 1L) * kcol + within
  quad <- function(Pk) {
    if (is.list(Pk)) penalties7::additive_penalty(Pk)
    else penalties7::quadratic_penalty(Pk)
  }
  if (!is.null(factory) && !one_var) {
    stop("a penalty factory is not available on a tensor product.",
         call. = FALSE)
  }

  # the fields an entry carries beyond its own three, filled in by
  # .entry_hyper() from what the caller held
  entry <- function(name, index, penalty) {
    list(name = name, index = index, penalty = penalty,
         fixed = list(), n_values = list(), values = list(),
         min_ratio = numeric(0), ids = character(0))
  }

  if (!per_level) {
    # ONE penalty over every level's penalized coordinates. Repeating a
    # separable penalty blockwise is the identity operation, so for the
    # families a factory is reached for -- the lasso, the elastic net, SCAD,
    # MCP, a heavy-tailed prior -- this is the same penalty a per-level
    # reading would build, under one smoothing parameter instead of m.
    idx <- unlist(lapply(seq_len(m), function(j) at(j, pen)))
    return(list(entry("", idx, .penalty_factory(factory)(length(idx)))))
  }
  lapply(seq_len(m), function(j) {
    if (is.null(factory)) {
      entry(lev[j], at(j, whole), quad(P_block))
    } else {
      entry(lev[j], at(j, pen), .penalty_factory(factory)(length(pen)))
    }
  })
}


#' @title The Penalties of a Smooth Term
#' @name term_penalties.SmoothTerm
#' @description
#' One entry over the whole block where a single roughness matrix covers it,
#' which is the default; one per level of a factor `by` under
#' `by_hyper = "level"`; and one over the penalized coordinates alone where
#' the smoother carries a penalty factory.
#' @details
#' The default answer is the base reading of [term_penalties()], returned by
#' calling that method rather than by repeating it, so the two cannot drift.
#' @param term A built [SmoothTerm].
#' @param ... Unused.
#' @return A list of entries, as [term_penalties()] documents.
#' @keywords internal
S7::method(term_penalties, SmoothTerm) <- function(term, ...) {
  ent <- term@blueprint$penalties
  if (is.null(ent)) return(S7::method(term_penalties, model_term)(term, ...))
  ent
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
#' b <- term_build(s(x, basis7::bspline_smooth(k = 8)), dd)
#' X <- term_matrix(b)
#'
#' # Reapplying is exact; rebuilding on the same rows is a different basis.
#' max(abs(term_predict(b, dd[1:10, ]) - X[1:10, ]))
#' max(abs(term_matrix(term_build(s(x, basis7::bspline_smooth(k = 8)), dd[1:10, ])) - X[1:10, ]))
#'
#' # A factor `by` keeps every level's columns at a subset that has two.
#' bf <- term_build(s(x, basis7::bspline_smooth(k = 5), by = g), dd)
#' nd <- droplevels(dd[dd$g %in% c("a", "b"), ])
#' c(levels_here = nlevels(nd$g), cols = ncol(term_predict(bf, nd)))
#'
#' @keywords internal
S7::method(term_predict, SmoothTerm) <- function(term, newdata, ...) {
  .assert_built(term)
  bp <- term@blueprint
  xs <- .smooth_x(bp$vars, newdata)
  core <- bp$core

  if (identical(core$kind, "smoother")) {
    Z <- basis7::smoother_apply(core$smoother, core$blueprint, xs[[1L]])
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
