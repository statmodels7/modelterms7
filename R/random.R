#' @include term_classes.R generics.R
NULL

#' @title S7 Class for Grouped Random-Effect Terms
#' @name RandomTerm
#'
#' @description
#' `RandomTerm` is the subclass of [additive_term()] that holds grouped
#' coefficients with a distribution attached to them: the within-group design
#' interacted with the group indicators, one coefficient per group and per
#' within-group column, and a penalty that carries that distribution's negative
#' log-density. [random()] constructs it.
#'
#' @details
#' # The properties of its own
#'
#' `formula` is the bar formula as given, `~ 1 | g` or `~ x | g`, kept with its
#' environment. `correlated` indicates whether the default Gaussian lets the
#' within-group effects depend on each other. It is read only where `distrib`
#' is `NULL`, since a distribution given in `distrib` determines the same
#' thing.
#'
#' `tag` is the covariance label written between two bars, `~ 1 | u | g`, or
#' `NA_character_` where the formula carries one bar. It is reported by
#' [term_tag()]. The `tag` and `correlated` differ in scope: `correlated`
#' concerns the columns of this term, and `tag` the columns of other terms.
#'
#' `distrib` is the effects' distribution as supplied, or `NULL` for the
#' default. The build turns it into a \pkg{penalties7} penalty, read through
#' [term_penalty()] or [term_penalties()], so the hyperparameter names and
#' their bounds all come from the distribution.
#'
#' # The block is sparse by construction
#'
#' A row belongs to one group, so the block has a density of \eqn{1/m} and is
#' always a `dgCMatrix`. [random()] takes no `sparse` argument for that reason,
#' and passing one is an error.
#'
#' @inheritParams additive_term
#' @inheritParams model_term
#' @param formula The bar formula, `~ 1 | g` or `~ x | g`, with the
#'   within-group design on the left and the grouping variable on the right.
#' @param correlated A single logical: whether the default Gaussian lets the
#'   within-group effects correlate. Read only when `distrib` is `NULL`.
#' @param distrib The effects' distribution, a \pkg{distributions7} object or
#'   a list of them with one per within-group column, or `NULL` for the
#'   default Gaussian.
#' @param tag The covariance label from the middle of two bars, or
#'   `NA_character_`.
#'
#' @return An S7 object of class `RandomTerm`, inheriting from
#'   [additive_term()] and [model_term()], with the four properties above in
#'   addition to those of the two parent classes.
#'
#' @seealso [random()], the constructor; [term_penalties()] for the entries
#'   that a built one declares; [edf()] for the effective degrees of freedom of
#'   a fitted random effect.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(9), g = factor(rep(c("a", "b", "c"), 3)))
#'
#' tm <- random(~ x | g)
#' S7::S7_inherits(tm, RandomTerm)
#' tm@formula
#' tm@correlated
#'
#' # The block is one diagonal block per group and is always sparse.
#' b <- term_build(tm, dd)
#' class(term_matrix(b))
#' term_coef_names(b)
#'
#' # The hyperparameters are the effects' distribution's own.
#' term_penalty(b)@params
#'
#' @export
RandomTerm <- S7::new_class(
  name = "RandomTerm",
  parent = additive_term,
  properties = list(
    formula = S7::class_any,
    correlated = S7::class_logical,
    distrib = S7::class_any,
    tag = S7::class_character
  )
)

#' Grouped Random-Effect Term
#'
#' @description
#' Random intercepts and slopes for a grouping factor:
#' `random(~ 1 | g)` builds one coefficient per level of `g`,
#' and `random(~ x | g)` one intercept and one slope per level. The
#' distribution of the effects is attached as the penalty on those
#' coefficients, which is how a random effect enters a penalized likelihood.
#'
#' @details
#' The left side of the bar is an ordinary one-sided formula for the
#' within-group design, with the usual intercept convention:
#' `~ x | g` carries an intercept and a slope per group and
#' `~ 0 + x | g` carries only the slope. The block interacts that design
#' with the group indicators, ordered group by group, so the coefficients
#' of one group are adjacent.
#'
#' The formula and the distribution of the effects determine the term. The
#' parametrization in which the hyperparameters are estimated, their names,
#' their number and the points where the log-density has a kink are all
#' properties of that distribution, read from it at build time.
#'
#' @section The covariance label, and the second bar:
#' A formula may carry two bars, `~ 1 + x | u | g`, following \pkg{brms}. The
#' last position is the grouping variable and the middle one is a label: terms
#' carrying the same label and the same grouping describe effects that are
#' correlated with each other, so an intercept in one equation and a slope in
#' another can share one covariance block. R nests bars to the left, so the
#' parser itself yields that reading.
#'
#' The label is a name and not data. It is read as a symbol and never
#' evaluated, and a column of the data that carries the same name is not read,
#' so `~ 1 | u | g` builds the same block as `~ 1 | g`, whatever `u` is.
#'
#' A labelled term declares no penalty of its own. The covariance block
#' spans columns to which the term has no access, possibly in another
#' equation, so the prior over it belongs to the class and is built by the
#' fitting layer.
#' [term_penalties()] returns an empty list and [term_penalty()] is `NULL`.
#' Consequently `hyper` and `id` are rejected when the term is built, because
#' the hyperparameters are those of the class, and a `distrib` that is given on
#' a labelled term is the class's joint prior: it must be multivariate, and
#' its dimension is checked where the class is known.
#'
#' The two arguments have different scopes. `correlated` indicates whether the
#' columns of this term may depend on each other, and the label indicates
#' that its columns depend on those of other terms. The combination of
#' `correlated = FALSE` with a label is contradictory and is rejected when the
#' term is built.
#'
#' @section The distribution of the effects:
#' `distrib` is `NULL`, a \pkg{distributions7} object, or a list of
#' them with one per within-group column.
#'
#' A multivariate distribution of the within-group dimension lets the effects
#' of one group depend on each other, and its matrix parameter carries the
#' dependence: `mvgaussian2_distrib(2, ar1(2))` is a prior whose
#' precision is autoregressive, `mvstudent_t1_distrib(2)` a heavy-tailed
#' one. Correlated effects require a family that carries a matrix parameter,
#' that is, a location block as long as the dimension together with a
#' covariance, precision or scale matrix, and that defines the mixed
#' response-parameter derivative block read by a marginal criterion. The term
#' reads these properties from the family, so a family added later is covered
#' without an edit here.
#'
#' A univariate distribution makes the effects independent, and the penalty is
#' the product of the densities. With more than one within-group column it is
#' a template: one copy per column, each with its own hyperparameters, since
#' an intercept and a slope are quantities of different units and a shared
#' scale would price them against each other. A list of distributions gives
#' one per column explicitly, when the columns need different priors.
#'
#' The default is Gaussian: `gaussian1_distrib` at one column, so the
#' hyperparameter is the standard deviation of the effects; the multivariate
#' Gaussian on an unstructured covariance when there are several columns and
#' `correlated = TRUE`; and the template of the first when
#' `correlated = FALSE`, with one standard deviation per column. The default
#' has its location parameters held at zero.
#'
#' A distribution supplied in `distrib` must have its location parameters held
#' with [distributions7::fixed()]. A free location in the effects is
#' confounded with the intercept of the equation in which the term sits, so a
#' free location is rejected at build time with a message naming it. The value
#' at which it is held is usually zero and is not checked, since the model is
#' identified whatever the value is. Where the prior is a transformation of
#' another family, the parameter is the mean on the original scale, so
#' `fixed(transformation(gamma2_distrib(), log_transform()), mu = 1)` is a
#' log-gamma prior with shape \eqn{a = 1/\sigma^2}, where \eqn{\sigma^2} is the
#' parameter `sigma2` of the gamma. The mean of that prior,
#' \eqn{\psi(a) - \log a}, lies between \eqn{-\sigma^2} and
#' \eqn{-\sigma^2/2} and tends to zero as \eqn{\sigma^2} tends to zero.
#'
#' A distribution used as a penalty gives joint-mode (penalized likelihood)
#' estimation of the effects. With Gaussian effects and a Gaussian response
#' the Laplace approximation behind a marginal criterion is exact, so the
#' variance component that the criterion returns is the marginal estimate.
#' With any other prior the Laplace approximation is not exact, and
#' \pkg{modelterms7} does not compute the marginal likelihood.
#'
#' @section The hyperparameters:
#' They are the distribution's own free parameters, and every one of them is
#' estimated unless it is held. There are two ways to hold one, and the fit is
#' the same either way, but the reported output differs. Holding it inside
#' the distribution, `fixed(pseudohuber_distrib(), mu = 0, nu = 2)`, removes
#' it: it becomes a constant of the prior and is not among the
#' model's hyperparameters. Naming it in `hyper` keeps it, reported as held at
#' the value given, as a penalized term's own hyperparameter argument does.
#'
#' The hyperparameters of a prior with a smooth log-density are estimated by
#' the outer criterion. A prior whose log-density has a kink, a Laplace or an
#' elastic net, has none that a marginal criterion can reach, and its
#' hyperparameter is chosen by a path over the sparse criterion of
#' [statmodels7::statmod()], [statmodels7::bic()] by default.
#'
#' An estimated hyperparameter, including a shape parameter, is reported with
#' a standard error and an interval. The standard error is absent where
#' the criterion at the point at which the run ended has a curvature of the
#' wrong sign, since such a point does not yield an interval. A shape
#' parameter that moves toward a limit is the common case.
#'
#' The parametrization of a family that is used as a prior matters. The
#' centered skew normal ([distributions7::skewnormal2_distrib()]) carries the
#' skewness itself, and its map to the direct parametrization is not twice
#' differentiable at zero skewness. The derivatives of the log-density with
#' respect to the skewness have a finite first-order limit there, and the
#' second-order ones grow like \eqn{\gamma_1^{-2/3}}. A marginal criterion
#' reads the second-order ones, and the symmetric bounds place the starting
#' value at that point, so the direct parametrization
#' ([distributions7::skewnormal1_distrib()]) is the one to use as a prior; its
#' derivatives at \eqn{\alpha = 0} are finite.
#'
#' How well a shape parameter is estimated depends on the number of groups,
#' since the parameter is read from that many latent values, which the prior
#' shrinks. With few groups the estimate of a shape parameter can move to the
#' Gaussian limit, and then only the scale is estimated. With more groups the
#' estimate stays finite, but the criterion is flat over the small values, so
#' the estimate locates the shape only to its order of magnitude. A shape
#' parameter can be held with [distributions7::fixed()] when the number of
#' groups is small, and an estimate from a moderate number of groups describes
#' the weight of the tail without fixing its value. The \eqn{\nu} of a
#' pseudo-Huber prior, the point at which the loss stops being quadratic, is
#' subject to the same limitation.
#'
#' Prediction maps new data onto the levels seen at build time. A level that
#' was not present at build time, or a missing grouping value, is rejected by
#' [term_predict()] unless `unseen = "zero"` is given, in which case those rows
#' have no effect of their own.
#'
#' A random effect is not standardized, and `random()` has no `standardize`
#' argument; passing one is an error. Its columns are grouping indicators and
#' its hyperparameter is a variance component with a meaning of its own.
#' Dividing each coefficient by the spread of its indicator would weight the
#' effects by the sizes of the groups, which changes the model itself.
#'
#' @section The block and its penalty:
#' With \eqn{m} levels and a within-group design \eqn{Z_i} of \eqn{d}
#' columns, the block is the interaction of that design with the group
#' indicators, ordered group by group. When the rows of the data are ordered by
#' group, the block and the coefficient vector are
#'
#' \deqn{Z = \operatorname{diag}(Z_1, \dots, Z_m),
#'   \qquad b = (b_1', \dots, b_m')',}
#'
#' and otherwise the rows of this matrix are permuted to follow the data. The
#' \eqn{d} coefficients of one group occupy adjacent positions and the
#' penalty reads them one block at a time,
#'
#' \deqn{\rho(b; \theta) = -\sum_{i=1}^m \log f(b_i; \theta),}
#'
#' for the effects' density \eqn{f}. Under the default Gaussian that is
#'
#' \deqn{\rho(b; \theta) = \tfrac{1}{2}\sum_i b_i'\Sigma(\theta)^{-1}b_i
#'   + \tfrac{m}{2}\log\lvert \Sigma(\theta)\rvert
#'   + \tfrac{md}{2}\log(2\pi),}
#'
#' and minimizing the penalized least squares in \eqn{(\beta, b)} is the
#' mixed-model equation at the variance ratio that \eqn{\Sigma} encodes, so the
#' minimizer is the best linear unbiased predictor. At \eqn{d = 1} and
#' \eqn{\Sigma = \sigma_b^2} it is the ridge, up to the constant that makes
#' \eqn{\sigma_b} estimable.
#'
#' @section Smooth random effects:
#' Where the left side of the bar is an [s()] or a [te()] call,
#' `random(~ s(x) | g)`, each level of `g` has a random curve of its own (a
#' random surface for [te()]), and the term is a [RandomSmoothTerm()]. When
#' the smooth's columns are split into the null space of its penalty
#' \eqn{N(x)} and the penalized coordinates \eqn{Z(x)}, group \eqn{i}
#' contributes
#'
#' \deqn{g_i(x) = a_{0i} + N(x)a_i + Z(x)\delta_i, \qquad
#'   (a_{0i}, a_i) \sim \mathrm{N}(0, \Sigma), \qquad
#'   \delta_i \sim \mathrm{N}(0, \sigma^2 I),}
#'
#' independently over the groups. The null space is whatever the smoother's
#' operator leaves unpenalized, less the constant: a straight line for a
#' second-derivative penalty, a sine and a cosine for the harmonic operator,
#' nothing for a cyclic smoother. The argument `null_space` of the smoother
#' determines what the null space contributes: under `"keep"` its columns carry
#' \eqn{\Sigma} with the random level, under `"shrink"` they join the penalized
#' coordinates, and under `"drop"` they are absent. `~ 0 + s(x) | g` removes the
#' random level.
#'
#' \eqn{\Sigma} is the prior that `distrib` and `correlated` describe for the
#' within-group columns of an ordinary random effect, here over the random
#' level and the null space. With an unstructured \eqn{\Sigma} the model does
#' not depend on the origin and the scale of `x`, nor on the phase of a
#' periodic one. The penalized coordinates are whitened, so that their prior
#' is independent with one standard deviation, reported as `sigma` of the
#' entry `smooth`. A smoother with several smoothing parameters (an adaptive
#' smoother, an anisotropic tensor product) keeps them: the precision is
#' \eqn{\sum_k \lambda_k S_k} in every group, and the summary reports
#' \eqn{\lambda_k^{-1/2}}. All the groups share these hyperparameters, so the
#' curves have one roughness. The population curve is a separate [s()] term
#' in the same formula, as in `y ~ s(x) + random(~ s(x) | g)`.
#'
#' The smooth is written without `by`, `sparse`, `hyper` or `id`. The
#' hyperparameters belong to the random effect, and `hyper` and `id` of this
#' function name them (`sigma.smooth`, `lambda1.smooth`, ...). A covariance
#' label between two bars is rejected.
#'
#' @param formula A bar formula, `~ lhs | g`, where `g`
#'   evaluates to the grouping variable in the data, or `~ lhs | tag | g`
#'   with a covariance label in the middle. `g` must not be continuous: a
#'   double carrying a value that is not a whole number is rejected, since its
#'   levels would be its formatted values.
#' @param distrib The distribution of the effects: `NULL` (the default
#'   Gaussian), a \pkg{distributions7} object, or a list of them with one per
#'   within-group column.
#' @param correlated Logical; whether the default Gaussian lets the
#'   within-group effects correlate. Giving it together with `distrib` is an
#'   error, since `distrib` determines the same thing.
#' @param label A single non-empty string prefixed to the coefficient
#'   names.
#' @param hyper The hyperparameters of the effects' distribution to hold, as
#'   a named vector or list; those not named are estimated. For a distribution
#'   whose log-density has a kink, several values for one name are the grid
#'   that the path visits, and for any other distribution they are rejected.
#'   The names are the distribution's own parameters, and the within-group
#'   column is appended where there is one copy per column. A name that the
#'   penalty does not carry is rejected when the term is built, which is where
#'   the penalty first exists.
#' @param id Labels under which the effects' hyperparameters are shared with
#'   the hyperparameters of other terms carrying the same labels, so that two
#'   random effects may be given one variance component. The names are
#'   `hyper`'s, checked at the same point and for the same reason; a single
#'   unnamed string is accepted where the effects' distribution carries one
#'   hyperparameter. `NULL`, the default, shares nothing. See [term_ids()].
#' @return An object of class [RandomTerm()] (a specification;
#'   see [term_build()]).
#'
#' @examples
#' dd <- data.frame(y = rnorm(9), x = rnorm(9),
#'                  g = factor(rep(c("a", "b", "c"), 3)))
#'
#' # The term has one variance component, which is reported as a standard
#' # deviation.
#' term_penalty(term_build(random(~ 1 | g), dd))@params
#'
#' # Intercepts and slopes are correlated, and the covariance of the effects
#' # is the hyperparameter.
#' built <- term_build(random(~ x | g), dd)
#' term_coef_names(built)
#' term_penalty(built)@params
#'
#' # The effects are independent, and each column has its own standard
#' # deviation.
#' vapply(term_penalties(term_build(random(~ x | g, correlated = FALSE), dd)),
#'        function(e) e$name, "")
#'
#' # A heavy-tailed prior is held at four degrees of freedom.
#' t4 <- distributions7::fixed(distributions7::student_t1_distrib(),
#'                             mu = 0, nu = 4)
#' term_penalty(term_build(random(~ 1 | g, distrib = t4), dd))@params
#'
#' # A covariance label is read from the middle of two bars and is reported by
#' # term_tag(). It is a name: the block is the one bar formula's.
#' term_tag(random(~ 1 | u | g))
#' identical(term_coef_names(term_build(random(~ 1 | u | g), dd)),
#'           term_coef_names(term_build(random(~ 1 | g), dd)))
#'
#' # A continuous grouping variable is rejected rather than turned into one
#' # level per row.
#' try(term_build(random(~ 1 | x), dd))
#'
#' # A smooth random effect has one random curve per level of g.
#' sd <- data.frame(x = runif(60), g = factor(rep(c("a", "b", "c"), 20)))
#' rs <- term_build(random(~ s(x, smoother = basis7::bspline_smooth(k = 5)) |
#'                           g), sd)
#' term_group(rs)$names
#'
#'
#' # The model is fitted to data simulated from a known truth, so the
#' # estimates below can be read against it.
#' if (requireNamespace("statmodels7", quietly = TRUE)) {
#'   set.seed(6)
#'   fd <- data.frame(gr = factor(rep(1:20, each = 15)), x = rnorm(300))
#'   bb <- rnorm(20, sd = 0.8)
#'   fd$y <- 1 + 0.5 * fd$x + bb[fd$gr] + rnorm(300, sd = 0.5)
#'   cf <- coef(statmodels7::statmod(y ~ x + random(~1 | gr),
#'                                   distributions7::gaussian1_distrib(), fd))$mu
#'   # The true slope is 0.5, and the shrunken effects track the ones drawn.
#'   round(c(slope = cf[["x"]],
#'           cor = cor(cf[grep("^random", names(cf))], bb)), 3)
#' }
#' @references
#' Laird, N. M. and Ware, J. H. (1982). Random-effects models for
#' longitudinal data. *Biometrics* 38, 963-974.
#'
#' @seealso [s()], [te()], [nl()]
#' @export
random <- function(formula, distrib = NULL, correlated = TRUE,
                   label = "random", hyper = NULL, id = NULL) {
  parts <- .random_parts(formula)
  if (!is.logical(correlated) || length(correlated) != 1L ||
      is.na(correlated)) {
    stop("'correlated' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(distrib) && !missing(correlated)) {
    stop(paste0(
      "'correlated' and 'distrib' both specify how the effects are\n",
      "  distributed. Say it once: a multivariate distribution lets them\n",
      "  correlate and a univariate one does not."), call. = FALSE)
  }
  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !nzchar(label)) {
    stop("'label' must be a single non-empty character string.",
         call. = FALSE)
  }
  sm <- .random_smooth_lhs(parts$within, environment(formula))
  if (!is.null(sm)) {
    if (!is.na(parts$tag)) {
      stop(paste0(
        "a smooth random effect carries no covariance label: its null space\n",
        "  and its penalized coordinates have priors of their own, which a\n",
        "  block shared with other terms would replace."), call. = FALSE)
    }
    return(RandomSmoothTerm(
      label = label, formula = formula, correlated = correlated,
      distrib = distrib, tag = parts$tag,
      hyper = as_hyper(hyper, label),
      ids = check_ids(id, NULL, label),
      X = NULL, coef_names = character(0),
      blueprint = list(), penalty = NULL,
      smooth = sm$term, intercept = sm$intercept))
  }
  RandomTerm(label = label, formula = formula, correlated = correlated,
             distrib = distrib, tag = parts$tag,
             hyper = as_hyper(hyper, label),
             ids = check_ids(id, NULL, label),
             X = NULL, coef_names = character(0),
             blueprint = list(), penalty = NULL)
}

#' The Three Pieces of a Bar Formula
#'
#' @description
#' Splits `~ lhs | g` or `~ lhs | tag | g` into the within-group design, the
#' grouping expression and the covariance label, and rejects the shapes that
#' are not one of the two.
#'
#' @details
#' R nests the bars to the left, so in `~ 1 + x | u | id` the outer bar has
#' `id` on its right and `1 + x | u` on its left. The last position is
#' therefore the grouping variable and the middle one the label, which is
#' \pkg{brms}'s convention and is the reading that the parser itself yields.
#'
#' The label is taken with `as.character()` on the symbol and never
#' evaluated, so a column of the data carrying that name is not read. It is
#' the one element of the formula that does not refer to the data.
#'
#' Two bars are the most that carry a meaning, so a third is rejected.
#'
#' The constructor and [term_build()] both call this function, so that the
#' split is defined in one place.
#'
#' @param formula The bar formula as given.
#'
#' @return A list with `within`, the left side of the bar (a call, a symbol or
#'   the constant `1`), `group`, a symbol or a call, and `tag`, a single string
#'   or `NA_character_`.
#'
#' @seealso [random()] and [term_build.RandomTerm()], its two callers;
#'   [term_tag()] for what reports the label afterwards.
#'
#' @keywords internal
.random_parts <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 2L) {
    stop("'formula' must be a one-sided bar formula, e.g. ~ 1 | g.",
         call. = FALSE)
  }
  e <- formula[[2L]]
  if (!is.call(e) || !identical(e[[1L]], as.name("|"))) {
    stop("'formula' must contain a grouping bar, e.g. ~ 1 | g.",
         call. = FALSE)
  }
  within <- e[[2L]]
  tag <- NA_character_
  if (is.call(within) && identical(within[[1L]], as.name("|"))) {
    tg <- within[[3L]]
    within <- within[[2L]]
    if (is.call(within) && identical(within[[1L]], as.name("|"))) {
      stop(paste0(
        "'formula' has three bars. Two bars are read as ~ lhs | tag | g,\n",
        "  the within-group design, the label of the covariance block and the\n",
        "  grouping variable, and a third bar would have nothing to name."),
        call. = FALSE)
    }
    if (!is.name(tg)) {
      stop(sprintf(paste0(
        "the covariance label between the two bars must be a name, and '%s'\n",
        "  is not. In ~ lhs | tag | g the middle position labels the effects\n",
        "  that share a covariance block; it is not evaluated and does not\n",
        "  refer to the data."), deparse1(tg)), call. = FALSE)
    }
    tag <- as.character(tg)
  }
  list(within = within, group = e[[3L]], tag = tag)
}

#' What a Grouping Variable Has to Be
#'
#' @description
#' Rejects a grouping expression whose value is continuous, naming it and the
#' spelling to use when the values really are labels.
#'
#' @details
#' A grouping variable's values are labels, and `factor()` turns any vector
#' into one by formatting its values, so a continuous covariate in that
#' position would give one level per distinct value and a random effect with as
#' many coefficients as there are rows.
#'
#' Only a double carrying a value that is not a whole number is rejected.
#' Integer codes are the ordinary way to write a grouping variable and are
#' accepted, and a caller who means the levels of a non-integer vector writes
#' `factor()` around it, which arrives here as a factor and passes.
#'
#' The check also applies to the two-bar form: `~ 1 | u | id` and
#' `~ 1 | id | u` have the same syntactic shape and differ only in which
#' position holds which element, so the last position is required to be usable
#' as a grouping variable.
#'
#' @param v The value of the grouping expression.
#' @param expr The expression, for the message.
#'
#' @return `v`, unchanged; called for its error.
#'
#' @seealso [random()], whose `formula` argument documents the rule;
#'   [term_build.RandomTerm()], which reaches it through the grouping
#'   variable's evaluation.
#'
#' @keywords internal
.random_check_group <- function(v, expr) {
  if (is.double(v) && any(is.finite(v) & v != floor(v))) {
    stop(sprintf(paste0(
      "the grouping variable '%s' is continuous, and a grouping variable's\n",
      "  values are labels rather than measurements: its levels would be the\n",
      "  formatted numbers, one group per distinct value. Write factor(%s) if\n",
      "  they really are labels."), deparse1(expr), deparse1(expr)),
      call. = FALSE)
  }
  v
}

.random_group <- function(expr, data, levels = NULL, unseen = "error") {
  v <- eval(expr, data, baseenv())
  if (is.null(levels)) return(factor(.random_check_group(v, expr)))
  f <- factor(v, levels = levels)
  if (identical(unseen, "error") && anyNA(v)) {
    stop(sprintf(paste("the grouping value is missing at row %d; use",
                       "unseen = \"zero\" for a row of zeros."),
                 which(is.na(v))[1L]), call. = FALSE)
  }
  if (identical(unseen, "error") && any(is.na(f) & !is.na(v))) {
    bad <- unique(as.character(v)[is.na(f) & !is.na(v)])
    stop(sprintf("grouping level '%s' was not present at build time.",
                 bad[1L]), call. = FALSE)
  }
  f
}

# The within-group design and the indicators, interacted group by group so
# the coefficients of one group are adjacent -- the order the blockwise
# penalty reads, and the order I_m %x% Sigma assumes.
#
# The block is SPARSE by construction and is built as such rather than
# assembled dense and converted: a row belongs to one group, so exactly d of
# its m*d entries are non-zero and the density is 1/m whatever the data. The
# dense form was quadratic in the wrong place -- at n = 20000 and m = 1000 it
# is 152.6 MB against 0.23 MB, built in 1.76 s against 0.0011 s, and the
# crossprod every iteration takes is 12.77 s against 0.0006 s. The
# intermediate `outer(g, levels(g), ==)` was itself a dense n x m.
.random_block <- function(g, W) {
  m <- nlevels(g)
  d <- ncol(W)
  n <- nrow(W)
  gi <- as.integer(g)
  # a row whose group is NA -- a level the fit never saw, read as having no
  # effect of its own -- lands in no column and is a row of zeros
  ok <- which(!is.na(gi))
  # one entry per (row, within-column) pair: row i lands in the d columns of
  # its own group, and in no other
  Matrix::sparseMatrix(
    i = rep(ok, each = d),
    j = as.vector(t(outer((gi[ok] - 1L) * d, seq_len(d), `+`))),
    x = as.vector(t(W[ok, , drop = FALSE])),
    dims = c(n, m * d))
}

.random_names <- function(label, glevels, wnames) {
  if (length(wnames) == 1L && wnames == "(Intercept)") {
    return(paste(label, glevels, sep = "."))
  }
  as.character(t(outer(glevels, wnames,
                       function(a, b) paste(label, a, b, sep = "."))))
}

#' The Family Under Any Wrappers
#'
#' @description
#' Follows `parent_distrib` to the distribution that a wrapper wraps, so that
#' a property of the family is read from the family.
#'
#' @details
#' The property is detected with `S7::prop_names()` instead of testing the
#' class, so a wrapper written later is followed without an edit here.
#'
#' @param d A \pkg{distributions7} object.
#'
#' @return A \pkg{distributions7} object.
#'
#' @keywords internal
.random_family <- function(d) {
  while ("parent_distrib" %in% S7::prop_names(d)) d <- d@parent_distrib
  d
}

#' Whether a Multivariate Family Can Carry Correlated Effects
#'
#' @description
#' Returns `TRUE` for a family with a location block as long as its dimension
#' and a matrix parameter, both of which a centered prior on \eqn{R^d}
#' requires.
#'
#' @details
#' This is a property of the family, so a multivariate family added later is
#' covered without an edit here. It is read from `params_interpretation`, the
#' same declaration from which a data-based starting value is built. The
#' condition excludes the simplex-valued families,
#' whose mean coordinates are one fewer than the dimension and which carry no
#' matrix parameter at all.
#'
#' @param d A \pkg{distributions7} object, under any wrappers.
#'
#' @return A single logical.
#'
#' @keywords internal
.random_mv_ok <- function(d) {
  interp <- d@params_interpretation
  sum(interp %in% c("location", "mean")) == d@n_dim &&
    any(interp %in% c("covariance", "precision", "scale"))
}

#' Whether a Multivariate Family Defines Its Mixed Block
#'
#' @description
#' Returns `TRUE` when `distrib_cross_y` comes from the family itself. The
#' multivariate base class rejects that generic, so the function returns
#' `FALSE` for a family that has not overridden it.
#'
#' @details
#' A marginal criterion reads that block to estimate the covariance of the
#' effects, so a family without one is rejected when the term is built. The
#' class that owns the method is read from its
#' signature and compared by name, since `identical()` compares object
#' identity and fails when a package's code is re-evaluated instead of loaded.
#'
#' @param d A \pkg{distributions7} object, under any wrappers.
#'
#' @return A single logical.
#'
#' @keywords internal
.random_has_cross <- function(d) {
  m <- tryCatch(S7::method(distributions7::distrib_cross_y, S7::S7_class(d)),
                error = function(e) NULL)
  if (is.null(m)) return(FALSE)
  !identical(attr(attr(m, "signature")[[1L]], "name"), "multivariate_distrib")
}

#' What a Prior on the Effects Has to Satisfy
#'
#' @description
#' Checks a caller-supplied effects distribution where the term is built and
#' names what is wrong with it.
#'
#' @details
#' A free location is rejected. It is confounded with the intercept of
#' the equation in which the term sits, leaving a flat direction in the
#' objective.
#'
#' A location held at a value is identified whatever that value is, and it is
#' not policed: it shrinks the effects toward that value, which is a modeling
#' statement. Nor could the value be policed in general. Where the prior is a
#' transformation of another family, the parameter interpreted as its mean is
#' the mean on the original scale, and holding the mean of a gamma at one
#' centers its logarithm.
#'
#' @param d A \pkg{distributions7} object.
#' @param dim_needed The number of within-group columns.
#' @param what How to name the argument in a message.
#'
#' @return Invisibly `TRUE`; called for its errors.
#'
#' @keywords internal
.random_check_prior <- function(d, dim_needed, what) {
  if (!S7::S7_inherits(d, distributions7::distrib)) {
    stop(sprintf("%s must be a distributions7 object.", what), call. = FALSE)
  }
  base <- .random_family(d)
  mv <- identical(d@dimension, "multivariate")
  if (!mv && !identical(d@dimension, "univariate")) {
    stop(sprintf("%s is %s, and effects are distributed on the line or on R^d.",
                 what, d@dimension), call. = FALSE)
  }
  if (S7::S7_inherits(base, distributions7::discrete_distrib)) {
    stop(sprintf("%s is discrete, and a random effect is not.", what),
         call. = FALSE)
  }
  if (mv) {
    if (!.random_mv_ok(base)) {
      stop(sprintf(paste0(
        "%s ('%s') carries no matrix parameter over R^%d, so it cannot express\n",
        "  how the effects of one group depend on each other. Correlated\n",
        "  effects need a family with a location block and a covariance,\n",
        "  precision or scale matrix, such as mvgaussian1_distrib() or\n",
        "  mvstudent_t1_distrib()."),
        what, base@distrib_name, base@n_dim), call. = FALSE)
    }
    if (d@n_dim != dim_needed) {
      stop(sprintf(paste0(
        "%s is %d-variate and the bar's left side has %d columns; a\n",
        "  multivariate prior describes the effects of one group."),
        what, d@n_dim, dim_needed), call. = FALSE)
    }
    if (!.random_has_cross(base)) {
      stop(sprintf(paste0(
        "%s ('%s') has no mixed response-parameter block, which a marginal\n",
        "  criterion reads to estimate the covariance of the effects. It is\n",
        "  rejected here rather than at the criterion, where the message\n",
        "  would name a generic. The gaussian and the Student t have one."),
        what, base@distrib_name), call. = FALSE)
    }
  }
  interp <- base@params_interpretation
  locs <- names(interp)[interp %in% c("location", "mean")]
  free <- intersect(locs, d@params)
  if (length(free)) {
    stop(sprintf(paste0(
      "%s has a free location ('%s'), which is confounded with the\n",
      "  intercept of the equation in which the term sits. Hold it with\n",
      "  distributions7::fixed(): at zero for a family on the line, and at\n",
      "  whatever centers the effects where the prior is a transformation of\n",
      "  another family, because the parameter is then the mean on the original\n",
      "  scale."),
      what, free[1L]), call. = FALSE)
  }
  invisible(TRUE)
}

#' The Centered Gaussian Defaults
#'
#' @description
#' Returns `gaussian1_distrib` at one within-group column, so the
#' hyperparameter is the standard deviation of the effects, and the
#' multivariate Gaussian on an unstructured covariance for several correlated
#' ones.
#'
#' @details
#' The matrix parametrization is the family's, and the family built here is
#' `mvgaussian1_distrib()`, so the structure is the effects' covariance
#' and its hyperparameters are the log-Cholesky coordinates of that covariance:
#' the logarithms of the diagonal entries of its Cholesky factor, and the
#' off-diagonal entries of that factor. The precision reading is a different
#' prior, because the two differ in the sign of the log-determinant term, and
#' it is written by passing `mvgaussian2_distrib()` to `distrib` instead.
#'
#' @param d The number of within-group columns.
#' @param correlated Whether the effects may correlate.
#'
#' @return A \pkg{distributions7} object.
#'
#' @keywords internal
.random_default <- function(d, correlated) {
  if (d == 1L || !correlated) {
    return(distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0))
  }
  mv <- distributions7::mvgaussian1_distrib(
    d, sigma = parameters7::log_cholesky(d))
  do.call(distributions7::fixed,
          c(list(mv), stats::setNames(as.list(rep(0, d)),
                                      paste0("mu", seq_len(d)))))
}

#' Penalty Entries of a Random-Effect Term
#'
#' @description
#' Returns one entry over the whole block when the prior reads a group's
#' effects together, and one per within-group column when they are independent.
#'
#' @details
#' The coefficients are ordered group by group, so column \eqn{j} is the
#' stride \eqn{j, d+j, 2d+j, \dots}. Those positions are named as a subset of
#' the term's own parameters, never selected with a map, which
#' keeps a kinked prior's proximal operator available: a separable penalty
#' under a selection map is the generalized-lasso problem, which has none.
#'
#' @param prior The effects' distribution, or a list of one per column.
#' @param d The number of within-group columns.
#' @param m The number of groups.
#' @param wnames The within-group column names.
#'
#' @return A list of entries, as [term_penalties()] documents.
#'
#' @keywords internal
.random_entries <- function(prior, d, m, wnames) {
  if (!is.list(prior)) {
    if (identical(prior@dimension, "multivariate") || d == 1L) {
      pen <- penalties7::distrib_penalty(prior, n_coef = m * d)
      return(list(.penalty_entry("", seq_len(m * d), pen)))
    }
    prior <- rep(list(prior), d)
  }
  lapply(seq_len(d), function(j) {
    .penalty_entry(wnames[j], seq(j, by = d, length.out = m),
                   penalties7::distrib_penalty(prior[[j]], n_coef = m))
  })
}

#' @title Build a Random-Effect Term
#' @name term_build.RandomTerm
#'
#' @description
#' Builds the within-group design from the left of the bar, interacts it with
#' the group indicators, and attaches the effects' distribution as the penalty
#' on the resulting coefficients. The levels of the grouping variable are
#' recorded, so [term_predict()] maps new rows onto the same ones.
#'
#' @details
#' # The block
#'
#' With \eqn{m} levels and a within-group design \eqn{Z_i} of \eqn{d} columns,
#' the block is \eqn{\mathrm{diag}(Z_1, \dots, Z_m)} for data ordered by group,
#' and its columns are ordered group-major, so the \eqn{d} coefficients of one
#' group are adjacent. It is built as a
#' `dgCMatrix`: a row belongs to one group, so the density is \eqn{1/m}.
#'
#' The coefficient names are `label.level.column`, or `label.level` when the
#' within-group design consists only of the intercept. The term `random(~ x |
#' g)` over three levels gives `random.a.(Intercept)`, `random.a.x`,
#' `random.b.(Intercept)` and so on, which is the group-major order, and
#' `random(~ 1 | g)` gives `random.a`, `random.b` and `random.c`. When the rows
#' of the data are not ordered by group, the rows of the block follow the data
#' and its columns keep this order.
#'
#' A term with a single random effect (one level and one within-group column)
#' is rejected. One level with several within-group columns is accepted: the
#' effects are then a single draw of a vector, such as a latent process over
#' the levels of a within-group factor.
#'
#' # The penalty, and what the build checks
#'
#' Where `distrib` is `NULL` the default is chosen here: a centered
#' `gaussian1_distrib` at one column or under `correlated = FALSE`, and a
#' centered multivariate Gaussian on an unstructured covariance for several
#' correlated ones.
#'
#' Whatever the distribution, its location parameters must be held. A free
#' location is confounded with the intercept of the equation in which the term
#' sits, and the build rejects it with a message naming the parameter and the
#' remedy. A multivariate distribution must also match the within-group
#' dimension, and one that carries no matrix parameter cannot express
#' correlated effects at all.
#'
#' Any value named in `hyper` is checked here too, against the penalty's own
#' names, because this is the first point at which the penalty exists.
#'
#' A term carrying a covariance label is finished differently: its block and
#' coefficient names are built as above, and it declares no penalty.
#' `hyper`, `id`, `correlated = FALSE` and a `distrib` that is a list or is
#' univariate are rejected for such a term.
#'
#' @param term An unbuilt or built [RandomTerm()].
#' @param data A data frame carrying the grouping variable and the
#'   within-group covariates.
#' @param ... Unused.
#'
#' @return The term with `X` (a `dgCMatrix` of \eqn{md} columns),
#'   `coef_names`, `blueprint` and `penalty` filled.
#'
#' @seealso [random()], [term_predict.RandomTerm()], [term_penalties()].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(9), g = factor(rep(c("a", "b", "c"), 3)))
#'
#' # Group-major: the two coefficients of one level are adjacent.
#' b <- term_build(random(~ x | g), dd)
#' term_coef_names(b)
#' as.matrix(term_matrix(b))
#'
#' # A free location is rejected, naming the parameter.
#' try(term_build(random(~ 1 | g,
#'                       distrib = distributions7::gaussian1_distrib()), dd))
#'
#' @keywords internal
S7::method(term_build, RandomTerm) <- function(term, data, ...) {
  parts <- .random_parts(term@formula)
  g <- .random_group(parts$group, data)
  m <- nlevels(g)

  wf <- stats::as.formula(call("~", parts$within),
                          env = environment(term@formula))
  mf <- stats::model.frame(wf, data, na.action = stats::na.pass,
                           drop.unused.levels = FALSE)
  tt <- attr(mf, "terms")
  W <- stats::model.matrix(tt, mf)
  contr <- attr(W, "contrasts")
  d <- ncol(W)
  if (d < 1L) {
    stop("the left side of the bar yields no columns.", call. = FALSE)
  }
  # one group is a legitimate term when it carries several effects: a latent
  # process over the levels of a within-group factor is a single draw of a
  # d-dimensional vector
  if (m * d < 2L) {
    stop(paste0(
      "the term has a single random effect: the grouping variable has one\n",
      "  level and the left side of the bar one column."), call. = FALSE)
  }

  Z <- .random_block(g, W)
  cn <- .random_names(term@label, levels(g), colnames(W))
  colnames(Z) <- cn

  # A LABELLED term does not own its penalty. Its coefficients share a
  # covariance block with those of the other terms carrying the label, and
  # that block spans columns this term cannot see -- it may not even be in
  # this equation -- so the prior belongs to the class and is built by the
  # fitting layer. Declaring one here would be a second answer to the same
  # question, and the fit would have to know which to ignore.
  if (!is.na(.random_tag(term))) {
    return(.random_labelled(term, Z, cn, parts, g, tt, mf, contr, d,
                            colnames(W)))
  }

  prior <- term@distrib
  if (is.null(prior)) {
    prior <- .random_default(d, term@correlated)
  } else if (is.list(prior)) {
    if (length(prior) != d) {
      stop(sprintf(paste0(
        "'distrib' is a list of %d and the bar's left side has %d columns;\n",
        "  a list gives one distribution per column."),
        length(prior), d), call. = FALSE)
    }
    for (j in seq_len(d)) {
      .random_check_prior(prior[[j]], d, sprintf("'distrib[[%d]]'", j))
      if (identical(prior[[j]]@dimension, "multivariate")) {
        stop(sprintf(paste0(
          "'distrib[[%d]]' is multivariate. A list gives one distribution\n",
          "  per within-group column, so each is univariate; pass a single\n",
          "  multivariate distribution to correlate the columns instead."),
          j), call. = FALSE)
      }
    }
  } else {
    .random_check_prior(prior, d, "'distrib'")
  }

  # the hyperparameter names are checked HERE, the first point at which the
  # penalties this term builds exist: which there are depends on the
  # distribution of the effects and on how many within-group columns there are
  entries <- .entry_hyper(.random_entries(prior, d, m, colnames(W)),
                           term@hyper, term@ids, term@label)

  term@X <- Z
  term@coef_names <- cn
  term@blueprint <- list(
    gexpr = parts$group, glevels = levels(g),
    terms = stats::delete.response(tt),
    xlev = stats::.getXlevels(tt, mf),
    contrasts = contr,
    entries = entries, dim = d, wnames = colnames(W)
  )
  term@hyper <- do.call(c, c(list(list()), lapply(entries, function(en) {
    if (!length(en$fixed) || !nzchar(en$name)) return(en$fixed)
    stats::setNames(en$fixed, paste0(names(en$fixed), ".", en$name))
  })))
  # term_penalty() answers for the common case, one penalty over the whole
  # block, and is NULL where there is one per column -- which is what a
  # reader of a partial penalty already has to handle
  term@penalty <- if (length(entries) == 1L && !nzchar(entries[[1L]]$name)) {
    entries[[1L]]$penalty
  } else {
    NULL
  }
  term
}

#' @title The Penalties of a Random-Effect Term
#' @name term_penalties.RandomTerm
#' @description
#' Returns one entry over the whole block where the effects of a group are
#' read together, and one per within-group column where they are independent.
#' A term carrying a covariance label returns an empty list.
#' @param term A built [RandomTerm()].
#' @param ... Unused.
#' @return A list of entries, as [term_penalties()] documents.
#' @keywords internal
S7::method(term_penalties, RandomTerm) <- function(term, ...) {
  ent <- term@blueprint$entries
  if (is.null(ent)) list() else ent
}

#' @title The Covariance Label of a Random-Effect Term
#' @name term_tag.RandomTerm
#' @description
#' Returns the label written between two bars, or `NA_character_` where the
#' formula carries one bar.
#' @details
#' It is recorded by the constructor, from the middle position of
#' `~ lhs | tag | g`, and is a name that is never evaluated.
#' @param term A [RandomTerm()], built or not.
#' @param ... Unused.
#' @return A single string, or `NA_character_`.
#' @seealso [term_tag()] for the generic, [random()] for the formula from which
#'   it is read.
#' @keywords internal
S7::method(term_tag, RandomTerm) <- function(term, ...) .random_tag(term)

#' @title The Grouping of a Random-Effect Term
#' @name term_group.RandomTerm
#' @description
#' Returns the grouping expression, its levels, the number of within-group
#' columns and their names, or `NULL` on an unbuilt term.
#' @details
#' These come from the blueprint, where the build recorded them. `dim` is
#' the within-group column count, so the block's level \eqn{i} occupies
#' positions \eqn{(i-1)d + 1} to \eqn{id}, which is the group-major order that
#' the blockwise penalty reads. `names` holds the column names of the
#' within-group design.
#' @param term A built [RandomTerm()].
#' @param ... Unused.
#' @return A list with `expr`, `levels`, `dim` and `names`, or `NULL`.
#' @seealso [term_group()] for the generic.
#' @keywords internal
S7::method(term_group, RandomTerm) <- function(term, ...) {
  bp <- term@blueprint
  if (!length(bp) || is.null(bp$glevels)) return(NULL)
  d <- bp$dim
  if (is.null(d)) d <- as.integer(term_npar(term) / length(bp$glevels))
  wn <- bp$wnames
  # the within-group design's own column names, in the block's order. A
  # consumer stacking two such blocks into one covariance needs them to say
  # what each coordinate of that covariance is about; a term built before
  # the field existed has none, and a positional name is what is left.
  if (length(wn) != d) wn <- paste0("column", seq_len(d))
  list(expr = bp$gexpr, levels = as.character(bp$glevels),
       dim = as.integer(d), names = as.character(wn))
}

#' The Label of a Random-Effect Term, Normalized
#'
#' @description
#' Returns the term's `tag` as a single string, `NA_character_` where it
#' carries none.
#'
#' @details
#' The property is a character vector and is empty on a term built before the
#' label existed, so the two spellings of "no label" are read in one place.
#'
#' @param term A [RandomTerm()].
#'
#' @return A single string, possibly `NA_character_`.
#'
#' @seealso [term_tag()], which returns it.
#'
#' @keywords internal
.random_tag <- function(term) {
  if (!length(term@tag)) NA_character_ else as.character(term@tag)[1L]
}

#' A Labelled Random-Effect Term, Built Without a Penalty
#'
#' @description
#' Finishes [term_build.RandomTerm()] for a term carrying a covariance label:
#' the block and its names as usual, and no penalty entry at all.
#'
#' @details
#' The coefficients of a labelled term share a covariance block with those of
#' every other term carrying the same label and the same grouping, and that
#' block may span another equation entirely. The prior over it is therefore a
#' property of the class and not of any of its members, and a penalty built
#' here would give two priors for the same block, with no rule for which the
#' fit should read. [term_penalties()] returns an empty list, and the fitting
#' layer collects the labelled terms and declares one entry over their stacked
#' columns.
#'
#' The blueprint records what the layer needs for that: the within-group
#' dimension and the grouping's levels, beside the pieces that
#' [term_predict.RandomTerm()] already reads.
#'
#' The function rejects `hyper`, `id` and `correlated = FALSE` on a labelled
#' term, and a `distrib` that is a list or is univariate. A `distrib` given
#' here is the class's joint prior. Its dimension is the total over the
#' class and cannot be checked from one term, so that check is left to the
#' layer. The remaining checks, which verify that it is a continuous
#' \pkg{distributions7} object whose location is held, are made here, where the
#' message can name the argument.
#'
#' @param term The term being built.
#' @param Z The block.
#' @param cn Its column names.
#' @param parts The formula's three pieces.
#' @param g The grouping factor.
#' @param tt,mf,contr The within-group terms, model frame and contrasts.
#' @param d The number of within-group columns.
#' @param wn Their names.
#'
#' @return The built term, with no penalty.
#'
#' @seealso [term_build.RandomTerm()], its caller; [term_tag()].
#'
#' @keywords internal
.random_labelled <- function(term, Z, cn, parts, g, tt, mf, contr, d, wn) {
  if (length(term@hyper)) {
    stop(sprintf(paste0(
      "'hyper' names a hyperparameter of this term's own penalty, and a term\n",
      "  carrying the label '%s' has none: its coefficients share a covariance\n",
      "  block with the other terms carrying that label, and the block's\n",
      "  hyperparameters belong to the class."), .random_tag(term)),
      call. = FALSE)
  }
  # and 'id' for the same reason, said the same way: it shares a
  # hyperparameter of this term's own penalty, which a labelled term has not
  # got. The two bars are already the sharing here, of the whole block.
  if (length(term@ids)) {
    stop(sprintf(paste0(
      "'id' shares a hyperparameter of this term's own penalty, and a term\n",
      "  carrying the label '%s' has none: the covariance block that it shares is\n",
      "  the class's, and so are the hyperparameters over it."),
      .random_tag(term)), call. = FALSE)
  }
  if (!isTRUE(term@correlated)) {
    stop(sprintf(paste0(
      "'correlated = FALSE' and the label '%s' specify two things at once:\n",
      "  that the columns of this term are independent of each other, and that\n",
      "  they are correlated with the columns of another term. The label puts\n",
      "  every column of the term into the shared block."),
      .random_tag(term)), call. = FALSE)
  }
  if (!is.null(term@distrib)) {
    if (is.list(term@distrib)) {
      stop(paste0(
        "'distrib' is a list, which gives one distribution per within-group\n",
        "  column and makes the columns independent. A labelled term's prior\n",
        "  is the class's joint one, so it is a single multivariate\n",
        "  distribution over the whole shared block."), call. = FALSE)
    }
    if (!identical(term@distrib@dimension, "multivariate")) {
      stop(paste0(
        "'distrib' is univariate, and the prior of a labelled term is the\n",
        "  joint one of its class: it describes the effects of one group over\n",
        "  every column that the label collects, so it is multivariate.\n",
        "  A joint law whose coordinates are of different families -- a\n",
        "  gaussian intercept correlated with a Student t slope -- is a\n",
        "  copula, which this toolkit does not have; an elliptical family\n",
        "  cannot give it, because a multivariate gaussian has gaussian margins\n",
        "  and a multivariate t has t margins with one shared nu."),
        call. = FALSE)
    }
    # the dimension is the CLASS's total and is checked by the fitting layer;
    # everything else about the prior is asked here
    .random_check_prior(term@distrib, term@distrib@n_dim, "'distrib'")
  }
  term@X <- Z
  term@coef_names <- cn
  term@blueprint <- list(
    gexpr = parts$group, glevels = levels(g),
    terms = stats::delete.response(tt),
    xlev = stats::.getXlevels(tt, mf),
    contrasts = contr,
    entries = list(), dim = d, wnames = wn
  )
  term@penalty <- NULL
  term
}

#' @title A Random-Effect Block at New Rows
#' @name term_predict.RandomTerm
#'
#' @description
#' Rebuilds the within-group design at `newdata` and interacts it with the
#' group indicators of the levels recorded at build time, so the block has
#' the same columns in the same order however few levels the new rows happen to
#' use. By default, a level that was not present at build time is rejected.
#'
#' @details
#' An unseen group has no fitted coefficient, so the prediction for it cannot
#' use one. The message names the level. A missing grouping value is
#' rejected in the same way. Predicting a new group's response means predicting
#' at the population value, which is the model without this term's
#' contribution.
#'
#' The block comes back sparse, as the fitted one is.
#'
#' @param term A built [RandomTerm()]. An unbuilt one signals the error
#'   `"the term has not been built; call term_build(term, data) first."`.
#' @param newdata A data frame carrying the grouping variable and the
#'   within-group covariates. Its grouping factor need carry only the levels
#'   its own rows use.
#' @param unseen How to treat a row with a level that was not present at build
#'   time, or with a missing grouping value: `"error"` (the default) signals the
#'   error described above, and `"zero"` gives the row zeros, which is the
#'   prediction with that group's effect at zero. The calls
#'   `predict(random = "zero")` and `predict(random = "marginal")` of
#'   `statmodels7` pass this value.
#' @param ... Unused.
#'
#' @return A `dgCMatrix` of `nrow(newdata)` rows and [term_npar()] columns,
#'   with the term's coefficient names as column names.
#'
#' @seealso [term_predict()] for the generic, [term_build.RandomTerm()] for
#'   what recorded the levels.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(9), g = factor(rep(c("a", "b", "c"), 3)))
#' b <- term_build(random(~ 1 | g), dd)
#'
#' # A subset using two levels still gets all three columns.
#' nd <- droplevels(dd[dd$g != "c", ])
#' c(levels_here = nlevels(nd$g), cols = ncol(term_predict(b, nd)))
#'
#' # On the fitting data it returns the block itself.
#' all.equal(term_predict(b, dd), term_matrix(b))
#'
#' # A level absent from the fitting data has no coefficient, so it is
#' # rejected.
#' bad <- dd
#' levels(bad$g) <- c("a", "b", "zz")
#' try(term_predict(b, bad))
#'
#' # With unseen = "zero", such a row gets no effect of its own.
#' as.matrix(term_predict(b, bad, unseen = "zero"))[bad$g == "zz", ]
#'
#' @keywords internal
S7::method(term_predict, RandomTerm) <- function(term, newdata,
                                                 unseen = c("error", "zero"),
                                                 ...) {
  unseen <- match.arg(unseen)
  .assert_built(term)
  bp <- term@blueprint
  g <- .random_group(bp$gexpr, newdata, levels = bp$glevels, unseen = unseen)
  Z <- .random_block(g, .random_within(bp, newdata))
  colnames(Z) <- term@coef_names
  Z
}

# The within-group design at new rows, from the blueprint the build recorded.
# A smooth random effect evaluates its smooth there and applies the map the
# build computed; every other one rebuilds its model matrix.
.random_within <- function(bp, newdata) {
  if (identical(bp$kind, "smooth")) {
    W <- as.matrix(term_predict(bp$inner, newdata)) %*% bp$transform
    if (isTRUE(bp$intercept)) W <- cbind(1, W)
    dimnames(W) <- list(NULL, bp$wnames)
    return(W)
  }
  mf <- stats::model.frame(bp$terms, newdata, na.action = stats::na.pass,
                           xlev = bp$xlev)
  stats::model.matrix(bp$terms, mf, contrasts.arg = bp$contrasts)
}

#' @title The Within-Group Design of a Random-Effect Term
#' @name term_within.RandomTerm
#' @description
#' Returns the design that one group's effect multiplies, at new rows: the
#' intercept, the covariates of a random slope, one column per coordinate of
#' the effect.
#'
#' @details
#' The design is read from the blueprint that the build recorded, so a factor
#' keeps its levels and contrasts. It does not depend on the grouping, so a
#' caller can integrate the effect of a group that is absent from the fitting
#' data: the effect is \eqn{z_i^\top b} with \eqn{z_i} this row and \eqn{b}
#' drawn from the term's prior.
#'
#' @param term A built [RandomTerm()].
#' @param newdata A data frame carrying the within-group covariates.
#' @param ... Unused.
#'
#' @return A numeric matrix of `nrow(newdata)` rows and `term_group(term)$dim`
#'   columns, named as `term_group(term)$names`.
#'
#' @seealso [term_within()] for the generic, [term_group()] for the layout.
#'
#' @examples
#' dd <- data.frame(x = c(0.5, 1, 2), g = factor(c("a", "b", "a")))
#' b <- term_build(random(~ 1 + x | g), dd)
#' term_within(b, data.frame(x = c(3, 4), g = c("new", "a")))
#'
#' @keywords internal
S7::method(term_within, RandomTerm) <- function(term, newdata, ...) {
  .assert_built(term)
  W <- .random_within(term@blueprint, newdata)
  gr <- term_group(term)
  if (!is.null(gr) && ncol(W) == length(gr$names)) colnames(W) <- gr$names
  attr(W, "assign") <- NULL
  attr(W, "contrasts") <- NULL
  W
}
