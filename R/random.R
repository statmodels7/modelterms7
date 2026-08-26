#' @include term_classes.R generics.R
NULL

#' @title S7 Class for Grouped Random-Effect Terms
#' @name RandomTerm
#'
#' @description
#' The subclass of [additive_term()] holding grouped coefficients with a
#' distribution attached to them: the within-group design interacted with the
#' group indicators, one coefficient per group and per within-group column,
#' and a penalty carrying that distribution's negative log-density. [random()]
#' constructs it.
#'
#' @details
#' # The three properties of its own
#'
#' `formula` is the bar formula as given, `~ 1 | g` or `~ x | g`, kept with its
#' environment. `correlated` says whether the **default** Gaussian lets the
#' within-group effects depend on each other; it is read only where `distrib`
#' is `NULL`, the two saying the same thing.
#'
#' `distrib` is the effects' distribution as supplied, or `NULL` for the
#' default. What the build turns it into is a \pkg{penalties7} penalty, read
#' through [term_penalty()] or [term_penalties()], so the hyperparameter names
#' and their bounds all come from the distribution.
#'
#' # The block is sparse by construction
#'
#' A row belongs to one group, so the block has a density of \eqn{1/m} and is
#' always a `dgCMatrix`. [random()] takes no `sparse` argument for that reason,
#' and passing one is an error.
#'
#' @inheritParams additive_term
#' @param formula The bar formula, `~ 1 | g` or `~ x | g`, with the
#'   within-group design on the left and the grouping variable on the right.
#' @param correlated A single logical: whether the default Gaussian lets the
#'   within-group effects correlate. Read only when `distrib` is `NULL`.
#' @param distrib The effects' distribution, a \pkg{distributions7} object or
#'   a list of them with one per within-group column, or `NULL` for the
#'   default Gaussian.
#'
#' @return An S7 object of class `RandomTerm`, inheriting from
#'   [additive_term()] and [model_term()], with the three properties above
#'   beside the ten they supply.
#'
#' @seealso [random()], the constructor; [term_penalties()] for the entries a
#'   built one declares; [edf()] for what a fitted random effect spends.
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
    distrib = S7::class_any
  )
)

#' Grouped Random-Effect Term
#'
#' @description
#' Random intercepts and slopes for a grouping factor:
#' `random(~ 1 | g)` builds one coefficient per level of `g`,
#' and `random(~ x | g)` one intercept and one slope per level, with
#' the distribution of the effects attached as the penalty on those
#' coefficients. That is what a random effect is under penalized likelihood.
#'
#' @details
#' The left side of the bar is an ordinary one-sided formula for the
#' within-group design, with the usual intercept convention:
#' `~ x | g` carries an intercept and a slope per group and
#' `~ 0 + x | g` the slope alone. The block interacts that design
#' with the group indicators, ordered group by group, so the coefficients
#' of one group are adjacent.
#'
#' The constructor asks for two things: the formula and the distribution of
#' the effects. Which chart the hyperparameters ride, what they are called,
#' how many there are and where the log-density has a kink are all properties
#' of that distribution, read off it at build time.
#'
#' @section The distribution of the effects:
#' `distrib` is `NULL`, a \pkg{distributions7} object, or a list of
#' them with one per within-group column.
#'
#' A multivariate distribution of the within-group dimension lets the effects
#' of one group depend on each other, its matrix parameter carrying the
#' dependence: `mvgaussian_distrib(2, omega = ar1(2))` is a prior whose
#' precision is autoregressive, `mvstudent_t_distrib(2)` a heavy-tailed
#' one. Correlation is available exactly for the families that carry a matrix
#' parameter: a location block as long as the dimension, together with a
#' covariance, precision or scale matrix. The term reads that property off the
#' family, so a family added later is covered without an edit here.
#'
#' A univariate distribution makes the effects independent, the penalty being
#' the product of the densities. With more than one within-group column it is
#' a template: one copy per column, each with its own hyperparameters, since
#' an intercept and a slope are quantities of different units and a shared
#' scale would price them against each other. A list of distributions gives
#' one per column explicitly, when the columns want different priors.
#'
#' The default is Gaussian: `gaussian1_distrib` at one column, so the
#' hyperparameter IS the standard deviation of the effects; the multivariate
#' Gaussian on an unstructured covariance when there are several and
#' `correlated = TRUE`; the template of the first when
#' `correlated = FALSE`, one standard deviation per column.
#'
#' Whatever it is, the distribution is centered, its location parameters held
#' with [distributions7::fixed()]. A free mean in the effects is
#' confounded with the intercept of the equation the term sits in, so a free
#' location is rejected at build time with a message naming it. The value it
#' is held at is usually zero and is not policed: the model is identified
#' whatever it is. Where the prior is a transformation of another family the
#' parameter is the mean on the original scale, so
#' `fixed(transformation(gamma2_distrib(), log_transform()), mu = 1)` is a
#' log-gamma prior whose own mean is \eqn{\psi(a) - \log a}, within
#' \eqn{\sigma^2/2} of zero and exactly zero in the limit.
#'
#' A distribution used as a penalty gives joint-mode (penalized likelihood)
#' estimation of the effects. With Gaussian effects and a Gaussian response
#' the Laplace approximation behind a marginal criterion is exact, so the
#' variance component it returns is the marginal estimate; with any other
#' prior it is an approximation, and the marginal likelihood is not computed
#' here.
#'
#' @section The hyperparameters:
#' They are the distribution's own free parameters, and every one of them is
#' estimated unless it is held. There are two ways to hold one, and the fit is
#' the same either way; what differs is what gets reported. Holding it inside
#' the distribution, `fixed(pseudohuber_distrib(), mu = 0, nu = 2)`, removes
#' it: it becomes a constant of the prior and appears nowhere among the
#' model's hyperparameters. Naming it in `hyper` keeps it, reported as held at
#' the value given, as a penalized term's own hyperparameter argument does.
#'
#' A smooth prior's hyperparameters are estimated by a marginal criterion. A
#' prior whose log-density has a kink, a Laplace or an elastic net, has none a
#' marginal criterion can reach, and its hyperparameter is chosen by a path on
#' a prediction criterion instead.
#'
#' Every estimated hyperparameter is reported with a standard error and an
#' interval, shape parameters included. Where one is absent the cause is the
#' point the run ended at: a criterion with no maximum there leaves a
#' curvature of the wrong sign, and no interval follows from it. A shape
#' escaping toward a limit is the common case.
#'
#' Which parametrization of a family is used matters here in a way it does
#' not elsewhere. The centered skew normal
#' ([distributions7::skewnormal2_distrib()]) carries the skewness
#' itself, and its map to the direct parametrization is not twice
#' differentiable at zero skewness: the first derivatives have a finite limit
#' there and the second ones grow like \eqn{\gamma_1^{-2/3}}. A marginal
#' criterion reads the second, and the symmetric bounds put the starting
#' value at exactly that point, so the direct parametrization
#' ([distributions7::skewnormal1_distrib()]) is the one to use as a
#' prior; its derivatives at \eqn{\alpha = 0} are ordinary numbers.
#'
#' How well a shape parameter is estimated depends on how many groups there
#' are, since it is read off that many latent values and the prior shrinks
#' them. Measured on effects drawn from a standard Student t with four
#' degrees of freedom, twelve observations per group and unit residual
#' standard deviation, the prior being a Student t with \eqn{\nu} free and
#' the criterion [statmodels7::reml()]:
#'
#' | groups | \eqn{\hat\nu} | \eqn{\hat\sigma} |
#' | --- | --- | --- |
#' | 20 | 5.97e+04 | 0.769 |
#' | 100 | 1.98 | 0.817 |
#' | 500 | 2.65 | 0.923 |
#'
#' At twenty groups the shape escapes to the Gaussian limit and only the
#' scale is really being fitted. From a hundred it stays finite, and the
#' profile is decisive about that much: with \eqn{\nu} held, the criterion is
#' -1922.4 at 3, -1923.7 at 4, -1924.9 at 5 and -1936.1 in the Gaussian
#' limit. What it is not decisive about is the value, the profile being flat
#' enough over the small integers that a single sample locates \eqn{\nu} to
#' little better than its order of magnitude. Estimate a shape from a hundred
#' groups or so, hold it below that, and read the estimate as a statement
#' about the tail rather than a measurement of it.
#'
#' A pseudo-Huber's \eqn{\nu} is the weaker case, being the point at which
#' the loss stops being quadratic; at 40 groups it escapes.
#'
#' Prediction maps new data onto the levels seen at build time; a level
#' the term has not seen is rejected.
#'
#' A random effect is not standardized, and there is no `standardize`
#' argument to ask for it with; passing one is an error. Its columns are
#' grouping indicators and its hyperparameter is a variance component with a
#' meaning of its own. Dividing each coefficient by the spread of its
#' indicator would weight the effects by the sizes of the groups, which
#' changes the model itself.
#'
#' @section The block and its penalty:
#' With \eqn{m} levels and a within-group design \eqn{Z_i} of \eqn{d}
#' columns, the block is the interaction of that design with the group
#' indicators, ordered group by group,
#'
#' \deqn{Z = \operatorname{diag}(Z_1, \dots, Z_m),
#'   \qquad b = (b_1', \dots, b_m')',}
#'
#' so the \eqn{d} coefficients of one group occupy adjacent positions and the
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
#' mixed-model equation at the variance ratio \eqn{\Sigma} encodes, so the
#' minimizer is the best linear unbiased predictor. At \eqn{d = 1} and
#' \eqn{\Sigma = \sigma_b^2} it is the ridge, up to the constant that makes
#' \eqn{\sigma_b} estimable.
#'
#' @param formula A bar formula, `~ lhs | g`, with `g`
#'   evaluating to the grouping variable in the data.
#' @param distrib The distribution of the effects: `NULL` (the default
#'   Gaussian), a \pkg{distributions7} object, or a list of them with one per
#'   within-group column.
#' @param correlated Logical; whether the default Gaussian lets the
#'   within-group effects correlate. It is an error together with
#'   `distrib`, which says the same thing and more.
#' @param label A single non-empty string prefixed to the coefficient
#'   names.
#' @param hyper The hyperparameters of the effects' distribution to hold, as
#'   a named vector or list; those not named are estimated. The names are the
#'   distribution's own parameters, with the within-group column appended
#'   where there is one copy per column. A name the penalty does not carry is
#'   reported when the term is built, which is where the penalty first exists.
#' @param ... Unused. A named argument here is reported by name, so a removed
#'   one such as `precision` or `kinks` gets a message saying what replaced
#'   it.
#'
#' @return An object of class [RandomTerm()] (a specification;
#'   see [term_build()]).
#'
#' @examples
#' dd <- data.frame(y = rnorm(9), x = rnorm(9),
#'                  g = factor(rep(c("a", "b", "c"), 3)))
#'
#' # one variance component, reported as a standard deviation
#' term_penalty(term_build(random(~ 1 | g), dd))@params
#'
#' # intercepts and slopes, correlated: the covariance of the effects
#' built <- term_build(random(~ x | g), dd)
#' term_coef_names(built)
#' term_penalty(built)@params
#'
#' # independent, one standard deviation per column
#' vapply(term_penalties(term_build(random(~ x | g, correlated = FALSE), dd)),
#'        function(e) e$name, "")
#'
#' # a heavy-tailed prior, held at four degrees of freedom
#' t4 <- distributions7::fixed(distributions7::student_t1_distrib(),
#'                             mu = 0, nu = 4)
#' term_penalty(term_build(random(~ 1 | g, distrib = t4), dd))@params
#'
#'
#' # Fitted. The data are simulated from a known truth, so the
#' # estimates below can be read against it.
#' if (requireNamespace("statmodels7", quietly = TRUE)) {
#'   set.seed(6)
#'   fd <- data.frame(gr = factor(rep(1:20, each = 15)), x = rnorm(300))
#'   bb <- rnorm(20, sd = 0.8)
#'   fd$y <- 1 + 0.5 * fd$x + bb[fd$gr] + rnorm(300, sd = 0.5)
#'   cf <- coef(statmodels7::statmod(y ~ x + random(~1 | gr),
#'                                   distributions7::gaussian1_distrib(), fd))$mu
#'   # truth: a slope of 0.5, and the shrunken effects track the ones drawn
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
                   label = "random", hyper = NULL, ...) {
  if (!inherits(formula, "formula") || length(formula) != 2L) {
    stop("'formula' must be a one-sided bar formula, e.g. ~ 1 | g.",
         call. = FALSE)
  }
  e <- formula[[2L]]
  if (!is.call(e) || !identical(e[[1L]], as.name("|"))) {
    stop("'formula' must contain a grouping bar, e.g. ~ 1 | g.",
         call. = FALSE)
  }
  if (!is.logical(correlated) || length(correlated) != 1L ||
      is.na(correlated)) {
    stop("'correlated' must be TRUE or FALSE.", call. = FALSE)
  }
  .random_retired(...)
  if (!is.null(distrib) && !missing(correlated)) {
    stop(paste0(
      "'correlated' and 'distrib' both say how the effects are\n",
      "  distributed. Say it once: a multivariate distribution lets them\n",
      "  correlate and a univariate one does not."), call. = FALSE)
  }
  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !nzchar(label)) {
    stop("'label' must be a single non-empty character string.",
         call. = FALSE)
  }
  RandomTerm(label = label, formula = formula, correlated = correlated,
             distrib = distrib,
             hyper = as_hyper(hyper, label),
             X = NULL, coef_names = character(0),
             blueprint = list(), penalty = NULL)
}

#' The Arguments random() No Longer Takes
#'
#' @description
#' Signals an error naming a removed argument, which the dots would otherwise
#' swallow in silence.
#'
#' @details
#' `precision` was a second spelling of a multivariate Gaussian whose
#' matrix parameter is a precision, and the first spelling did not say which
#' matrix the structure was. `kinks` was derived from the effects'
#' distribution by \pkg{penalties7} all along, and the default here overrode
#' that derivation, so a Laplace prior declared none and a fitting layer sent
#' its block to the scheme that cannot solve it.
#'
#' @param ... The dots of [random()].
#'
#' @return Invisibly `NULL`; called for its error.
#'
#' @keywords internal
.random_retired <- function(...) {
  dots <- list(...)
  if (!length(dots)) return(invisible(NULL))
  nm <- names(dots)
  if (is.null(nm) || !any(nzchar(nm))) {
    stop("random() takes no unnamed arguments beyond 'formula'.",
         call. = FALSE)
  }
  if ("precision" %in% nm) {
    stop(paste0(
      "'precision' has been removed: a structured precision is the matrix\n",
      "  parameter of a multivariate Gaussian, so it is written as one --\n",
      "  distrib = fixed(mvgaussian_distrib(d, omega = <structure>),\n",
      "                  mu1 = 0, ..., mud = 0)\n",
      "  which also says which of the two matrices the structure is."),
      call. = FALSE)
  }
  if ("kinks" %in% nm) {
    stop(paste0(
      "'kinks' has been removed: they are derived from the effects'\n",
      "  distribution by penalties7::distrib_kinks(), which is what every\n",
      "  other penalized term relies on."), call. = FALSE)
  }
  stop(sprintf("random() has no argument '%s'.", nm[nzchar(nm)][1L]),
       call. = FALSE)
}

.random_group <- function(expr, data, levels = NULL) {
  v <- eval(expr, data, baseenv())
  if (is.null(levels)) return(factor(v))
  f <- factor(v, levels = levels)
  if (any(is.na(f) & !is.na(v))) {
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
  # one entry per (row, within-column) pair: row i lands in the d columns of
  # its own group, and in no other
  Matrix::sparseMatrix(
    i = rep(seq_len(n), each = d),
    j = as.vector(t(outer((gi - 1L) * d, seq_len(d), `+`))),
    x = as.vector(t(W)),
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
#' Follows `parent_distrib` to the distribution a wrapper wraps, so that
#' a question about the family is asked of the family.
#'
#' @details
#' The property is asked for with `S7::prop_names()` instead of testing the
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
#' `TRUE` for a family with a location block as long as its dimension and a
#' matrix parameter, which together are what a centered prior on \eqn{R^d}
#' needs.
#'
#' @details
#' The question is a property of the family, so a multivariate family added
#' later is covered without an edit here. It is
#' read off `params_interpretation`, the same declaration a data-based
#' starting value is built from. It excludes the simplex-valued families,
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

#' Whether a Multivariate Family Answers Its Mixed Block
#'
#' @description
#' `TRUE` when `distrib_cross_y` comes from the family itself. The
#' multivariate base class rejects that generic, so a family that has not
#' overridden it answers `FALSE`.
#'
#' @details
#' A marginal criterion reads that block to estimate the covariance of the
#' effects, so a family without one can be fitted at held hyperparameters and
#' not at estimated ones. The owning class of a method is read through its
#' signature and compared by name, never with `identical()`, which is
#' object identity and fails whenever a package's code is re-evaluated rather
#' than loaded.
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
#' What is rejected is a free location. It is confounded with the intercept of
#' the equation the term sits in, leaving a flat direction along which the fit
#' has no answer.
#'
#' A location held at a value is identified whatever that value is, and it is
#' not policed: it shrinks the effects toward that value, which is a modeling
#' statement. Nor could the value be policed in general. Where the prior is a
#' transformation of another family, the parameter interpreted as its mean is
#' the mean on the original scale, and holding the mean of a gamma at one is
#' what centers its logarithm.
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
        "%s ('%s') carries no matrix parameter over R^%d, so it cannot say\n",
        "  how the effects of one group depend on each other. Correlated\n",
        "  effects need a family with a location block and a covariance,\n",
        "  precision or scale matrix, such as mvgaussian_distrib() or\n",
        "  mvstudent_t_distrib()."),
        what, base@distrib_name, base@n_dim), call. = FALSE)
    }
    if (d@n_dim != dim_needed) {
      stop(sprintf(paste0(
        "%s is %d-variate and the bar's left side has %d columns; a\n",
        "  multivariate prior describes the effects of ONE group."),
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
      "  intercept of the equation the term sits in. Hold it with\n",
      "  distributions7::fixed(): at zero for a family on the line, and at\n",
      "  whatever centers the effects where the prior is a transformation of\n",
      "  another family, the parameter being the mean on the original scale."),
      what, free[1L]), call. = FALSE)
  }
  invisible(TRUE)
}

#' The Centered Gaussian Defaults
#'
#' @description
#' `gaussian1_distrib` at one within-group column, so the hyperparameter
#' IS the standard deviation of the effects, and the multivariate Gaussian on
#' an unstructured covariance for several correlated ones.
#'
#' @details
#' The structure's role is declared, `"covariance"` here. A structure left at
#' `"either"` does not say which matrix of the prior it is, and the two cannot
#' be read interchangeably: they differ in the sign of the log-determinant
#' term.
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
  mv <- distributions7::mvgaussian_distrib(
    d, sigma = parameters7::log_cholesky(d, role = "covariance"))
  do.call(distributions7::fixed,
          c(list(mv), stats::setNames(as.list(rep(0, d)),
                                      paste0("mu", seq_len(d)))))
}

#' The Penalty Entries a Random-Effect Term Declares
#'
#' @description
#' One entry over the whole block when the prior reads a group's effects
#' together, and one per within-group column when they are independent.
#'
#' @details
#' The coefficients are ordered group by group, so column \eqn{j} is the
#' stride \eqn{j, d+j, 2d+j, \dots}. Those positions are named as a subset of
#' the term's own parameters, never selected with a map, and that is what
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
      return(list(list(name = "", index = seq_len(m * d), penalty = pen,
                       fixed = list(), n_values = list(), values = list(),
                       min_ratio = numeric(0))))
    }
    prior <- rep(list(prior), d)
  }
  lapply(seq_len(d), function(j) {
    list(name = wnames[j], index = seq(j, by = d, length.out = m),
         penalty = penalties7::distrib_penalty(prior[[j]], n_coef = m),
         fixed = list(), n_values = list(), values = list(),
         min_ratio = numeric(0))
  })
}

#' Where a Held Hyperparameter Belongs
#'
#' @description
#' Splits the term's `hyper` over its penalty entries and checks every
#' name against the penalty that carries it.
#'
#' @details
#' A name is qualified by the within-group column where there is one penalty
#' per column. An unqualified one is an error that lists what there is, not a
#' value recycled over every column: a caller who wants the same value
#' everywhere writes it into the distribution, where it stops being a
#' hyperparameter at all, and silent recycling is the trap this file's
#' history records for `ifelse`.
#'
#' @param entries The entries from [.random_entries()].
#' @param hyper The term's `hyper`, already normalized.
#' @param label The term's label, for the message.
#'
#' @return The entries, with `fixed` filled in and checked.
#'
#' @keywords internal
.random_hyper <- function(entries, hyper, label) {
  qual <- function(en) {
    if (nzchar(en$name)) paste0(en$penalty@params, ".", en$name)
    else en$penalty@params
  }
  avail <- unlist(lapply(entries, qual), use.names = FALSE)
  unknown <- setdiff(names(hyper), avail)
  if (length(unknown)) {
    stop(sprintf(paste0(
      "'hyper$%s' is not a hyperparameter of the effects' distribution in\n",
      "  '%s'. It carries: %s."),
      unknown[1L], label, paste(avail, collapse = ", ")), call. = FALSE)
  }
  lapply(entries, function(en) {
    keep <- hyper[intersect(names(hyper), qual(en))]
    if (nzchar(en$name) && length(keep)) {
      names(keep) <- sub(sprintf("\\.\\Q%s\\E$", en$name), "", names(keep))
    }
    # several values are a grid for a PATH to visit, and only a penalty with a
    # kink is swept along one -- the same three checks a penalized constructor
    # runs, in the same order
    vals <- check_values(keep, en$penalty, label)
    reject_pathless_values(vals, en$penalty, label)
    en$fixed <- check_hyper(keep, en$penalty, label)
    en$values <- vals
    en
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
#' the block is \eqn{\mathrm{diag}(Z_1, \dots, Z_m)}, ordered **group-major**,
#' so the \eqn{d} coefficients of one group are adjacent. It is built as a
#' `dgCMatrix`: a row belongs to one group, so the density is \eqn{1/m}.
#'
#' The coefficient names are `label.level.column`, so `random(~ x | g)` over
#' three levels gives `random.a.(Intercept)`, `random.a.x`,
#' `random.b.(Intercept)` and so on, which is the group-major order read off.
#'
#' # The penalty, and what the build checks
#'
#' Where `distrib` is `NULL` the default is chosen here: a centered
#' `gaussian1_distrib` at one column or under `correlated = FALSE`, and a
#' centered multivariate Gaussian on an unstructured covariance for several
#' correlated ones.
#'
#' Whatever the distribution, its location parameters must be **held**. A free
#' location is confounded with the intercept of the equation the term sits in,
#' and the build rejects it with a message naming the parameter and the fix. A
#' multivariate distribution must also match the within-group dimension, and
#' one that carries no matrix parameter cannot express correlated effects at
#' all.
#'
#' Any value named in `hyper` is checked here too, against the penalty's own
#' names, this being the first point at which the penalty exists.
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
#' # A free location is refused, naming the parameter.
#' try(term_build(random(~ 1 | g,
#'                       distrib = distributions7::gaussian1_distrib()), dd))
#'
#' @keywords internal
S7::method(term_build, RandomTerm) <- function(term, data, ...) {
  e <- term@formula[[2L]]
  g <- .random_group(e[[3L]], data)
  m <- nlevels(g)
  if (m < 2L) {
    stop("the grouping variable must have at least two levels.",
         call. = FALSE)
  }

  wf <- stats::as.formula(call("~", e[[2L]]),
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

  Z <- .random_block(g, W)
  cn <- .random_names(term@label, levels(g), colnames(W))
  colnames(Z) <- cn

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
          "'distrib[[%d]]' is multivariate. A list gives ONE distribution\n",
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
  entries <- .random_hyper(.random_entries(prior, d, m, colnames(W)),
                           term@hyper, term@label)

  term@X <- Z
  term@coef_names <- cn
  term@blueprint <- list(
    gexpr = e[[3L]], glevels = levels(g),
    terms = stats::delete.response(tt),
    xlev = stats::.getXlevels(tt, mf),
    contrasts = contr,
    entries = entries
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
#' One entry over the whole block where the effects of a group are read
#' together, and one per within-group column where they are independent.
#' @param term A built [RandomTerm()].
#' @param ... Unused.
#' @return A list of entries, as [term_penalties()] documents.
#' @keywords internal
S7::method(term_penalties, RandomTerm) <- function(term, ...) {
  ent <- term@blueprint$entries
  if (is.null(ent)) list() else ent
}

#' @title A Random-Effect Block at New Rows
#' @name term_predict.RandomTerm
#'
#' @description
#' Rebuilds the within-group design at `newdata` and interacts it with the
#' group indicators **of the levels recorded at build time**, so the block has
#' the same columns in the same order however few levels the new rows happen to
#' use. A level the term never saw is refused.
#'
#' @details
#' The refusal is the right answer rather than a limitation: a coefficient was
#' never fitted for an unseen group, so there is nothing to predict with. The
#' message names the level. Predicting a new group's response means predicting
#' at the population value, which is the model without this term's
#' contribution.
#'
#' The block comes back sparse, as the fitted one is.
#'
#' @param term A built [RandomTerm()]. An unbuilt one throws
#'   `"the term has not been built; call term_build(term, data) first."`.
#' @param newdata A data frame carrying the grouping variable and the
#'   within-group covariates. Its grouping factor need carry only the levels
#'   its own rows use.
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
#' # A level the fit never saw has no coefficient, so it is refused.
#' bad <- dd
#' levels(bad$g) <- c("a", "b", "zz")
#' try(term_predict(b, bad))
#'
#' @keywords internal
S7::method(term_predict, RandomTerm) <- function(term, newdata, ...) {
  .assert_built(term)
  bp <- term@blueprint
  g <- .random_group(bp$gexpr, newdata, levels = bp$glevels)
  mf <- stats::model.frame(bp$terms, newdata, na.action = stats::na.pass,
                           xlev = bp$xlev)
  W <- stats::model.matrix(bp$terms, mf, contrasts.arg = bp$contrasts)
  Z <- .random_block(g, W)
  colnames(Z) <- term@coef_names
  Z
}
