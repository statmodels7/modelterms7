#' @include term_classes.R generics.R structural.R regime.R segmented.R
NULL

#' @title The Levels of a Likelihood-Shaped Structural Term
#'
#' @description
#' The shifts a term of the likelihood shape adds to its equation's predictor,
#' one per mixture component, in the order the columns of [term_posterior()]
#' carry the components. A fitting layer reads the two together to assemble
#' Fisher's identity.
#'
#' @details
#' By Fisher's identity the derivative of a likelihood mixed over latent
#' states, in **any** predictor the model carries, is the posterior-weighted
#' derivative of the ordinary one, each component read at the predictor
#' shifted by its own level:
#'
#' \deqn{\frac{\partial L}{\partial \eta_{q,t}}
#'   = \sum_k \gamma_t(k)\,
#'     \frac{\partial \ell(y_t; \theta_t(k))}{\partial \eta_q}.}
#'
#' [term_posterior()] supplies the \eqn{\gamma_t(k)} and this supplies the
#' shifts, so a fitting layer assembles the identity without reading the
#' term's internals.
#'
#' For [regime()] the levels are the ordered regime means, one number per
#' component. For a marginal break-point term of the step kind they are the
#' sums of the changes of level over the active break-points, one number per
#' side pattern.
#'
#' # A shift may vary by observation
#'
#' The quadrature nodes of a marginal [seg()] or [jseg()] term shift each
#' observation by its own hinge value, so the method may return a **matrix**
#' of one row per observation and one column per component. A caller reads a
#' column of it wherever it would read a constant level, and must accept both
#' shapes.
#'
#' The method on [structural_term()] throws, naming the class: a term of the
#' filter shape reports a predictor instead of components, and answers
#' [term_filter()].
#'
#' @param term A built structural term of the likelihood shape.
#' @param psi The term's parameters on the parameter scale, named as
#'   [term_params()].
#' @param ... Passed to methods.
#'
#' @return A numeric vector of one level per component, or a numeric matrix of
#'   `n` rows and one column per component where a shift varies by
#'   observation. The order matches [term_posterior()]'s columns.
#'
#' @seealso [term_posterior()] for the weights, [term_loglik()] for the
#'   likelihood, [regime()] and [jump()] for the two implementations.
#'
#' @examples
#' # A two-regime chain: the first level and the cumulated gap.
#' term_levels(regime(2), list(level1 = 0, gap2 = 3, alr1.1 = 2, alr2.1 = -2))
#'
#' # A one-break-point step term: no shift where no break-point is active,
#' # and the change of level where one is.
#' set.seed(1)
#' dd <- data.frame(id = rep(1:3, each = 8), x = rep(1:8, 3))
#' dd$y <- rnorm(24, 2 * (dd$x >= 4.5), 0.4)
#' tm <- term_build(jump(x, psi ~ random(~ 1 | id), marginal = TRUE), dd)
#' term_levels(tm, list(m1 = 4.5, tau1 = 0.5, delta1 = 2))
#'
#' @export
#' @aliases term_levels.structural_term
term_levels <- S7::new_generic("term_levels", "term",
  function(term, psi, ...) S7::S7_dispatch())

S7::method(term_levels, structural_term) <- function(term, psi, ...) {
  stop(sprintf("the term class '%s' does not implement term_levels().",
               attr(S7::S7_class(term), "name")), call. = FALSE)
}

# Registered here rather than in regime.R because the generic is: this file
# loads after that one, so the class exists when the method is written.
#' @title Levels of a Regime Term
#' @name term_levels.RegimeTerm
#' @description
#' The regime means, ordered by construction: the first level and the
#' cumulative sums of the positive gaps.
#' @param term A [RegimeTerm()].
#' @param psi The term's parameters.
#' @param ... Unused.
#' @return A numeric vector of length `k`.
#' @keywords internal
S7::method(term_levels, RegimeTerm) <- function(term, psi, ...) {
  v <- unlist(psi[term_params(term)])
  k <- term@k
  unname(cumsum(c(v[["level1"]],
                  if (k > 1L) v[paste0("gap", seq.int(2L, k))]
                  else numeric(0))))
}


#' @title The Posterior of a Structural Term's Latent Variable
#'
#' @description
#' A summary of the latent variable a structural term integrates over, given
#' the whole sample. For a marginal break-point term it is the posterior mean
#' and standard deviation of each group's break-points. That is what a reader
#' wants from such a fit: where each group's change happened, and how sure the
#' data are about it.
#'
#' @details
#' [term_posterior()] answers the fitting layer's question, the component
#' weights Fisher's identity needs at every observation. This one answers the
#' reader's. For the marginal break-point term the two come from the same
#' decomposition: the mean and variance within an interval are those of the
#' prior truncated to it, and under quadrature the moments of the node
#' posterior.
#'
#' # What a heavy-tailed prior can refuse
#'
#' The moments are the prior's, so a prior without them has none to report. A
#' Student t below one degree of freedom has no mean on an edge interval and
#' below two no variance, and the quadrature returns `NA` there instead of a
#' number. That is a property of the prior the caller chose.
#'
#' The method on [structural_term()] throws, naming the class: [gas()] and
#' [regime()] have no continuous latent to summarize this way, a regime's
#' latent being the discrete state [term_posterior()] already reports.
#'
#' @param term A built structural term.
#' @param eta The static predictor of the equation the term sits in, one value
#'   per observation.
#' @param y The response.
#' @param logdens A function `(e, i)` returning the log-density of observation
#'   `i` at predictor `e`, as [term_loglik()] takes it.
#' @param psi The term's parameters on the parameter scale, named as
#'   [term_params()].
#' @param ... Passed to methods.
#'
#' @return A data frame with one row per group and break-point and four
#'   columns: `group`, the grouping level; `psi`, which break-point;
#'   `mean` and `sd`, the posterior moments of its position. `NA` in a moment
#'   the prior does not possess.
#'
#' @seealso [term_posterior()] for the component weights a fit reads,
#'   [jump()] and [seg()] for the terms that implement it.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(id = rep(1:3, each = 8), x = rep(1:8, 3))
#' dd$y <- rnorm(24, 2 * (dd$x >= 4.5), 0.4)
#' tm <- term_build(jump(x, psi ~ random(~ 1 | id), marginal = TRUE), dd)
#'
#' # Every group's step is at 4.5, and the posterior finds it there,
#' # with a spread well inside the prior's own 0.5.
#' term_latent(tm, rep(0, 24), dd$y,
#'             logdens = function(e, i) dnorm(dd$y[i], e, 0.4, log = TRUE),
#'             psi = list(m1 = 4.5, tau1 = 0.5, delta1 = 2))
#'
#' @export
#' @aliases term_latent.structural_term
term_latent <- S7::new_generic("term_latent", "term",
  function(term, eta, y, logdens, psi, ...) S7::S7_dispatch())

S7::method(term_latent, structural_term) <- function(term, eta, y, logdens,
                                                     psi, ...) {
  stop(sprintf("the term class '%s' does not implement term_latent().",
               attr(S7::S7_class(term), "name")), call. = FALSE)
}


#' @title S7 Class for Marginal Break-Point Terms
#' @name MarginalBreakTerm
#'
#' @description
#' The subclass of [structural_term()] holding break-points that vary by group
#' as latent variables **integrated out** of the likelihood.
#' [jump()], [seg()] and [jseg()] construct it when given `marginal = TRUE`
#' together with a `psi ~ random(~ 1 | g)` subformula.
#'
#' Its contribution is a likelihood, not a predictor, so it implements
#' [term_loglik()] and the prior over the positions is part of that likelihood
#' itself. [term_penalties()] declares nothing, and the prior's parameters
#' are estimated by plain maximum likelihood.
#'
#' @details
#' # The eight properties of its own
#'
#' `kind` is `"jump"`, `"seg"` or `"jseg"`; `var` the covariate expression;
#' `npsi` the number of break-points per group; `linear` whether the term
#' carries the linear effect as a parameter of its own, which `seg` and `jseg`
#' do.
#'
#' `group` is the grouping expression, taken from the break-point's `random()`
#' subformula. `prior` is the latent's distribution: `NULL` for the Gaussian,
#' or a \pkg{distributions7} object where `random(distrib = )` named one, and
#' its location must be fixed at zero, `m1` carrying the position.
#'
#' `spec` holds the resolved construction settings and `blueprint` the
#' grouping and the interval structure [term_build()] worked out.
#'
#' # The parameters
#'
#' They are numbered, one set per break-point: `m1`, `tau1`, `delta1` for a
#' one-break-point step term, with `m` the prior's location, `tau` its scale on
#' a log chart, and `delta` the change of level. A continuous kind adds `beta`
#' for the linear effect and `gamma1` for the change of slope.
#'
#' # What it costs
#'
#' At most **eight** break-points. The forward recursion of the step kind costs
#' \eqn{n K 2^K} and stays cheap well past that; what does not is the fitting
#' layer, which reads a posterior over the \eqn{2^K} side patterns and
#' evaluates the family once per pattern. The continuous kinds stay at one
#' break-point, a product quadrature over more being far dearer still.
#'
#' @inheritParams model_term
#' @param kind One of `"jump"`, `"seg"` or `"jseg"`.
#' @param var The covariate expression, unevaluated.
#' @param npsi The number of break-points per group, an integer between 1 and
#'   8 for `"jump"` and exactly 1 for the other two.
#' @param linear Whether the term carries the linear effect as its own
#'   parameter, `TRUE` for `"seg"` and `"jseg"`.
#' @param group The grouping expression, from the break-point's `random()`
#'   subformula.
#' @param prior The latent's distribution: `NULL` for the Gaussian, or a
#'   \pkg{distributions7} object with its location held at zero.
#' @param spec A named list of the resolved construction settings.
#' @param blueprint A named list of the resolved grouping and interval
#'   structure, empty until [term_build()] fills it.
#'
#' @return An S7 object of class `MarginalBreakTerm`, inheriting from
#'   [structural_term()] and [model_term()], with the eight properties above
#'   beside [model_term()]'s six.
#'
#' @seealso [jump()], [seg()] and [jseg()] for the constructors;
#'   [term_loglik()] for the likelihood; [term_latent()] for the posterior
#'   positions a reader wants.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(id = rep(1:3, each = 8), x = rep(1:8, 3))
#' dd$y <- rnorm(24, 2 * (dd$x >= 4.5), 0.4)
#'
#' tm <- term_build(jump(x, psi ~ random(~ 1 | id), marginal = TRUE), dd)
#' S7::S7_inherits(tm, MarginalBreakTerm)
#' c(kind = tm@kind, npsi = tm@npsi, linear = tm@linear)
#'
#' # Numbered parameters: the prior's location and scale, and the change.
#' term_params(tm)
#' vapply(term_links(tm), function(l) l@link_name, character(1))
#'
#' # The prior is part of the likelihood, so nothing is declared penalized.
#' length(term_penalties(tm))
#'
#' # A continuous kind adds the linear effect and the change of slope.
#' term_params(term_build(seg(x, psi ~ random(~ 1 | id), marginal = TRUE), dd))
#'
#' @export
MarginalBreakTerm <- S7::new_class(
  name = "MarginalBreakTerm",
  parent = structural_term,
  properties = list(
    kind = S7::class_character,
    var = S7::class_any,
    npsi = S7::class_integer,
    linear = S7::class_logical,
    group = S7::class_any,
    prior = S7::class_any,
    spec = S7::class_list,
    blueprint = S7::class_list
  )
)

# The constructor branch of the three break-point terms under
# marginal = TRUE. The prior is read off the break-point's subformula,
# which must be a single random(~1 | g) call: the latent IS the random
# effect, and under the marginal its distribution is part of the
# likelihood rather than a penalty.
.marg_spec <- function(kind, var, npsi, psi, by, dots, smoothed, c0_set,
                       n_boot_set, label, linear = FALSE) {
  if (!is.null(smoothed)) {
    message(paste("'smoothed' is ignored under marginal = TRUE: the marginal",
                  "likelihood is smooth in every parameter already, and a",
                  "mollifier would only add its bias."))
  }
  if (c0_set) {
    message(paste("'c0' is the scaling schedule of the working",
                  "parametrization, and the marginal construction has none:",
                  "it is ignored."))
  }
  if (n_boot_set) {
    message(paste("'n_boot' is the restart count of the working",
                  "construction; the marginal term starts from the exact",
                  "profile and declares none."))
  }
  if (!is.null(by)) {
    stop(paste("'by' does not combine with marginal = TRUE: the grouping",
               "comes from the break-point's subformula, psi ~ random(~1 |",
               "g), and the term's own coefficients are scalars."),
         call. = FALSE)
  }
  if (!is.numeric(npsi) || length(npsi) != 1L || is.na(npsi) ||
      npsi != round(npsi) || npsi < 1) {
    stop("'npsi' must be a single integer of at least 1.", call. = FALSE)
  }
  npsi <- as.integer(npsi)
  if (kind == "jump" && npsi > 8L) {
    stop(paste("marginal = TRUE covers up to eight break-points: the side",
               "chain's forward recursion costs n K 2^K, which stays cheap,",
               "but a fitting layer reads the posterior over the 2^K side",
               "patterns and evaluates the family once per pattern, which",
               "past 256 components no longer is."), call. = FALSE)
  }
  if (kind != "jump" && npsi > 1L) {
    stop(sprintf(paste("marginal = TRUE covers one break-point for a %s",
                       "term: its conditional is smooth in the positions, so",
                       "several latents need a product quadrature whose",
                       "component count no fitting layer can carry."), kind),
         call. = FALSE)
  }
  if (!is.null(psi) && (!is.numeric(psi) || length(psi) != npsi ||
                        anyNA(psi))) {
    stop("'psi' must give one starting position per break-point.",
         call. = FALSE)
  }
  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !nzchar(label)) {
    stop("'label' must be a single non-empty character string.",
         call. = FALSE)
  }

  need <- paste("marginal = TRUE requires the break-point to carry the",
                "latent, psi ~ random(~1 | g): the marginal likelihood",
                "integrates over it, and without one there is nothing to",
                "integrate.")
  if (!length(dots)) stop(need, call. = FALSE)
  nms <- names(dots)
  sub <- NULL
  for (i in seq_along(dots)) {
    d <- dots[[i]]
    if (!is.null(nms) && nzchar(nms[i])) {
      stop(sprintf("'%s' is not an argument of a break-point term.", nms[i]),
           call. = FALSE)
    }
    if (!inherits(d, "formula") || length(d) != 3L || !is.name(d[[2L]])) {
      stop(paste("arguments in '...' must be two-sided formulas developing",
                 "a coefficient of the term, e.g. psi ~ random(~1 | g)."),
           call. = FALSE)
    }
    lhs <- as.character(d[[2L]])
    if (!lhs %in% c("psi", paste0("psi", seq_len(npsi)))) {
      stop(sprintf(paste("'%s' cannot be developed under marginal = TRUE:",
                         "the term's own coefficients are scalars, estimated",
                         "by maximum likelihood with the prior's",
                         "parameters."), lhs), call. = FALSE)
    }
    if (!identical(lhs, "psi") && npsi > 1L) {
      stop(paste("under marginal = TRUE the latent development is shared:",
                 "name the stem, psi ~ random(~1 | g), so that every",
                 "break-point carries it. A fixed break-point beside a",
                 "random one would put the step back into an ordinary",
                 "parameter."), call. = FALSE)
    }
    if (!is.null(sub)) stop("'psi' is developed twice.", call. = FALSE)
    sub <- d
  }
  if (is.null(sub)) stop(need, call. = FALSE)

  rt <- tryCatch(eval(sub[[3L]], environment(sub)), error = function(e) NULL)
  if (is.null(rt) || !S7::S7_inherits(rt, RandomTerm)) {
    stop(paste("under marginal = TRUE the break-point's subformula must be a",
               "single random(...) call: the latent IS the random effect,",
               "and its distribution is the prior the likelihood integrates",
               "against."), call. = FALSE)
  }
  e <- rt@formula[[2L]]
  if (!(is.numeric(e[[2L]]) && length(e[[2L]]) == 1L && e[[2L]] == 1)) {
    stop(paste("the latent break-point is one value per group: the random()",
               "of a marginal term takes an intercept-only formula,",
               "~ 1 | g."), call. = FALSE)
  }
  prior <- rt@distrib
  if (!is.null(prior)) {
    if (kind != "jump" || npsi > 1L) {
      stop(paste("a prior other than the gaussian is carried by the single",
                 "break-point of a jump term, whose interval masses are",
                 "differences of the prior's cdf; the quadrature",
                 "constructions read the prior's density and are built for",
                 "the gaussian."), call. = FALSE)
    }
    if (!S7::S7_inherits(prior, distributions7::continuous_distrib)) {
      stop("the latent's prior must be a continuous univariate distribution.",
           call. = FALSE)
    }
    ip <- tryCatch(prior@params_interpretation, error = function(e) NULL)
    loc <- names(ip)[unlist(ip) %in% c("location", "mean")]
    if (any(loc %in% prior@params)) {
      stop(sprintf(paste("the prior's location ('%s') must be fixed at zero",
                         "-- fixed(%s, %s = 0) -- since the population",
                         "position m1 carries it."),
                   loc[1L], "student_t1_distrib()", loc[1L]), call. = FALSE)
    }
  }
  if (length(rt@hyper)) {
    stop(paste("'hyper' has nothing to hold under marginal = TRUE: the prior",
               "is part of the likelihood and its parameters are ordinary",
               "ones estimated with the rest."), call. = FALSE)
  }

  MarginalBreakTerm(label = label, kind = kind, var = var, npsi = npsi,
                    linear = linear, group = e[[3L]], prior = prior,
                    spec = list(psi = psi), blueprint = list())
}

# The names of the prior's own parameters for one break-point: the
# gaussian's are the position and the scale, named by the break-point; an
# explicit prior contributes its own free names beside the position.
.marg_prior_names <- function(term, k) {
  if (is.null(term@prior)) c(paste0("m", k), paste0("tau", k))
  else c("m1", term@prior@params)
}

#' @title The Parameters of a Marginal Break-Point Term
#' @name term_params.MarginalBreakTerm
#'
#' @description
#' Numbered, one set per break-point, in a fixed order: the linear effect
#' `beta` where the kind carries one, then the prior's parameters for each
#' break-point, then the changes of slope `gamma1` ... and the changes of
#' level `delta1` ....
#'
#' @details
#' The prior's parameters are `mk` and `tauk` under the default Gaussian, the
#' location and the scale of break-point \eqn{k}. Where `random(distrib = )`
#' named another family the names are that family's own, its location fixed at
#' zero and `mk` carrying the position, so a Student t prior adds `nuk`.
#'
#' Which of `gamma` and `delta` appear is the kind: `"seg"` has the changes of
#' slope, `"jump"` the changes of level, `"jseg"` both. Only `"seg"` and
#' `"jseg"` carry `beta`.
#'
#' @param term A [MarginalBreakTerm()], built or not.
#' @param ... Unused.
#'
#' @return A character vector, of length [term_npar()].
#'
#' @seealso [term_links()] for the chart each rides, [MarginalBreakTerm()] for
#'   what they mean.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(id = rep(1:3, each = 8), x = rep(1:8, 3))
#' dd$y <- rnorm(24, 2 * (dd$x >= 4.5), 0.4)
#'
#' # A step term: the prior's location and scale, and the change of level.
#' term_params(term_build(jump(x, psi ~ random(~ 1 | id), marginal = TRUE), dd))
#'
#' # A continuous one adds the linear effect and the change of slope.
#' term_params(term_build(seg(x, psi ~ random(~ 1 | id), marginal = TRUE), dd))
#'
#' @keywords internal
S7::method(term_params, MarginalBreakTerm) <- function(term, ...) {
  K <- term@npsi
  c(if (term@linear) "beta" else character(0),
    unlist(lapply(seq_len(K), function(k) .marg_prior_names(term, k))),
    if (term@kind %in% c("seg", "jseg")) paste0("gamma", seq_len(K))
    else character(0),
    if (term@kind %in% c("jump", "jseg")) paste0("delta", seq_len(K))
    else character(0))
}

#' @title The Charts of a Marginal Break-Point Term's Parameters
#' @name term_links.MarginalBreakTerm
#'
#' @description
#' The **log** link on every `tauk`, a prior's own link on any parameter it
#' contributes, and the identity on everything else. A prior scale must be
#' positive; a position, a change of level and a change of slope are already
#' unconstrained.
#'
#' @details
#' Where `random(distrib = )` named a family other than the Gaussian, that
#' family's parameters keep the links it declares in `link_params`, so a
#' Student t prior's degrees of freedom ride whatever chart
#' \pkg{distributions7} gives them. Nothing about the chart is restated here.
#'
#' @param term A [MarginalBreakTerm()], built or not.
#' @param ... Unused.
#'
#' @return A named list of \pkg{linkfunctions7} links, one per entry of
#'   [term_params()].
#'
#' @seealso [term_params()], [term_start()].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(id = rep(1:3, each = 8), x = rep(1:8, 3))
#' dd$y <- rnorm(24, 2 * (dd$x >= 4.5), 0.4)
#' tm <- term_build(jump(x, psi ~ random(~ 1 | id), marginal = TRUE), dd)
#' vapply(term_links(tm), function(l) l@link_name, character(1))
#'
#' @keywords internal
S7::method(term_links, MarginalBreakTerm) <- function(term, ...) {
  nm <- term_params(term)
  out <- stats::setNames(vector("list", length(nm)), nm)
  for (p in nm) {
    out[[p]] <- if (grepl("^tau[0-9]+$", p)) {
      linkfunctions7::log_link()
    } else if (!is.null(term@prior) && p %in% term@prior@params) {
      term@prior@link_params[[p]]
    } else {
      linkfunctions7::identity_link()
    }
  }
  out
}

#' @title Build a Marginal Break-Point Term
#' @name term_build.MarginalBreakTerm
#'
#' @description
#' Evaluates the covariate and the grouping variable, works out the interval
#' structure the marginal likelihood is summed or integrated over, and records
#' both in the blueprint together with a data-based starting point for the
#' term's parameters.
#'
#' @details
#' The covariate must vary; a constant one has no break-point to place and
#' throws. The grouping variable must give one value per row.
#'
#' What the blueprint records is the group of each observation, the covariate
#' sorted within each group, and the labels of the groups. The intervals the
#' likelihood decomposes over are the gaps between a group's ordered
#' observations, so they follow from that ordering.
#'
#' The start is computed here because zero is degenerate for this term: with `delta = 0` the intervals are
#' indistinguishable, every mass derivative sums to the derivative of a
#' constant, and the surface is exactly flat in the prior's location and scale.
#' [term_start()] returns what this computed.
#'
#' @param term A [MarginalBreakTerm()].
#' @param data A data frame carrying the covariate and the grouping variable.
#' @param ... Unused.
#'
#' @return The term with `blueprint` filled. [term_is_built()] stays `FALSE`,
#'   that predicate testing for a design block.
#'
#' @seealso [jump()], [term_loglik()], [term_start()].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(id = rep(1:3, each = 8), x = rep(1:8, 3))
#' dd$y <- rnorm(24, 2 * (dd$x >= 4.5), 0.4)
#'
#' b <- term_build(jump(x, psi ~ random(~ 1 | id), marginal = TRUE), dd)
#' names(b@blueprint)
#' b@blueprint$labels
#'
#' # The start is data-based, zero being degenerate here.
#' round(term_start(b), 4)
#'
#' # A constant covariate has no break-point to place.
#' flat <- transform(dd, x = 1)
#' try(term_build(jump(x, psi ~ random(~ 1 | id), marginal = TRUE), flat))
#'
#' @keywords internal
S7::method(term_build, MarginalBreakTerm) <- function(term, data, ...) {
  xv <- .seg_x(term@var, data)
  n <- length(xv)
  if (diff(range(xv)) <= 0) {
    stop("the covariate of a break-point term must vary.", call. = FALSE)
  }
  g <- eval(term@group, data, baseenv())
  if (length(g) != n) {
    stop(paste("the grouping of a marginal break-point term must give one",
               "value per row."), call. = FALSE)
  }
  if (anyNA(g)) {
    stop(paste("the grouping of a marginal break-point term must not contain",
               "missing values."), call. = FALSE)
  }
  f <- factor(g)
  rows <- split(seq_len(n), f)
  term@blueprint <- list(
    n = n, x = xv,
    # each group's rows sorted by the covariate: the intervals between them
    # are where the sum (or the quadrature) runs
    groups = lapply(rows, function(r) r[order(xv[r])]),
    labels = names(rows))
  term
}

# log(Phi(zu) - Phi(zl)), stable in either tail: the difference is taken on
# whichever side keeps both terms away from one, through the log-scale tail
# probabilities.
.marg_lmass <- function(zl, zu) {
  out <- numeric(length(zl))
  lo_inf <- !is.finite(zl) & zl < 0
  hi_inf <- !is.finite(zu) & zu > 0
  i <- lo_inf & hi_inf
  out[i] <- 0
  i <- lo_inf & !hi_inf
  out[i] <- stats::pnorm(zu[i], log.p = TRUE)
  i <- hi_inf & !lo_inf
  out[i] <- stats::pnorm(zl[i], lower.tail = FALSE, log.p = TRUE)
  i <- !lo_inf & !hi_inf
  if (any(i)) {
    a <- zl[i]
    b <- zu[i]
    lower <- (a + b) < 0
    la <- ifelse(lower, stats::pnorm(b, log.p = TRUE),
                 stats::pnorm(a, lower.tail = FALSE, log.p = TRUE))
    lb <- ifelse(lower, stats::pnorm(a, log.p = TRUE),
                 stats::pnorm(b, lower.tail = FALSE, log.p = TRUE))
    d <- lb - la
    out[i] <- la + ifelse(d > -log(2), log(-expm1(d)), log1p(-exp(d)))
  }
  out
}

.marg_lse <- function(a) {
  mx <- max(a)
  if (!is.finite(mx)) return(-Inf)
  mx + log(sum(exp(a - mx)))
}

# The interval decomposition of one group under the GAUSSIAN prior:
# boundaries between its sorted covariate values, the log masses and the
# derivatives of those masses in (m, tau) on the parameter scale. Interval
# j (0-based) is (x_(j), x_(j+1)] with x_(0) = -Inf and x_(n+1) = +Inf: a
# break-point in interval j leaves the observations of sorted rank above j
# on the shifted side. An interval whose mass underflows carries zero
# posterior weight, and its derivative rows are set to zero so that the
# products they enter stay finite.
.marg_intervals <- function(xs, m, tau, d2 = FALSE) {
  b <- c(-Inf, xs, Inf)
  nb <- length(b)
  zl <- (b[-nb] - m) / tau
  zu <- (b[-1L] - m) / tau
  lm <- .marg_lmass(zl, zu)
  mass <- exp(lm)
  ok <- mass > 0
  fd <- function(z) {
    phi <- ifelse(is.finite(z), stats::dnorm(z), 0)
    zf <- ifelse(is.finite(z), z, 0)
    list(phi = phi, zp = zf * phi, z2p = zf^2 * phi, z3p = zf^3 * phi)
  }
  L <- fd(zl)
  U <- fd(zu)
  # dF/dm = -phi(z)/tau and dF/dtau = -z phi(z)/tau at each boundary; a mass
  # derivative is the difference of the two boundary values
  dm <- ifelse(ok, ((L$phi - U$phi) / tau) / mass, 0)
  dt <- ifelse(ok, ((L$zp - U$zp) / tau) / mass, 0)
  out <- list(lm = lm, mass = mass, ok = ok, dm = dm, dt = dt,
              zl = zl, zu = zu, phi_l = L$phi, phi_u = U$phi,
              zp_l = L$zp, zp_u = U$zp)
  if (d2) {
    # second derivatives of F = Phi((b - m)/tau):
    #   F_mm = -z phi/tau^2, F_mtau = phi(1 - z^2)/tau^2,
    #   F_tautau = z phi (2 - z^2)/tau^2
    mm <- (L$zp - U$zp) / tau^2
    mt <- ((U$phi - U$z2p) - (L$phi - L$z2p)) / tau^2
    tt <- ((2 * U$zp - U$z3p) - (2 * L$zp - L$z3p)) / tau^2
    out$dmm <- ifelse(ok, mm / mass, 0) - dm * dm
    out$dmt <- ifelse(ok, mt / mass, 0) - dm * dt
    out$dtt <- ifelse(ok, tt / mass, 0) - dt * dt
  }
  out
}

# The same decomposition under an EXPLICIT prior (one break-point): masses
# are differences of the prior's cdf at the boundaries shifted by the
# position, and their derivatives come from the density (for the position)
# and from the cdf surface (for the prior's own parameters), the route the
# censored likelihoods use. Columns of dlm follow .marg_prior_names().
.marg_prior_intervals <- function(prior, xs, m, theta) {
  b <- c(-Inf, xs, Inf)
  nb <- length(b)
  l <- b[-nb]
  u <- b[-1L]
  th <- stats::setNames(as.list(theta[prior@params]), prior@params)
  cdf_at <- function(q) {
    out <- numeric(length(q))
    fin <- is.finite(q)
    if (any(fin)) {
      out[fin] <- as.numeric(distributions7::distrib_cdf(prior, q[fin] - m,
                                                         th))
    }
    out[!fin] <- ifelse(q[!fin] > 0, 1, 0)
    out
  }
  pdf_at <- function(q) {
    out <- numeric(length(q))
    fin <- is.finite(q)
    if (any(fin)) {
      out[fin] <- as.numeric(distributions7::distrib_pdf(prior, q[fin] - m,
                                                         th))
    }
    out
  }
  gcdf_at <- function(q) {
    out <- matrix(0, length(q), length(prior@params))
    fin <- is.finite(q)
    if (any(fin)) {
      g <- distributions7::distrib_grad_cdf(prior, q[fin] - m, th,
                                            lower.tail = TRUE, log = FALSE)
      for (j in seq_along(prior@params)) {
        out[fin, j] <- rep_len(as.numeric(g[[prior@params[j]]]), sum(fin))
      }
    }
    out
  }
  Fl <- cdf_at(l)
  Fu <- cdf_at(u)
  mass <- pmax(Fu - Fl, 0)
  ok <- mass > 0
  lm <- ifelse(ok, log(mass), -Inf)
  # d/dm F(b - m) = -f(b - m), so the position's column is the density
  # difference; the prior's own columns are the cdf-gradient differences
  dm_col <- ifelse(ok, (pdf_at(l) - pdf_at(u)) / mass, 0)
  Gd <- gcdf_at(u) - gcdf_at(l)
  Gd[!ok, ] <- 0
  Gd[ok, ] <- Gd[ok, , drop = FALSE] / mass[ok]
  dlm <- cbind(dm_col, Gd)
  colnames(dlm) <- c("m1", prior@params)
  list(lm = lm, mass = mass, ok = ok, dlm = dlm, l = l, u = u)
}

.marg_check_psi <- function(term, psi) {
  nm <- term_params(term)
  v <- unlist(psi[nm])
  if (length(v) != length(nm) || anyNA(v)) {
    stop(sprintf("'psi' must supply %s.", paste(nm, collapse = ", ")),
         call. = FALSE)
  }
  taus <- grep("^tau[0-9]+$", nm, value = TRUE)
  if (length(taus) && any(v[taus] <= 0)) {
    stop("every prior scale must be positive.", call. = FALSE)
  }
  v
}

.marg_built <- function(term) {
  bp <- term@blueprint
  if (!length(bp)) {
    stop("the term has not been built; call term_build(term, data) first.",
         call. = FALSE)
  }
  bp
}

# ---------------------------------------------------------------------------
# The step kind. With K latent positions the conditional is constant on the
# product partition of each coordinate's intervals -- (n+1)^K cells per
# group -- but the sum is never taken over the cells: the side process
# S_t = {k : psi_k <= x_(t)} over the SORTED observations is monotone on
# the subset lattice with independent coordinates, so it is a hidden
# Markov chain on the 2^K side patterns whose transition factors over the
# coordinates, and the forward recursion costs n K 2^K where the cell sum
# costs (n+1)^K. A coordinate that flips at step t contributes its
# interval's prior mass as the transition weight; one that never flips
# contributes its upper-tail mass at the end, through the survival factor
# of the final state.

# the per-coordinate interval masses with their derivatives on the
# parameter scale, and the survivals the per-step totals read -- the
# reverse cumulated masses, so the telescoping is exact by construction.
# d2 adds the second derivatives, for the propagated Hessian.
.marg_hmm_masses <- function(term, xs, v, d2 = FALSE) {
  K <- term@npsi
  ng <- length(xs)
  J <- ng + 1L
  pr <- .marg_jump_prior(term, xs, v)
  out <- list(pcols = lapply(pr$dlm, colnames), ivs = pr$ivs, J = J)
  out$M <- vapply(pr$lm, exp, numeric(J))
  out$dM <- vector("list", K)
  out$d2M <- if (d2) vector("list", K) else NULL
  for (k in seq_len(K)) {
    mass <- out$M[, k]
    out$dM[[k]] <- mass * pr$dlm[[k]]
    if (d2) {
      npc <- ncol(pr$dlm[[k]])
      A <- array(0, c(J, npc, npc))
      if (is.null(term@prior)) {
        iv <- .marg_intervals(xs, v[[paste0("m", k)]],
                              v[[paste0("tau", k)]], d2 = TRUE)
        # d2 mass = mass (d2 log mass + d log mass squared)
        dl <- pr$dlm[[k]]
        A[, 1L, 1L] <- mass * (iv$dmm + dl[, 1L]^2)
        A[, 1L, 2L] <- mass * (iv$dmt + dl[, 1L] * dl[, 2L])
        A[, 2L, 1L] <- A[, 1L, 2L]
        A[, 2L, 2L] <- mass * (iv$dtt + dl[, 2L]^2)
      } else {
        A <- .marg_prior_d2mass(term@prior, xs, v)
      }
      out$d2M[[k]] <- A
    }
  }
  # survivals over the steps: SV[t + 1, k] is the mass beyond x_(t), the
  # reverse cumulated masses, with SV[1, ] = 1
  out$SV <- apply(out$M, 2L, function(m) rev(cumsum(rev(m))))[-1L, ,
                                                              drop = FALSE]
  out$SV <- rbind(1, out$SV)
  out$dSV <- lapply(seq_len(K), function(k) {
    D <- apply(out$dM[[k]], 2L, function(m) rev(cumsum(rev(m))))
    if (is.null(dim(D))) D <- matrix(D, nrow = J)
    rbind(0, D[-1L, , drop = FALSE])
  })
  if (d2) {
    out$d2SV <- lapply(seq_len(K), function(k) {
      A <- out$d2M[[k]]
      npc <- dim(A)[2L]
      B <- array(0, c(J, npc, npc))
      for (i in seq_len(npc)) {
        for (j in seq_len(npc)) {
          rc <- rev(cumsum(rev(A[, i, j])))
          B[, i, j] <- c(0, rc[-1L])
        }
      }
      B
    })
  }
  out
}

# the second derivatives of an explicit prior's interval masses, from the
# surfaces the censored likelihoods use: the position's through the
# density's own response derivative, the mixed block through the density
# times the score, and the prior's own through the cdf Hessian
.marg_prior_d2mass <- function(prior, xs, v) {
  b <- c(-Inf, xs, Inf)
  nb <- length(b)
  l <- b[-nb]
  u <- b[-1L]
  m <- v[["m1"]]
  th <- stats::setNames(as.list(v[prior@params]), prior@params)
  npc <- 1L + length(prior@params)
  J <- length(l)
  np <- length(prior@params)
  pdfgy_at <- function(q) {
    out <- numeric(length(q))
    fin <- is.finite(q)
    if (any(fin)) {
      qq <- q[fin] - m
      pd <- as.numeric(distributions7::distrib_pdf(prior, qq, th))
      gy <- as.numeric(distributions7::distrib_grad_y(prior, qq, th))
      out[fin] <- pd * gy
    }
    out
  }
  pdfg_at <- function(q) {
    out <- matrix(0, length(q), np)
    fin <- is.finite(q)
    if (any(fin)) {
      qq <- q[fin] - m
      pd <- as.numeric(distributions7::distrib_pdf(prior, qq, th))
      g <- distributions7::distrib_gradient(prior, qq, th)
      for (j in seq_len(np)) {
        out[fin, j] <- pd * rep_len(as.numeric(g[[prior@params[j]]]),
                                    sum(fin))
      }
    }
    out
  }
  hcdf_at <- function(q) {
    out <- matrix(0, length(q), np * np)
    fin <- is.finite(q)
    if (any(fin)) {
      qq <- q[fin] - m
      h <- distributions7::distrib_hess_cdf(prior, qq, th,
                                            lower.tail = TRUE, log = FALSE)
      for (i in seq_len(np)) {
        for (j in seq_len(np)) {
          key <- paste(prior@params[sort(c(i, j))], collapse = "_")
          out[fin, (i - 1L) * np + j] <- rep_len(as.numeric(h[[key]]),
                                                 sum(fin))
        }
      }
    }
    out
  }
  A <- array(0, c(J, npc, npc))
  # d2/dm2 through the density's response derivative at the endpoints
  A[, 1L, 1L] <- -pdfgy_at(l) + pdfgy_at(u)
  # mixed: the density times the prior's own score at the endpoints
  fg_l <- pdfg_at(l)
  fg_u <- pdfg_at(u)
  for (j in seq_len(np)) {
    A[, 1L, 1L + j] <- fg_l[, j] - fg_u[, j]
    A[, 1L + j, 1L] <- A[, 1L, 1L + j]
  }
  # the prior's own block: cdf Hessian differences
  H_l <- hcdf_at(l)
  H_u <- hcdf_at(u)
  for (i in seq_len(np)) {
    for (j in seq_len(np)) {
      A[, 1L + i, 1L + j] <- H_u[, (i - 1L) * np + j] -
        H_l[, (i - 1L) * np + j]
    }
  }
  A
}

# the transition index sets of the side chain: for each coordinate, the
# states without its bit and their partners with it
.marg_hmm_idx <- function(K) {
  bits <- .marg_bits(K, numeric(K))$bits
  lapply(seq_len(K), function(k) {
    i0 <- which(!bits[, k])
    list(i0 = i0, i1 = i0 + 2^(k - 1L))
  })
}

# the survival product of every state at one step: the product over the
# INACTIVE coordinates of their survivals
.marg_hmm_sp <- function(idx, SVt, P) {
  sp <- rep(1, P)
  for (k in seq_along(idx)) sp[idx[[k]]$i0] <- sp[idx[[k]]$i0] * SVt[k]
  sp
}

# the per-coordinate interval structures and the prior-parameter columns
# they own, shared by every method of the step kind
.marg_jump_prior <- function(term, xs, v) {
  K <- term@npsi
  if (is.null(term@prior)) {
    ivs <- lapply(seq_len(K), function(k)
      .marg_intervals(xs, v[[paste0("m", k)]], v[[paste0("tau", k)]]))
    list(
      lm = lapply(ivs, `[[`, "lm"),
      dlm = lapply(seq_len(K), function(k) {
        M <- cbind(ivs[[k]]$dm, ivs[[k]]$dt)
        M[!ivs[[k]]$ok, ] <- 0
        colnames(M) <- c(paste0("m", k), paste0("tau", k))
        M
      }),
      ivs = ivs)
  } else {
    pv <- .marg_prior_intervals(term@prior, xs, v[["m1"]],
                                v[term@prior@params])
    list(lm = list(pv$lm), dlm = list(pv$dlm), ivs = list(pv))
  }
}

# the 2^K side patterns: their bit matrix and the shift each adds
.marg_bits <- function(K, deltas) {
  P <- 2^K
  bits <- matrix(0L, P, K)
  for (k in seq_len(K)) bits[, k] <- bitwAnd(seq_len(P) - 1L, 2^(k - 1L)) > 0
  list(bits = bits, shifts = as.numeric(bits %*% deltas))
}

#' @title Log-Likelihood of a Marginal Break-Point Term
#' @name term_loglik.MarginalBreakTerm
#' @description
#' The exact marginal likelihood of a break-point model with latent
#' positions per group, with the exact derivatives in the term's own
#' parameters on the parameter scale.
#' @details
#' For the step kind the conditional is constant on the product partition
#' of the intervals between a group's ordered observations, and the exact
#' sum over it is taken by the forward recursion of the side chain: the
#' monotone process of active break-points is a hidden Markov chain on the
#' \eqn{2^K} side patterns whose transition factors over the coordinates,
#' each flip weighted by its interval's prior mass, so the cost is
#' \eqn{n K 2^K}, against the \eqn{(n+1)^K} of the cells. The
#' derivatives ride the same recursion: the masses' in the prior's
#' parameters, the emissions' in the changes of level.
#' For the continuous kinds the conditional is smooth within an interval
#' and the integral runs on a fixed Gauss-Kronrod panel per interval
#' ([numericals7::gauss_kronrod15()]), the interior nodes fixed
#' points of the data so that the derivatives in the prior's parameters
#' read the prior alone; the region below the data, where the hinge keeps
#' moving, is covered by panels that follow the prior's bulk, whose node
#' motion the derivatives carry, and the region above it, where the
#' conditional is constant, contributes its closed tail mass.
#'
#' The per-observation contributions are the one-step predictive densities
#' given the group's observations at smaller covariate values, which sum
#' to each group's marginal log-likelihood; producing them costs one pass
#' over the component weights per observation, so the decomposition is a
#' factor of the group's size dearer than the total alone, which the
#' per-group cell or node count already prices.
#' @param term A built [MarginalBreakTerm()].
#' @param eta The static predictor.
#' @param y The response, reaching the sum through the callbacks.
#' @param logdens,score The log-density and its derivative in the predictor.
#' @param psi The parameters, named as [term_params()].
#' @param ... Unused.
#' @return A list with `loglik` and `jacobian`, the latter on the
#'   parameter scale.
#' @keywords internal
S7::method(term_loglik, MarginalBreakTerm) <- function(term, eta, y, logdens,
                                                       score, psi, ...) {
  if (term@kind == "jump") {
    .marg_jump_loglik(term, eta, y, logdens, score, psi)
  } else {
    .marg_seg_loglik(term, eta, y, logdens, score, psi)
  }
}

.marg_jump_loglik <- function(term, eta, y, logdens, score, psi) {
  bp <- .marg_built(term)
  v <- .marg_check_psi(term, psi)
  nm <- term_params(term)
  K <- term@npsi
  n <- bp$n
  idx <- seq_len(n)
  pb <- .marg_bits(K, v[paste0("delta", seq_len(K))])
  P <- 2^K
  LD <- matrix(0, n, P)
  SC <- matrix(0, n, P)
  for (p in seq_len(P)) {
    LD[, p] <- as.numeric(logdens(eta + pb$shifts[p], idx))
    SC[, p] <- as.numeric(score(eta + pb$shifts[p], idx))
  }
  tix <- .marg_hmm_idx(K)
  dcol <- match(paste0("delta", seq_len(K)), nm)

  ll <- numeric(n)
  jac <- matrix(0, n, length(nm), dimnames = list(NULL, nm))
  for (rs in bp$groups) {
    ng <- length(rs)
    hm <- .marg_hmm_masses(term, bp$x[rs], v)
    pcol <- lapply(hm$pcols, function(pc) match(pc, nm))

    alpha <- numeric(P)
    alpha[1L] <- 1
    Dal <- matrix(0, P, length(nm))
    ls <- 0
    lnum_prev <- 0
    rat_prev <- numeric(length(nm))
    for (t in seq_len(ng)) {
      row <- rs[t]
      # the transition: each coordinate's flip moves mass onto its bit,
      # weighted by the interval's prior mass; the factors commute, so a
      # sequential in-place application is the product
      for (k in seq_len(K)) {
        i0 <- tix[[k]]$i0
        i1 <- tix[[k]]$i1
        q <- hm$M[t, k]
        Dal[i1, ] <- Dal[i1, ] + q * Dal[i0, ]
        Dal[i1, pcol[[k]]] <- Dal[i1, pcol[[k]]] +
          alpha[i0] %o% hm$dM[[k]][t, ]
        alpha[i1] <- alpha[i1] + q * alpha[i0]
      }
      # the emission: the state IS the side pattern, so the per-pattern
      # log-density indexes directly
      lf <- LD[row, ]
      mx <- max(lf)
      w <- exp(lf - mx)
      ls <- ls + mx
      alpha <- alpha * w
      Dal <- Dal * w
      for (k in seq_len(K)) {
        Dal[, dcol[k]] <- Dal[, dcol[k]] + alpha * SC[row, ] * pb$bits[, k]
      }
      # the per-step total: every inactive coordinate carries its survival
      sp <- .marg_hmm_sp(tix, hm$SV[t + 1L, ], P)
      num <- sum(alpha * sp)
      dnum <- as.numeric(crossprod(Dal, sp))
      for (k in seq_len(K)) {
        spk <- rep(1, P)
        for (k2 in seq_len(K)) {
          if (k2 == k) next
          spk[tix[[k2]]$i0] <- spk[tix[[k2]]$i0] * hm$SV[t + 1L, k2]
        }
        a0 <- sum(alpha[tix[[k]]$i0] * spk[tix[[k]]$i0])
        dnum[pcol[[k]]] <- dnum[pcol[[k]]] + a0 * hm$dSV[[k]][t + 1L, ]
      }
      lnum <- log(num) + ls
      ll[row] <- lnum - lnum_prev
      jac[row, ] <- dnum / num - rat_prev
      lnum_prev <- lnum
      rat_prev <- dnum / num
      # a constant rescale keeps the state representable and moves no
      # derivative
      cs <- max(alpha)
      if (cs > 0 && (cs > 1e100 || cs < 1e-100)) {
        alpha <- alpha / cs
        Dal <- Dal / cs
        ls <- ls + log(cs)
      }
    }
  }
  list(loglik = ll, jacobian = jac)
}

# ---------------------------------------------------------------------------
# The continuous kinds: nodes. The conditional is smooth within an interval,
# so each interval carries a fixed Gauss-Kronrod panel; the interval
# boundaries are the data, so interior nodes do not move with the prior's
# parameters. Below the data the hinge keeps varying, and the panels there
# follow the prior's bulk: their nodes move with (m, tau) and the
# derivatives carry that motion. Above the data every device is saturated
# at zero, the conditional is constant, and the tail contributes its
# closed mass through a node at +Inf.

# the node set of one group: positions, log quadrature weights (prior
# density included), and the geometry the derivatives need
.marg_nodes <- function(xs, m, tau) {
  gk <- numericals7::gauss_kronrod15()
  x1 <- xs[1L]
  xn <- xs[length(xs)]
  lo_pts <- numeric(0)
  hi_pts <- numeric(0)
  edge <- logical(0)
  # the lower region, from 8.5 prior sds below the bulk up to the first
  # observation; its panel edges are affine in the lower limit, which is
  # what the node-motion derivatives read
  a <- min(m, x1) - 8.5 * tau
  wid <- x1 - a
  ns <- min(8L, max(4L, ceiling(wid / (2.5 * tau))))
  ed <- seq(a, x1, length.out = ns + 1L)
  for (s in seq_len(ns)) {
    lo_pts <- c(lo_pts, ed[s])
    hi_pts <- c(hi_pts, ed[s + 1L])
    edge <- c(edge, TRUE)
  }
  # the interior intervals, split so that no panel spans more than a couple
  # of prior sds
  for (j in seq_len(length(xs) - 1L)) {
    wj <- xs[j + 1L] - xs[j]
    if (wj <= 0) next
    nsj <- min(6L, max(1L, ceiling(wj / (2.5 * tau))))
    edj <- seq(xs[j], xs[j + 1L], length.out = nsj + 1L)
    for (s in seq_len(nsj)) {
      lo_pts <- c(lo_pts, edj[s])
      hi_pts <- c(hi_pts, edj[s + 1L])
      edge <- c(edge, FALSE)
    }
  }
  nn <- length(gk$nodes)
  np <- length(lo_pts)
  h <- rep((hi_pts - lo_pts) / 2, each = nn)
  mid <- rep((hi_pts + lo_pts) / 2, each = nn)
  p <- mid + h * gk$nodes
  wk <- rep(gk$wk, np)
  z <- (p - m) / tau
  lw <- log(wk * h) + stats::dnorm(z, log = TRUE) - log(tau)
  is_edge <- rep(edge, each = nn)
  # the relative position of a node inside the lower region: its motion is
  # (1 - t) times the motion of the lower limit
  tfrac <- ifelse(is_edge, (p - a) / pmax(x1 - a, .Machine$double.eps), 0)
  # the closed upper tail: every device saturates at zero above the data,
  # so the conditional there is the one a node at +Inf evaluates
  zn <- (xn - m) / tau
  mr <- numericals7::mills_ratio(-zn)$r
  p <- c(p, Inf)
  lw <- c(lw, stats::pnorm(zn, lower.tail = FALSE, log.p = TRUE))
  z <- c(z, Inf)
  is_edge <- c(is_edge, FALSE)
  tfrac <- c(tfrac, 0)
  # d log(weight)/d(m, tau) at fixed node: the prior's own derivatives, the
  # tail's through the Mills ratio; the log panel width of the lower region
  # scales with the limit as well
  glw_m <- c(z[-length(p)] / tau, mr / tau)
  glw_t <- c((z[-length(p)]^2 - 1) / tau, zn * mr / tau)
  da_m <- as.numeric(m < x1)
  da_t <- -8.5
  glw_m <- glw_m + is_edge * (-da_m / pmax(x1 - a, .Machine$double.eps))
  glw_t <- glw_t + is_edge * (-da_t / pmax(x1 - a, .Machine$double.eps))
  dpsi_m <- ifelse(is_edge, (1 - tfrac) * da_m, 0)
  dpsi_t <- ifelse(is_edge, (1 - tfrac) * da_t, 0)
  # the node-motion part of the prior weight's derivative: d log phi/d psi
  dlphi_dpsi <- ifelse(is.finite(z), -z / tau, 0)
  # The second derivatives of the log weight in (m, tau), node motion
  # included: the motion is affine, so the chain rule closes at
  # z = (psi(theta) - m)/tau and its own first two derivatives, plus the
  # panel widths of the lower region, which scale with the moving limit,
  # and the closed tail through the Mills ratio.
  fin <- is.finite(z)
  zf <- ifelse(fin, z, 0)
  dzm <- (dpsi_m - 1) / tau
  dzt <- (dpsi_t - zf) / tau
  d2zmt <- -(dpsi_m - 1) / tau^2
  d2ztt <- -2 * (dpsi_t - zf) / tau^2
  alw_mm <- -dzm^2
  alw_mt <- -dzm * dzt - zf * d2zmt
  alw_tt <- -dzt^2 - zf * d2ztt + 1 / tau^2
  denom <- pmax(x1 - a, .Machine$double.eps)
  alw_mm <- alw_mm - is_edge * (da_m * da_m) / denom^2
  alw_mt <- alw_mt - is_edge * (da_m * da_t) / denom^2
  alw_tt <- alw_tt - is_edge * (da_t * da_t) / denom^2
  i_t <- length(p)
  lamp <- mr * (mr - zn)
  alw_mm[i_t] <- -lamp / tau^2
  alw_mt[i_t] <- -lamp * zn / tau^2 - mr / tau^2
  alw_tt[i_t] <- -lamp * zn^2 / tau^2 - 2 * mr * zn / tau^2
  list(p = p, lw = lw, z = z,
       glw_m = glw_m + dlphi_dpsi * dpsi_m,
       glw_t = glw_t + dlphi_dpsi * dpsi_t,
       dpsi_m = dpsi_m, dpsi_t = dpsi_t,
       alw_mm = alw_mm, alw_mt = alw_mt, alw_tt = alw_tt)
}

# the shift each node adds to each of the group's observations, and its
# derivatives in the term's own coefficients and in the position
.marg_seg_shift <- function(term, xg, nd, v) {
  hinge <- outer(xg, nd$p, function(x, p) pmax(x - ifelse(is.finite(p), p,
                                                          Inf), 0))
  step <- outer(xg, nd$p, function(x, p) as.numeric(is.finite(p) & x >= p))
  shift <- 0
  if (term@linear) shift <- shift + v[["beta"]] * xg
  shift <- shift + v[["gamma1"]] * hinge
  if (term@kind == "jseg") shift <- shift + v[["delta1"]] * step
  # d shift/d psi = -gamma 1(x > psi); the step's derivative is zero away
  # from the point
  dshift_dpsi <- -v[["gamma1"]] * (hinge > 0)
  list(shift = shift, hinge = hinge, step = step, dshift_dpsi = dshift_dpsi)
}

.marg_seg_loglik <- function(term, eta, y, logdens, score, psi) {
  bp <- .marg_built(term)
  v <- .marg_check_psi(term, psi)
  nm <- term_params(term)
  n <- bp$n
  ll <- numeric(n)
  jac <- matrix(0, n, length(nm), dimnames = list(NULL, nm))
  for (rs in bp$groups) {
    ng <- length(rs)
    nd <- .marg_nodes(bp$x[rs], v[["m1"]], v[["tau1"]])
    C <- length(nd$p)
    sh <- .marg_seg_shift(term, bp$x[rs], nd, v)
    ei <- rep(eta[rs], C) + as.numeric(sh$shift)
    ii <- rep(rs, C)
    LD <- matrix(as.numeric(logdens(ei, ii)), ng, C)
    SC <- matrix(as.numeric(score(ei, ii)), ng, C)

    A <- nd$lw
    tot <- .marg_lse(A)
    u <- exp(A - tot)
    acc <- list()
    own <- c(if (term@linear) "beta", "gamma1",
             if (term@kind == "jseg") "delta1")
    dmat <- list(beta = if (term@linear) matrix(bp$x[rs], ng, C) else NULL,
                 gamma1 = sh$hinge,
                 delta1 = if (term@kind == "jseg") sh$step else NULL)
    for (p in own) acc[[p]] <- numeric(C)
    accp <- numeric(C)   # the psi-motion accumulator, sum_t s * dshift/dpsi
    for (t in seq_len(ng)) {
      row <- rs[t]
      A2 <- A + LD[t, ]
      tot2 <- .marg_lse(A2)
      w <- exp(A2 - tot2)
      ll[row] <- tot2 - tot
      for (p in own) {
        a2 <- acc[[p]] + SC[t, ] * dmat[[p]][t, ]
        jac[row, p] <- sum(w * a2) - sum(u * acc[[p]])
        acc[[p]] <- a2
      }
      ap2 <- accp + SC[t, ] * sh$dshift_dpsi[t, ]
      jac[row, "m1"] <- sum(w * (nd$glw_m + ap2 * nd$dpsi_m)) -
        sum(u * (nd$glw_m + accp * nd$dpsi_m))
      jac[row, "tau1"] <- sum(w * (nd$glw_t + ap2 * nd$dpsi_t)) -
        sum(u * (nd$glw_t + accp * nd$dpsi_t))
      accp <- ap2
      A <- A2
      tot <- tot2
      u <- w
    }
  }
  list(loglik = ll, jacobian = jac)
}

#' @title Posterior Components of a Marginal Break-Point Term
#' @name term_posterior.MarginalBreakTerm
#' @description
#' The component weights Fisher's identity takes at every observation: for
#' the step kind, the posterior probability of each side pattern given the
#' whole sample, one column per pattern of active break-points; for the
#' continuous kinds, each group's posterior over its quadrature nodes,
#' repeated down the group's rows and zero-padded to the widest group.
#' @param term A built [MarginalBreakTerm()].
#' @param eta The static predictor.
#' @param y The response.
#' @param logdens The log-density.
#' @param psi The term's parameters.
#' @param ... Unused.
#' @return A numeric matrix, one row per observation, rows summing to one.
#' @keywords internal
S7::method(term_posterior, MarginalBreakTerm) <- function(term, eta, y,
                                                          logdens, psi, ...) {
  if (term@kind == "jump") {
    .marg_jump_posterior(term, eta, y, logdens, psi)
  } else {
    .marg_seg_posterior(term, eta, y, logdens, psi)$gamma
  }
}

# the plain forward of one group's side chain, storing the state after
# every step; the per-step rescales are constants and live in a log-scale
# vector
.marg_hmm_forward <- function(hm, LD, rs, tix, K) {
  ng <- length(rs)
  P <- 2^K
  A <- matrix(0, ng, P)
  lsA <- numeric(ng)
  alpha <- numeric(P)
  alpha[1L] <- 1
  ls <- 0
  for (t in seq_len(ng)) {
    for (k in seq_len(K)) {
      i0 <- tix[[k]]$i0
      i1 <- tix[[k]]$i1
      alpha[i1] <- alpha[i1] + hm$M[t, k] * alpha[i0]
    }
    lf <- LD[rs[t], ]
    mx <- max(lf)
    alpha <- alpha * exp(lf - mx)
    ls <- ls + mx
    cs <- max(alpha)
    if (cs > 0) {
      alpha <- alpha / cs
      ls <- ls + log(cs)
    }
    A[t, ] <- alpha
    lsA[t] <- ls
  }
  list(A = A, lsA = lsA)
}

# the matching backward: beta_t is defined so that the likelihood is
# sum_S alpha_t(S) beta_t(S) at every t, the final one the survival product
.marg_hmm_backward <- function(hm, LD, rs, tix, K) {
  ng <- length(rs)
  P <- 2^K
  B <- matrix(0, ng, P)
  lsB <- numeric(ng)
  beta <- .marg_hmm_sp(tix, hm$SV[ng + 1L, ], P)
  ls <- 0
  B[ng, ] <- beta
  lsB[ng] <- 0
  if (ng > 1L) {
    for (t in seq.int(ng, 2L)) {
      lf <- LD[rs[t], ]
      mx <- max(lf)
      tmp <- beta * exp(lf - mx)
      ls <- ls + mx
      for (k in seq_len(K)) {
        i0 <- tix[[k]]$i0
        i1 <- tix[[k]]$i1
        tmp[i0] <- tmp[i0] + hm$M[t, k] * tmp[i1]
      }
      cs <- max(tmp)
      if (cs > 0) {
        tmp <- tmp / cs
        ls <- ls + log(cs)
      }
      beta <- tmp
      B[t - 1L, ] <- beta
      lsB[t - 1L] <- ls
    }
  }
  list(B = B, lsB = lsB)
}

.marg_jump_posterior <- function(term, eta, y, logdens, psi) {
  bp <- .marg_built(term)
  v <- .marg_check_psi(term, psi)
  K <- term@npsi
  P <- 2^K
  n <- bp$n
  idx <- seq_len(n)
  pb <- .marg_bits(K, v[paste0("delta", seq_len(K))])
  LD <- matrix(0, n, P)
  for (p in seq_len(P)) {
    LD[, p] <- as.numeric(logdens(eta + pb$shifts[p], idx))
  }
  tix <- .marg_hmm_idx(K)
  out <- matrix(0, n, P)
  for (rs in bp$groups) {
    ng <- length(rs)
    hm <- .marg_hmm_masses(term, bp$x[rs], v)
    fw <- .marg_hmm_forward(hm, LD, rs, tix, K)
    bw <- .marg_hmm_backward(hm, LD, rs, tix, K)
    # the emission at step t already sits in alpha_t, and beta_t carries
    # everything after it, so the product is the state posterior; the
    # observation's side pattern IS the state, and the scales cancel in
    # the normalization
    for (t in seq_len(ng)) {
      w <- fw$A[t, ] * bw$B[t, ]
      out[rs[t], ] <- w / sum(w)
    }
  }
  out
}

# the posterior of each coordinate's flip interval: the probability, given
# the data, that the latent sits between two adjacent observations, read
# off the forward-backward pair -- the numerator forces the coordinate to
# flip at one step and lets everything else run
.marg_hmm_flip <- function(term, rs, hm, LD, tix, K) {
  ng <- length(rs)
  P <- 2^K
  J <- ng + 1L
  fw <- .marg_hmm_forward(hm, LD, rs, tix, K)
  bw <- .marg_hmm_backward(hm, LD, rs, tix, K)
  logL <- log(sum(fw$A[ng, ] * bw$B[ng, ])) + fw$lsA[ng] + bw$lsB[ng]
  wm <- matrix(0, J, K)
  for (t in seq_len(ng)) {
    a_prev <- if (t == 1L) {
      z <- numeric(P)
      z[1L] <- 1
      z
    } else fw$A[t - 1L, ]
    ls_prev <- if (t == 1L) 0 else fw$lsA[t - 1L]
    lf <- LD[rs[t], ]
    mx <- max(lf)
    w_t <- exp(lf - mx)
    for (k in seq_len(K)) {
      af <- numeric(P)
      af[tix[[k]]$i1] <- hm$M[t, k] * a_prev[tix[[k]]$i0]
      for (k2 in seq_len(K)) {
        if (k2 == k) next
        i0 <- tix[[k2]]$i0
        i1 <- tix[[k2]]$i1
        af[i1] <- af[i1] + hm$M[t, k2] * af[i0]
      }
      contrib <- sum(af * w_t * bw$B[t, ])
      wm[t, k] <- contrib * exp(ls_prev + mx + bw$lsB[t] - logL)
    }
  }
  # never flipping is the last interval: the final states without the
  # coordinate, under the survival product they already carry
  for (k in seq_len(K)) {
    i0 <- tix[[k]]$i0
    wm[J, k] <- sum(fw$A[ng, i0] * bw$B[ng, i0]) *
      exp(fw$lsA[ng] + bw$lsB[ng] - logL)
  }
  wm[wm < 0] <- 0
  sw <- colSums(wm)
  for (k in seq_len(K)) if (sw[k] > 0) wm[, k] <- wm[, k] / sw[k]
  wm
}

# the node posterior of the continuous kinds, with the shift matrix the
# levels method reads: both padded to the widest group
.marg_seg_posterior <- function(term, eta, y, logdens, psi) {
  bp <- .marg_built(term)
  v <- .marg_check_psi(term, psi)
  n <- bp$n
  states <- list()
  Cmax <- 0L
  for (g in seq_along(bp$groups)) {
    rs <- bp$groups[[g]]
    nd <- .marg_nodes(bp$x[rs], v[["m1"]], v[["tau1"]])
    sh <- .marg_seg_shift(term, bp$x[rs], nd, v)
    C <- length(nd$p)
    ei <- rep(eta[rs], C) + as.numeric(sh$shift)
    LD <- matrix(as.numeric(logdens(ei, rep(rs, C))), length(rs), C)
    la <- nd$lw + colSums(LD)
    w <- exp(la - .marg_lse(la))
    states[[g]] <- list(rs = rs, w = w, sh = sh, nd = nd)
    Cmax <- max(Cmax, C)
  }
  gamma <- matrix(0, n, Cmax)
  shift <- matrix(0, n, Cmax)
  for (st in states) {
    C <- length(st$w)
    gamma[st$rs, seq_len(C)] <- matrix(st$w, length(st$rs), C, byrow = TRUE)
    shift[st$rs, seq_len(C)] <- st$sh$shift
  }
  list(gamma = gamma, shift = shift, states = states)
}

#' @title Levels of a Marginal Break-Point Term
#' @name term_levels.MarginalBreakTerm
#' @description
#' For the step kind, the constant shift of each side pattern, the sums of
#' the changes of level over the active break-points. For the continuous
#' kinds the shift varies by observation, each node contributing its own
#' hinge value, and
#' a matrix is returned, aligned with [term_posterior()]'s
#' columns; it takes the callbacks because the node set is theirs to
#' rebuild.
#' @param term A built [MarginalBreakTerm()].
#' @param psi The term's parameters.
#' @param eta,y,logdens For the continuous kinds, the quantities the node
#'   set is built from; ignored by the step kind.
#' @param ... Unused.
#' @return A numeric vector, or a matrix with one row per observation.
#' @keywords internal
S7::method(term_levels, MarginalBreakTerm) <- function(term, psi, eta = NULL,
                                                       y = NULL,
                                                       logdens = NULL, ...) {
  v <- .marg_check_psi(term, psi)
  if (term@kind == "jump") {
    return(.marg_bits(term@npsi,
                      v[paste0("delta", seq_len(term@npsi))])$shifts)
  }
  if (is.null(eta) || is.null(logdens)) {
    stop(paste("the continuous kinds shift each observation by its node's",
               "hinge value: pass eta, y and logdens so the node set can be",
               "rebuilt."), call. = FALSE)
  }
  .marg_seg_posterior(term, eta, y, logdens, psi)$shift
}

#' @title Posterior Break-Points of a Marginal Term
#' @name term_latent.MarginalBreakTerm
#' @description
#' The posterior mean and standard deviation of each group's break-points.
#' For the step kind under the gaussian prior these are mixtures of
#' truncated-normal moments over the interval posterior, the edge
#' intervals read through [numericals7::mills_ratio()]; under an
#' explicit prior the truncated moments come from
#' [distributions7::truncated()] and
#' [distributions7::expectation()], one interval at a time. For
#' the continuous kinds they are the moments of the node posterior, the
#' closed upper tail entering through its truncated-normal moments.
#' @param term A built [MarginalBreakTerm()].
#' @param eta The static predictor.
#' @param y The response.
#' @param logdens The log-density.
#' @param psi The term's parameters.
#' @param ... Unused.
#' @return A data frame with `group`, `psi`, `mean` and
#'   `sd`.
#' @keywords internal
S7::method(term_latent, MarginalBreakTerm) <- function(term, eta, y, logdens,
                                                       psi, ...) {
  bp <- .marg_built(term)
  v <- .marg_check_psi(term, psi)
  K <- term@npsi
  n <- bp$n
  idx <- seq_len(n)
  rows <- list()
  if (term@kind == "jump") {
    P <- 2^K
    pb <- .marg_bits(K, v[paste0("delta", seq_len(K))])
    LD <- matrix(0, n, P)
    for (p in seq_len(P)) {
      LD[, p] <- as.numeric(logdens(eta + pb$shifts[p], idx))
    }
    tix <- .marg_hmm_idx(K)
    for (g in seq_along(bp$groups)) {
      rs <- bp$groups[[g]]
      hm <- .marg_hmm_masses(term, bp$x[rs], v)
      wm_all <- .marg_hmm_flip(term, rs, hm, LD, tix, K)
      for (k in seq_len(K)) {
        mo <- .marg_interval_moments(term, bp$x[rs], v, k, wm_all[, k],
                                     list(ivs = hm$ivs))
        rows[[length(rows) + 1L]] <- data.frame(
          group = bp$labels[g], psi = k, mean = mo$mean, sd = mo$sd,
          stringsAsFactors = FALSE)
      }
    }
  } else {
    ps <- .marg_seg_posterior(term, eta, y, logdens, psi)
    for (g in seq_along(ps$states)) {
      st <- ps$states[[g]]
      w <- st$w
      p <- st$nd$p
      C <- length(p)
      fin <- is.finite(p)
      m1 <- sum(w[fin] * p[fin])
      m2 <- sum(w[fin] * p[fin]^2)
      # the closed upper tail: the truncated normal above the last
      # observation
      wt <- w[C]
      if (wt > 0) {
        xn <- max(bp$x[st$rs])
        zn <- (xn - v[["m1"]]) / v[["tau1"]]
        lam <- numericals7::mills_ratio(-zn)$r
        mu_t <- v[["m1"]] + v[["tau1"]] * lam
        var_t <- v[["tau1"]]^2 * (1 + zn * lam - lam^2)
        m1 <- m1 + wt * mu_t
        m2 <- m2 + wt * (var_t + mu_t^2)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        group = bp$labels[g], psi = 1L, mean = m1,
        sd = sqrt(max(m2 - m1^2, 0)), stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

# the truncated moments of one break-point's prior over its intervals,
# mixed by the posterior marginal
.marg_interval_moments <- function(term, xs, v, k, wm, pr) {
  J <- length(wm)
  if (is.null(term@prior)) {
    m <- v[[paste0("m", k)]]
    tau <- v[[paste0("tau", k)]]
    iv <- pr$ivs[[k]]
    r1 <- ifelse(iv$ok, (iv$phi_l - iv$phi_u) / iv$mass, 0)
    r2 <- ifelse(iv$ok, (iv$zp_l - iv$zp_u) / iv$mass, 0)
    mr_u <- numericals7::mills_ratio(iv$zu[1L])$r
    mr_l <- numericals7::mills_ratio(-iv$zl[J])$r
    r1[1L] <- -mr_u
    r2[1L] <- -iv$zu[1L] * mr_u
    r1[J] <- mr_l
    r2[J] <- iv$zl[J] * mr_l
    mean_j <- m + tau * r1
    var_j <- tau^2 * (1 + r2 - r1^2)
    keep <- wm > 0
    mu <- sum(wm[keep] * mean_j[keep])
    m2 <- sum(wm[keep] * (var_j[keep] + mean_j[keep]^2))
    return(list(mean = mu, sd = sqrt(max(m2 - mu^2, 0))))
  }
  # an explicit prior: the moments of the prior truncated to each interval
  # that carries appreciable weight, through the batched engines. A moment
  # the engine cannot deliver is reported NA rather than approximated: a
  # heavy-tailed prior's edge intervals keep no mean below one degree of
  # freedom and no variance below two, and the quadrature's refusal is the
  # honest reading of that.
  iv <- pr$ivs[[k]]
  th <- stats::setNames(as.list(v[term@prior@params]), term@prior@params)
  mu <- 0
  m2 <- 0
  for (j in seq_len(J)) {
    if (wm[j] < 1e-9) next
    # an infinite endpoint is omitted rather than passed: truncating at the
    # support's own bound removes no mass and truncated() rejects it
    tr_args <- list(term@prior)
    if (is.finite(iv$l[j])) tr_args$lower <- iv$l[j] - v[["m1"]]
    if (is.finite(iv$u[j])) tr_args$upper <- iv$u[j] - v[["m1"]]
    tr <- do.call(distributions7::truncated, tr_args)
    e1 <- tryCatch(as.numeric(distributions7::expectation(
      tr, function(y, theta) y, th)), error = function(e) NA_real_)
    e2 <- tryCatch(as.numeric(distributions7::expectation(
      tr, function(y, theta) y^2, th)), error = function(e) NA_real_)
    mu <- mu + wm[j] * (v[["m1"]] + e1)
    m2 <- m2 + wm[j] * (e2 + 2 * v[["m1"]] * e1 + v[["m1"]]^2)
  }
  sc <- sum(wm[wm >= 1e-9])
  list(mean = mu / sc,
       sd = if (is.finite(m2) && is.finite(mu)) {
         sqrt(max(m2 / sc - (mu / sc)^2, 0))
       } else NA_real_)
}

#' @title Where a Marginal Break-Point Term Starts
#' @name term_start.MarginalBreakTerm
#' @description
#' Given a target (the response on the scale of the predictor, which the
#' fitting layer supplies), the exact profile: each group's best split by
#' least squares on its own observations under the term's own design, the
#' population positions and the prior scales read off those, and the
#' changes their pooled means. Without one, the covariate's quantiles and
#' a scale from its spread, with the changes at zero.
#' @param term A built [MarginalBreakTerm()].
#' @param ... Unused.
#' @param target The response on the predictor's scale, or `NULL`.
#' @return A named numeric vector on the unconstrained scale.
#' @keywords internal
S7::method(term_start, MarginalBreakTerm) <- function(term, ...,
                                                      target = NULL) {
  nm <- term_params(term)
  links <- term_links(term)
  out <- stats::setNames(numeric(length(nm)), nm)
  bp <- term@blueprint
  if (!length(bp)) return(out)
  xv <- bp$x
  K <- term@npsi
  qs <- stats::quantile(xv, c(0.05, 0.95), names = FALSE)
  m0 <- if (!is.null(term@spec$psi)) sort(term@spec$psi) else {
    as.numeric(stats::quantile(xv, seq_len(K) / (K + 1), names = FALSE))
  }
  tau0 <- rep((qs[2L] - qs[1L]) / 6, K)
  ch <- list(delta = rep(0, K), gamma = rep(0, K), beta = 0)
  if (!is.null(target) && length(target) == bp$n) {
    prof <- .marg_profile_start(term, bp, as.numeric(target))
    if (!is.null(prof)) {
      if (is.null(term@spec$psi)) m0 <- prof$m
      tau0 <- prof$tau
      ch <- prof$ch
    }
  }
  tau0 <- pmax(tau0, .marg_gap_floor(bp))
  for (k in seq_len(K)) {
    if (is.null(term@prior)) {
      out[[paste0("m", k)]] <- m0[k]
      out[[paste0("tau", k)]] <- log(tau0[k])
    }
  }
  if (!is.null(term@prior)) {
    out[["m1"]] <- m0[1L]
    # a prior parameter read as a scale starts at the profile's spread,
    # carried through its own chart; the others keep the chart's zero
    ip <- tryCatch(term@prior@params_interpretation, error = function(e) NULL)
    for (p in term@prior@params) {
      if (!is.null(ip) && !is.null(ip[[p]]) &&
          ip[[p]] %in% c("scale", "standard deviation")) {
        z <- tryCatch(linkfunctions7::linkfun(links[[p]], tau0[1L]),
                      error = function(e) NA_real_)
        if (is.finite(z)) out[[p]] <- z
      }
    }
  }
  if (term@linear) out[["beta"]] <- ch$beta
  for (k in seq_len(K)) {
    if (term@kind %in% c("jump", "jseg")) {
      out[[paste0("delta", k)]] <- ch$delta[k]
    }
    if (term@kind %in% c("seg", "jseg")) {
      out[[paste0("gamma", k)]] <- ch$gamma[k]
    }
  }
  out
}

# Half the median gap between a group's consecutive distinct covariate
# values, pooled over the groups: the finest scale at which the data can
# tell two positions apart, and the floor under a starting prior scale.
.marg_gap_floor <- function(bp) {
  gaps <- unlist(lapply(bp$groups, function(rs) {
    u <- sort(unique(bp$x[rs]))
    if (length(u) > 1L) diff(u) else numeric(0)
  }))
  if (!length(gaps)) gaps <- diff(range(bp$x))
  max(stats::median(gaps) / 2, sqrt(.Machine$double.eps))
}

# The exact least-squares profile of the target under the term's own
# design, in two stages. The POOLED stage estimates the population
# positions and the changes on all the groups at once, per-group
# intercepts absorbing whatever levels the target carries, the positions
# taken greedily off a quantile grid -- one group's own profile is not
# used for the coefficients, because a per-group fit with the full design
# on a noisy target overfits its own noise (measured on a Poisson panel:
# per-group linear effects with median -4.7 against a truth of 0.15, and
# the search then started in the wrong basin). The PER-GROUP stage moves
# one position at a time with everything else held at the pooled
# estimates, which is one free value per group and cannot overfit; the
# spread of those positions is the prior scale's start.
.marg_profile_start <- function(term, bp, target) {
  K <- term@npsi
  n <- bp$n
  xv <- bp$x
  keep <- is.finite(target)
  if (sum(keep) < K + 3L) return(NULL)
  G <- length(bp$groups)
  gvec <- integer(n)
  for (g in seq_len(G)) gvec[bp$groups[[g]]] <- g
  Gm <- matrix(0, n, G)
  Gm[cbind(seq_len(n), gvec)] <- 1

  cols_at <- function(cs) {
    X <- if (term@linear) cbind(Gm, xv) else Gm
    for (c0 in cs) {
      if (term@kind %in% c("seg", "jseg")) X <- cbind(X, pmax(xv - c0, 0))
      if (term@kind %in% c("jump", "jseg")) {
        X <- cbind(X, as.numeric(xv >= c0))
      }
    }
    X
  }
  cand <- unique(as.numeric(stats::quantile(
    xv[keep], seq(0.06, 0.94, length.out = 23L), names = FALSE)))
  pooled_rss <- function(cs) {
    X <- cols_at(cs)
    fit <- tryCatch(stats::lm.fit(X[keep, , drop = FALSE], target[keep]),
                    error = function(e) NULL)
    if (is.null(fit)) return(list(rss = Inf, cf = NULL))
    list(rss = sum(fit$residuals^2), cf = fit$coefficients)
  }

  # the columns of the final scoring fit, with PER-GROUP positions: shared
  # coefficients over each group's own hinge and step
  cols_per_group <- function(psis) {
    X <- if (term@linear) cbind(Gm, xv) else Gm
    pg <- psis[gvec, , drop = FALSE]
    for (k in seq_len(K)) {
      if (term@kind %in% c("seg", "jseg")) {
        X <- cbind(X, pmax(xv - pg[, k], 0))
      }
      if (term@kind %in% c("jump", "jseg")) {
        X <- cbind(X, as.numeric(xv >= pg[, k]))
      }
    }
    X
  }
  read_ch <- function(cf) {
    cf[!is.finite(cf)] <- 0
    pos <- G + as.integer(term@linear)
    beta0 <- if (term@linear) cf[G + 1L] else 0
    gamma0 <- numeric(K)
    delta0 <- numeric(K)
    for (k in seq_len(K)) {
      if (term@kind == "seg") gamma0[k] <- cf[pos + k]
      if (term@kind == "jump") delta0[k] <- cf[pos + k]
      if (term@kind == "jseg") {
        gamma0[k] <- cf[pos + 2L * k - 1L]
        delta0[k] <- cf[pos + 2L * k]
      }
    }
    list(beta = unname(beta0), gamma = gamma0, delta = delta0)
  }
  refine <- function(found, ch) {
    contrib <- function(cs) {
      v <- if (term@linear) ch$beta * xv else numeric(n)
      for (k in seq_along(cs)) {
        if (term@kind %in% c("seg", "jseg")) {
          v <- v + ch$gamma[k] * pmax(xv - cs[k], 0)
        }
        if (term@kind %in% c("jump", "jseg")) {
          v <- v + ch$delta[k] * as.numeric(xv >= cs[k])
        }
      }
      v
    }
    psis <- matrix(NA_real_, G, K)
    for (g in seq_len(G)) {
      rs <- bp$groups[[g]]
      okr <- rs[is.finite(target[rs])]
      if (length(okr) < 3L) next
      xs <- xv[okr]
      dist <- which(xs[-length(xs)] < xs[-1L])
      if (!length(dist)) next
      mids <- (xs[dist] + xs[dist + 1L]) / 2
      for (k in seq_len(K)) {
        best <- Inf
        bc <- found[k]
        for (c0 in mids) {
          cs <- found
          cs[k] <- c0
          r <- target[okr] - contrib(cs)[okr]
          rss <- sum((r - mean(r))^2)
          if (rss < best) {
            best <- rss
            bc <- c0
          }
        }
        psis[g, k] <- bc
      }
      psis[g, ] <- sort(psis[g, ])
    }
    psis
  }

  # The pooled single-position profile is nearly flat when the positions
  # vary by group -- their spread smears the break -- and its global
  # minimum is not reliably in the right basin (measured on a Poisson
  # panel: 103.6 at c = 8.3 against 104.8 at c = 4.4, the truth's basin
  # the runner-up). Every local minimum is therefore carried through the
  # per-group refinement, and the winner is the candidate whose FINAL fit,
  # shared coefficients over per-group positions, leaves the least behind:
  # at per-group positions the smearing is gone and the true basin wins.
  scan <- vapply(cand, function(c0) pooled_rss(c0)$rss, numeric(1))
  nc <- length(scan)
  locmin <- which(scan <= c(Inf, scan[-nc]) & scan <= c(scan[-1L], Inf) &
                    is.finite(scan))
  if (!length(locmin)) return(NULL)
  locmin <- locmin[order(scan[locmin])]
  locmin <- locmin[seq_len(min(4L, length(locmin)))]

  best <- NULL
  for (i0 in locmin) {
    found <- cand[i0]
    ok_seed <- TRUE
    if (K > 1L) {
      for (k in seq.int(2L, K)) {
        rssk <- vapply(setdiff(cand, found), function(c0)
          pooled_rss(c(found, c0))$rss, numeric(1))
        if (!any(is.finite(rssk))) {
          ok_seed <- FALSE
          break
        }
        found <- c(found, setdiff(cand, found)[which.min(rssk)])
      }
    }
    if (!ok_seed) next
    found <- sort(found)
    pf <- pooled_rss(found)
    if (is.null(pf$cf)) next
    ch <- read_ch(pf$cf)
    psis <- refine(found, ch)
    okg <- stats::complete.cases(psis)
    if (!any(okg)) next
    for (g in which(!okg)) psis[g, ] <- found
    Xf <- cols_per_group(psis)
    ff <- tryCatch(stats::lm.fit(Xf[keep, , drop = FALSE], target[keep]),
                   error = function(e) NULL)
    if (is.null(ff)) next
    rssf <- sum(ff$residuals^2)
    if (is.null(best) || rssf < best$rss) {
      best <- list(rss = rssf, psis = psis, ch = read_ch(ff$coefficients),
                   okg = okg)
    }
  }
  if (is.null(best)) return(NULL)
  psis <- best$psis[best$okg, , drop = FALSE]
  tau <- apply(psis, 2L, stats::mad)
  bad <- !is.finite(tau) | tau <= 0
  if (any(bad)) tau[bad] <- apply(psis, 2L, stats::sd, na.rm = TRUE)
  list(m = apply(psis, 2L, stats::median),
       tau = pmax(ifelse(is.finite(tau), tau, 0), 0),
       ch = best$ch)
}

# The observed Hessian, analytic throughout. The step kind propagates first
# and second derivatives through the side chain's forward recursion: the
# transition weights carry the prior's interval masses, whose derivatives
# are closed for the gaussian and come from the cdf surface for an explicit
# prior, and the emissions carry the family's derivatives at the pattern
# shifts. The continuous kinds differentiate the node sum twice, the node
# motion of the panels below the data being AFFINE in the prior's
# parameters, so the chain rule closes with no curvature term from the
# nodes themselves. The one-break-point gaussian step keeps the interval-sum
# route, which shares none of this arithmetic and is the twin the tests
# hold the propagation to.

#' @title Observed Hessian of a Marginal Break-Point Term
#' @name term_hessian.MarginalBreakTerm
#' @description
#' The exact observed Hessian of the marginal log-likelihood over a
#' caller's unknowns, the coefficients of every equation together with the
#' term's own parameters on the unconstrained scale.
#' @details
#' Analytic throughout. The step kind propagates first and second
#' derivatives through the side chain's forward recursion, the prior's
#' interval-mass derivatives closed for the gaussian and read off the cdf
#' surface for an explicit prior (whose own degrees-of-freedom column
#' carries that surface's documented single stencil, the one non-closed
#' piece anywhere). The continuous kinds differentiate the node sum twice;
#' the moving panels below the data are affine in the prior's parameters,
#' so their motion enters the chain rule with no curvature of its own. The
#' one-break-point gaussian step keeps the interval-sum route, independent
#' arithmetic the tests hold the propagation to.
#'
#' The marginal likelihood of a group does not factorize over its
#' observations, so an observation weight has a reading only when it is
#' constant within each group; anything else is rejected.
#' @param term A built [MarginalBreakTerm()].
#' @param eta The static predictor of the level equation.
#' @param y The response.
#' @param logdens,grad,hess The log-density and its first two derivatives
#'   in the predictors.
#' @param psi The term's parameters.
#' @param seed The derivative of each predictor in the caller's unknowns.
#' @param cols The columns the term's own parameters occupy.
#' @param level The distribution parameter the term shifts.
#' @param weights Observation weights, constant within each group.
#' @param ... Unused.
#' @return A list with `loglik`, `gradient` and `hessian`.
#' @keywords internal
S7::method(term_hessian, MarginalBreakTerm) <- function(term, eta, y, logdens,
                                                        grad, hess, psi, seed,
                                                        cols, level,
                                                        weights = NULL, ...) {
  bp <- .marg_built(term)
  v <- .marg_check_psi(term, psi)
  nm <- term_params(term)
  n <- bp$n
  if (!is.list(seed) || !length(seed)) {
    stop("'seed' must be a list with one matrix per distribution parameter.",
         call. = FALSE)
  }
  seed <- lapply(seed, as.matrix)
  mm <- ncol(seed[[1L]])
  npar <- length(seed)
  if (any(vapply(seed, nrow, 1L) != n) ||
      any(vapply(seed, ncol, 1L) != mm)) {
    stop(sprintf("every element of 'seed' must be %d by %d.", n, mm),
         call. = FALSE)
  }
  cols <- as.integer(cols)
  level <- as.integer(level)
  if (length(cols) != length(nm) || any(cols < 1L) || any(cols > mm)) {
    stop(sprintf("'cols' must give the %d columns the term's parameters",
                 length(nm)), call. = FALSE)
  }
  if (length(level) != 1L || level < 1L || level > npar) {
    stop("'level' must index one of the distribution parameters.",
         call. = FALSE)
  }
  w <- if (is.null(weights)) rep(1, n) else rep_len(as.numeric(weights), n)
  for (rs in bp$groups) {
    if (length(unique(w[rs])) > 1L) {
      stop(paste("a marginal likelihood does not factorize over a group's",
                 "observations, so an observation weight has a reading only",
                 "when it is constant within each group."), call. = FALSE)
    }
  }

  if (term@kind == "jump" && term@npsi == 1L && is.null(term@prior)) {
    return(.marg_jump1_hessian(term, eta, y, logdens, grad, hess, psi, seed,
                               cols, level, w))
  }
  if (term@kind == "jump") {
    .marg_jump_hessian(term, eta, y, logdens, grad, hess, psi, seed, cols,
                       level, w)
  } else {
    .marg_seg_hessian(term, eta, y, logdens, grad, hess, psi, seed, cols,
                      level, w)
  }
}

# The propagated Hessian of the step kind: the side chain's forward
# recursion differentiated twice in the caller's unknowns. Transitions are
# linear with the prior's interval masses as weights, emissions multiply by
# the pattern density, and the survival products of the final states carry
# the tails; each is differentiated by its own product rule, and the mass
# and survival derivative columns are chained onto the unconstrained scale
# before the propagation so the result needs no chart correction.
.marg_jump_hessian <- function(term, eta, y, logdens, grad, hess, psi, seed,
                               cols, level, w) {
  bp <- term@blueprint
  v <- .marg_check_psi(term, psi)
  nm <- term_params(term)
  links <- term_links(term)
  K <- term@npsi
  P <- 2^K
  n <- bp$n
  npar <- length(seed)
  mm <- ncol(seed[[1L]])
  idx <- seq_len(n)
  pb <- .marg_bits(K, v[paste0("delta", seq_len(K))])
  LD <- matrix(0, n, P)
  G <- vector("list", P)
  H2 <- vector("list", P)
  for (p in seq_len(P)) {
    LD[, p] <- as.numeric(logdens(eta + pb$shifts[p], idx))
    G[[p]] <- as.matrix(grad(eta + pb$shifts[p], idx))
    H2[[p]] <- array(as.numeric(hess(eta + pb$shifts[p], idx)),
                     c(n, npar, npar))
  }
  dcol <- vapply(seq_len(K), function(k)
    cols[match(paste0("delta", k), nm)], integer(1))
  # per-pattern per-observation scores in the unknowns, the pattern's
  # active break-points adding the level score to their delta columns
  gm <- vector("list", P)
  for (p in seq_len(P)) {
    g0 <- matrix(0, n, mm)
    for (q in seq_len(npar)) g0 <- g0 + G[[p]][, q] * seed[[q]]
    for (k in seq_len(K)) {
      if (pb$bits[p, k]) g0[, dcol[k]] <- g0[, dcol[k]] + G[[p]][, level]
    }
    gm[[p]] <- g0
  }
  tix <- .marg_hmm_idx(K)
  zch <- vapply(nm, function(j)
    linkfunctions7::linkfun(links[[j]], v[[j]]), numeric(1))
  J1 <- vapply(nm, function(j)
    linkfunctions7::dlinkinv(links[[j]], zch[[j]]), numeric(1))
  J2 <- vapply(nm, function(j)
    linkfunctions7::d2linkinv(links[[j]], zch[[j]]), numeric(1))

  loglik <- numeric(n)
  gradient <- numeric(mm)
  hessian <- matrix(0, mm, mm)
  for (rs in bp$groups) {
    ng <- length(rs)
    wi <- w[rs][1L]
    hm <- .marg_hmm_masses(term, bp$x[rs], v, d2 = TRUE)
    pcol <- lapply(hm$pcols, function(pc) cols[match(pc, nm)])
    # the masses' and survivals' derivative columns onto the unconstrained
    # scale, second order carrying the chart's own curvature
    for (k in seq_len(K)) {
      ii <- match(hm$pcols[[k]], nm)
      Jk <- J1[ii]
      Dk <- J2[ii]
      dpar <- hm$dM[[k]]
      hm$dM[[k]] <- sweep(dpar, 2L, Jk, `*`)
      A <- hm$d2M[[k]]
      for (i in seq_along(Jk)) {
        for (j in seq_along(Jk)) {
          A[, i, j] <- A[, i, j] * Jk[i] * Jk[j] +
            (i == j) * dpar[, i] * Dk[i]
        }
      }
      hm$d2M[[k]] <- A
      dspar <- hm$dSV[[k]]
      hm$dSV[[k]] <- sweep(dspar, 2L, Jk, `*`)
      B <- hm$d2SV[[k]]
      for (i in seq_along(Jk)) {
        for (j in seq_along(Jk)) {
          B[, i, j] <- B[, i, j] * Jk[i] * Jk[j] +
            (i == j) * dspar[, i] * Dk[i]
        }
      }
      hm$d2SV[[k]] <- B
    }

    alpha <- numeric(P)
    alpha[1L] <- 1
    Dal <- matrix(0, P, mm)
    D2al <- array(0, c(P, mm, mm))
    ls <- 0
    lnum_prev <- 0
    for (t in seq_len(ng)) {
      row <- rs[t]
      for (k in seq_len(K)) {
        i0 <- tix[[k]]$i0
        i1 <- tix[[k]]$i1
        q <- hm$M[t, k]
        pc <- pcol[[k]]
        npc <- length(pc)
        dq <- hm$dM[[k]][t, ]
        d2q <- matrix(hm$d2M[[k]][t, , ], npc, npc)
        D2al[i1, , ] <- D2al[i1, , ] + q * D2al[i0, , ]
        for (i in seq_len(npc)) {
          D2al[i1, pc[i], ] <- D2al[i1, pc[i], ] + dq[i] * Dal[i0, ]
          D2al[i1, , pc[i]] <- D2al[i1, , pc[i]] + dq[i] * Dal[i0, ]
          for (j in seq_len(npc)) {
            D2al[i1, pc[i], pc[j]] <- D2al[i1, pc[i], pc[j]] +
              alpha[i0] * d2q[i, j]
          }
        }
        Dal[i1, ] <- Dal[i1, ] + q * Dal[i0, ]
        Dal[i1, pc] <- Dal[i1, pc] + alpha[i0] %o% dq
        alpha[i1] <- alpha[i1] + q * alpha[i0]
      }
      lf <- LD[row, ]
      mx <- max(lf)
      ls <- ls + mx
      As0 <- matrix(0, npar, mm)
      for (q in seq_len(npar)) As0[q, ] <- seed[[q]][row, ]
      for (s in seq_len(P)) {
        e <- exp(lf[s] - mx)
        glf <- gm[[s]][row, ]
        As <- As0
        for (k in seq_len(K)) {
          if (pb$bits[s, k]) As[level, dcol[k]] <- As[level, dcol[k]] + 1
        }
        Hlf <- crossprod(As, matrix(H2[[s]][row, , ], npar, npar) %*% As)
        d <- Dal[s, ]
        D2al[s, , ] <- e * (D2al[s, , ] + outer(d, glf) + outer(glf, d) +
                              alpha[s] * (outer(glf, glf) + Hlf))
        Dal[s, ] <- e * (d + alpha[s] * glf)
        alpha[s] <- e * alpha[s]
      }
      sp <- .marg_hmm_sp(tix, hm$SV[t + 1L, ], P)
      lnum <- log(sum(alpha * sp)) + ls
      loglik[row] <- lnum - lnum_prev
      lnum_prev <- lnum
      cs <- max(alpha)
      if (cs > 0 && (cs > 1e100 || cs < 1e-100)) {
        alpha <- alpha / cs
        Dal <- Dal / cs
        D2al <- D2al / cs
        ls <- ls + log(cs)
      }
    }
    # the total and its two derivatives, the tails entering through the
    # survival products of the final states; the per-coordinate products
    # are built by exclusion, so a zero survival divides nothing
    Ft <- 0
    DF <- numeric(mm)
    D2F <- matrix(0, mm, mm)
    SVv <- hm$SV[ng + 1L, ]
    for (s in seq_len(P)) {
      inact <- which(!pb$bits[s, ])
      sp_s <- prod(SVv[inact])
      dsp <- numeric(mm)
      d2sp <- matrix(0, mm, mm)
      for (k in inact) {
        pk <- prod(SVv[setdiff(inact, k)])
        dsp[pcol[[k]]] <- dsp[pcol[[k]]] + pk * hm$dSV[[k]][ng + 1L, ]
        npc <- length(pcol[[k]])
        d2sp[pcol[[k]], pcol[[k]]] <- d2sp[pcol[[k]], pcol[[k]]] +
          pk * matrix(hm$d2SV[[k]][ng + 1L, , ], npc, npc)
        for (k2 in inact) {
          if (k2 == k) next
          pkk <- prod(SVv[setdiff(inact, c(k, k2))])
          d2sp[pcol[[k]], pcol[[k2]]] <- d2sp[pcol[[k]], pcol[[k2]]] +
            pkk * hm$dSV[[k]][ng + 1L, ] %o% hm$dSV[[k2]][ng + 1L, ]
        }
      }
      Ft <- Ft + alpha[s] * sp_s
      DF <- DF + Dal[s, ] * sp_s + alpha[s] * dsp
      D2F <- D2F + sp_s * matrix(D2al[s, , ], mm, mm) +
        outer(Dal[s, ], dsp) + outer(dsp, Dal[s, ]) + alpha[s] * d2sp
    }
    g0 <- DF / Ft
    gradient <- gradient + wi * g0
    hessian <- hessian + wi * (D2F / Ft - outer(g0, g0))
  }
  hessian <- (hessian + t(hessian)) / 2
  list(loglik = loglik, gradient = gradient, hessian = hessian)
}

# The Hessian of the continuous kinds: the node sum differentiated twice.
# The conditional's second derivatives collapse into per-node
# per-observation Hessians whose level rows carry the node's own shift
# partials, node motion included; the first-derivative products run over
# the node posterior; and the weights' own curvature is closed, the node
# motion being affine in the prior's parameters so the chain rule ends at
# first order in it. What remains is the chart's curvature on the log
# scale of the prior's spread, added where the identity ends.
.marg_seg_hessian <- function(term, eta, y, logdens, grad, hess, psi, seed,
                              cols, level, w) {
  bp <- term@blueprint
  v <- .marg_check_psi(term, psi)
  nm <- term_params(term)
  links <- term_links(term)
  n <- bp$n
  npar <- length(seed)
  mm <- ncol(seed[[1L]])
  chain <- vapply(nm, function(j)
    linkfunctions7::dlinkinv(links[[j]],
                             linkfunctions7::linkfun(links[[j]], v[[j]])),
    numeric(1))
  ps <- .marg_seg_posterior(term, eta, y, logdens, as.list(v))
  loglik <- numeric(n)
  gradient <- numeric(mm)
  hessian <- matrix(0, mm, mm)

  for (st in ps$states) {
    rs <- st$rs
    ng <- length(rs)
    wi <- w[rs][1L]
    C <- length(st$w)
    ei <- rep(eta[rs], C) + as.numeric(st$sh$shift)
    ii <- rep(rs, C)
    LD <- matrix(as.numeric(logdens(ei, ii)), ng, C)
    Gf <- as.matrix(grad(ei, ii))
    Hf <- array(as.numeric(hess(ei, ii)), c(ng * C, npar, npar))
    GL <- matrix(Gf[, level], ng, C)

    # the contributions, from the same states
    A <- st$nd$lw
    tot <- .marg_lse(A)
    for (t in seq_len(ng)) {
      A <- A + LD[t, ]
      tot2 <- .marg_lse(A)
      loglik[rs[t]] <- tot2 - tot
      tot <- tot2
    }

    # the per-node partials of the shift in the term's own coefficients,
    # and the level column each node's shift augments
    ownm <- list(beta = if (term@linear) matrix(bp$x[rs], ng, C),
                 gamma1 = st$sh$hinge,
                 delta1 = if (term@kind == "jseg") st$sh$step)
    ownm <- ownm[!vapply(ownm, is.null, logical(1))]

    # the per-node gradient in the unknowns
    DA <- matrix(0, C, mm)
    for (q in seq_len(npar)) {
      Gq <- matrix(Gf[, q], ng, C)
      DA <- DA + crossprod(Gq, seed[[q]][rs, , drop = FALSE])
    }
    for (p in names(ownm)) {
      dc <- cols[match(p, nm)]
      DA[, dc] <- DA[, dc] + colSums(GL * ownm[[p]]) * chain[match(p, nm)]
    }
    cpsi <- colSums(GL * st$sh$dshift_dpsi)
    im <- cols[match("m1", nm)]
    it <- cols[match("tau1", nm)]
    DA[, im] <- DA[, im] +
      (st$nd$glw_m + cpsi * st$nd$dpsi_m) * chain[match("m1", nm)]
    DA[, it] <- DA[, it] +
      (st$nd$glw_t + cpsi * st$nd$dpsi_t) * chain[match("tau1", nm)]
    gbar <- as.numeric(crossprod(DA, st$w))
    gradient <- gradient + wi * gbar
    hessian <- hessian + wi * (crossprod(DA, st$w * DA) - outer(gbar, gbar))

    # the conditional's second derivatives at each node, the level row of
    # the per-observation Jacobian augmented by the node's own partials
    for (c0 in seq_len(C)) {
      wc <- wi * st$w[c0]
      if (wc == 0) next
      rows0 <- (c0 - 1L) * ng + seq_len(ng)
      seedc <- lapply(seed, function(s) s[rs, , drop = FALSE])
      aug <- matrix(0, ng, mm)
      for (p in names(ownm)) {
        aug[, cols[match(p, nm)]] <- ownm[[p]][, c0] * chain[match(p, nm)]
      }
      aug[, im] <- aug[, im] +
        st$sh$dshift_dpsi[, c0] * st$nd$dpsi_m[c0] * chain[match("m1", nm)]
      aug[, it] <- aug[, it] +
        st$sh$dshift_dpsi[, c0] * st$nd$dpsi_t[c0] * chain[match("tau1", nm)]
      seedc[[level]] <- seedc[[level]] + aug
      for (aq in seq_len(npar)) {
        for (bq in seq_len(npar)) {
          hessian <- hessian +
            crossprod(seedc[[aq]], (wc * Hf[rows0, aq, bq]) * seedc[[bq]])
        }
      }
    }

    # the weights' own curvature, closed because the node motion is affine
    Jm <- chain[match("m1", nm)]
    Jt <- chain[match("tau1", nm)]
    a_mm <- sum(st$w * st$nd$alw_mm)
    a_mt <- sum(st$w * st$nd$alw_mt)
    a_tt <- sum(st$w * st$nd$alw_tt)
    hessian[im, im] <- hessian[im, im] + wi * a_mm * Jm * Jm
    hessian[im, it] <- hessian[im, it] + wi * a_mt * Jm * Jt
    hessian[it, im] <- hessian[it, im] + wi * a_mt * Jm * Jt
    hessian[it, it] <- hessian[it, it] + wi * a_tt * Jt * Jt
    # the shift's mixed second derivative, -1(x > psi) between the change
    # of slope and a moving node
    ig <- cols[match("gamma1", nm)]
    crossg <- colSums(GL * -(st$sh$hinge > 0))
    vm <- sum(st$w * crossg * st$nd$dpsi_m) * Jm
    vt <- sum(st$w * crossg * st$nd$dpsi_t) * Jt
    hessian[ig, im] <- hessian[ig, im] + wi * vm
    hessian[im, ig] <- hessian[im, ig] + wi * vm
    hessian[ig, it] <- hessian[ig, it] + wi * vt
    hessian[it, ig] <- hessian[it, ig] + wi * vt
    # the chart's own curvature, on the parameter-scale gradient
    for (pown in nm) {
      d2l <- linkfunctions7::d2linkinv(
        links[[pown]], linkfunctions7::linkfun(links[[pown]], v[[pown]]))
      if (identical(d2l, 0) || d2l == 0) next
      cpn <- cols[match(pown, nm)]
      hessian[cpn, cpn] <- hessian[cpn, cpn] +
        wi * (gbar[cpn] / chain[match(pown, nm)]) * d2l
    }
  }
  hessian <- (hessian + t(hessian)) / 2
  list(loglik = loglik, gradient = gradient, hessian = hessian)
}

# the fully propagated Hessian of the one-break-point gaussian step model:
# the interval sum differentiated twice
.marg_jump1_hessian <- function(term, eta, y, logdens, grad, hess, psi, seed,
                                cols, level, w) {
  bp <- term@blueprint
  v <- .marg_check_psi(term, psi)
  m <- v[["m1"]]
  tau <- v[["tau1"]]
  delta <- v[["delta1"]]
  n <- bp$n
  npar <- length(seed)
  mm <- ncol(seed[[1L]])
  idx <- seq_len(n)
  LF0 <- as.numeric(logdens(eta, idx))
  LF1 <- as.numeric(logdens(eta + delta, idx))
  G0 <- as.matrix(grad(eta, idx))
  G1 <- as.matrix(grad(eta + delta, idx))
  H0 <- array(as.numeric(hess(eta, idx)), c(n, npar, npar))
  H1 <- array(as.numeric(hess(eta + delta, idx)), c(n, npar, npar))
  if (nrow(G0) != n || ncol(G0) != npar) {
    stop(sprintf("'grad' must return an %d by %d matrix.", n, npar),
         call. = FALSE)
  }

  # the per-observation emission scores in the caller's unknowns, one per
  # side of the break-point: the shifted side's delta column carries the
  # level equation's own score
  gm0 <- matrix(0, n, mm)
  gm1 <- matrix(0, n, mm)
  for (q in seq_len(npar)) {
    gm0 <- gm0 + G0[, q] * seed[[q]]
    gm1 <- gm1 + G1[, q] * seed[[q]]
  }
  ic <- cols[match(c("m1", "tau1", "delta1"), term_params(term))]
  gm1[, ic[3L]] <- gm1[, ic[3L]] + G1[, level]

  loglik <- numeric(n)
  gradient <- numeric(mm)
  hessian <- matrix(0, mm, mm)
  P1 <- numeric(n)

  for (rs in bp$groups) {
    ng <- length(rs)
    wi <- w[rs][1L]
    iv <- .marg_intervals(bp$x[rs], m, tau, d2 = TRUE)
    cvec <- c(0, cumsum(LF0[rs] - LF1[rs])) + sum(LF1[rs])
    a <- iv$lm + cvec
    logL <- .marg_lse(a)
    wj <- exp(a - logL)
    P1[rs] <- cumsum(wj)[seq_len(ng)]

    # the interval-level gradients, on the unconstrained scale: the
    # conditional part by the same one-observation-per-interval swap as the
    # values, the mass part chained onto tau's log chart
    M <- gm0[rs, , drop = FALSE] - gm1[rs, , drop = FALSE]
    CS <- apply(M, 2L, cumsum)
    if (is.null(dim(CS))) CS <- matrix(CS, nrow = 1L)
    base <- colSums(gm1[rs, , drop = FALSE])
    DA <- matrix(base, ng + 1L, mm, byrow = TRUE)
    DA[-1L, ] <- DA[-1L, ] + CS
    DA[, ic[1L]] <- DA[, ic[1L]] + iv$dm
    DA[, ic[2L]] <- DA[, ic[2L]] + iv$dt * tau
    gbar <- colSums(wj * DA)
    gradient <- gradient + wi * gbar
    hessian <- hessian + wi * (crossprod(DA, wj * DA) - outer(gbar, gbar))

    # the masses' own curvature, in the (m, zeta_tau) block: the log chart
    # contributes tau^2 times the second derivative plus tau times the first
    smm <- sum(wj * iv$dmm)
    smt <- sum(wj * iv$dmt) * tau
    stt <- sum(wj * iv$dtt) * tau^2 + sum(wj * iv$dt) * tau
    i1 <- ic[1L]
    i2 <- ic[2L]
    hessian[i1, i1] <- hessian[i1, i1] + wi * smm
    hessian[i1, i2] <- hessian[i1, i2] + wi * smt
    hessian[i2, i1] <- hessian[i2, i1] + wi * smt
    hessian[i2, i2] <- hessian[i2, i2] + wi * stt

    # the predictive contributions, from the same states
    A <- iv$lm
    tot <- .marg_lse(A)
    for (t in seq_len(ng)) {
      row <- rs[t]
      A <- A + c(rep(LF1[row], t), rep(LF0[row], ng + 1L - t))
      tot2 <- .marg_lse(A)
      loglik[row] <- tot2 - tot
      tot <- tot2
    }
  }

  # sum_j w_j d2c_j collapses over the intervals: each observation's second
  # derivative enters weighted by the posterior probability of its side
  seed1 <- seed
  seed1[[level]][, ic[3L]] <- seed1[[level]][, ic[3L]] + 1
  for (s in 0:1) {
    Hs <- if (s == 1L) H1 else H0
    sd_s <- if (s == 1L) seed1 else seed
    ws <- w * (if (s == 1L) P1 else 1 - P1)
    for (aq in seq_len(npar)) {
      for (bq in seq_len(npar)) {
        hessian <- hessian +
          crossprod(sd_s[[aq]], (ws * Hs[, aq, bq]) * sd_s[[bq]])
      }
    }
  }

  hessian <- (hessian + t(hessian)) / 2
  list(loglik = loglik, gradient = gradient, hessian = hessian)
}

#' @title Print a Marginal Break-Point Term
#' @name print.MarginalBreakTerm
#'
#' @description
#' Prints the label and the kind, how many latent break-points each group
#' carries, and, for a built term, over how many groups. A second line lists
#' the parameters, and a third names the prior where one was given.
#'
#' @details
#' The form is
#'
#' ```
#' <MarginalBreakTerm> 'jump' (jump): 1 latent break-point per group,
#'                     integrated out (3 groups)
#'   parameters: m1, tau1, delta1
#' ```
#'
#' The prior line appears only where `random(distrib = )` named a family; under
#' the default Gaussian there is nothing to name. A built structural term is
#' never described as "built", the group count being the tell:
#' [term_is_built()] tests for a design block, which this branch does not have.
#'
#' @param x A [MarginalBreakTerm()], built or not.
#' @param ... Unused, and accepted so that the signature matches [print()]'s.
#'
#' @return `x`, invisibly. Called for the lines it writes.
#'
#' @seealso [jump()], [term_params()].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(id = rep(1:3, each = 8), x = rep(1:8, 3))
#' dd$y <- rnorm(24, 2 * (dd$x >= 4.5), 0.4)
#'
#' # A specification, and the same term built over three groups.
#' jump(x, psi ~ random(~ 1 | id), marginal = TRUE)
#' term_build(jump(x, psi ~ random(~ 1 | id), marginal = TRUE), dd)
#'
#' @keywords internal
S7::method(print, MarginalBreakTerm) <- function(x, ...) {
  built <- length(x@blueprint) > 0L
  cat(sprintf(paste0("<MarginalBreakTerm> '%s' (%s): %d latent break-point%s",
                     " per group, integrated out%s\n"), x@label, x@kind,
              x@npsi, if (x@npsi > 1L) "s" else "",
              if (built) sprintf(" (%d groups)", length(x@blueprint$groups))
              else " (specification)"))
  cat("  parameters: ", paste(term_params(x), collapse = ", "), "\n",
      sep = "")
  if (!is.null(x@prior)) {
    cat("  prior: ", attr(S7::S7_class(x@prior), "name"), "\n", sep = "")
  }
  invisible(x)
}


#' @name term_simulate.MarginalBreakTerm
#'
#' @title Drawing Break-Points From Their Prior
#'
#' @description
#' One set of latent positions per group, drawn from the prior the term
#' declares, and the predictor each observation gets from them.
#'
#' @details
#' The latent positions ARE the model here, integrated out of the likelihood
#' and never estimated, so simulating from the model means
#' drawing them, once per group, and then evaluating the term at what was
#' drawn. Under the gaussian prior that is \eqn{N(m_k, \tau_k)}; under an
#' explicit prior it is a draw from that family with its location fixed at
#' zero, shifted by \eqn{m_1}, which is the same convention the likelihood
#' is written with.
#'
#' The shift is the term's own construction read at the drawn positions: a
#' change of level at each break-point for the step kind, a change of slope
#' for the continuous one, both for the joint one, and the linear term
#' beside them where the term carries it.
#'
#' The response is not drawn, the positions not reading it, so the
#' caller draws at the returned predictor.
#'
#' @param term A built [MarginalBreakTerm()].
#' @param psi The term's parameters, on the parameter scale.
#' @param eta The static part of the predictor.
#' @param draw Ignored: the positions do not read the response.
#' @param ... Ignored.
#'
#' @return A list with `eta`, `y` (`NULL`) and
#'   `latent`, a data frame of the drawn positions by group.
#'
#' @seealso [term_simulate()], [seg()]
#'
#' @keywords internal
S7::method(term_simulate, MarginalBreakTerm) <- function(term, psi, eta,
                                                         draw, ...) {
  bp <- .marg_built(term)
  v <- .marg_check_psi(term, psi)
  K <- term@npsi
  x <- bp$x
  out <- as.numeric(eta)
  if (term@linear) out <- out + v[["beta"]] * x
  rows <- list()
  for (gi in seq_along(bp$groups)) {
    rs <- bp$groups[[gi]]
    ps <- numeric(K)
    for (k in seq_len(K)) {
      ps[[k]] <- if (is.null(term@prior)) {
        stats::rnorm(1, v[[paste0("m", k)]], v[[paste0("tau", k)]])
      } else {
        # the prior's location is fixed at zero and m1 carries it, which is
        # the convention the likelihood is written with
        th <- as.list(v[term@prior@params])
        v[["m1"]] + as.numeric(
          distributions7::distrib_rng(term@prior, 1L, th))
      }
      if (term@kind %in% c("jump", "jseg")) {
        out[rs] <- out[rs] + v[[paste0("delta", k)]] * (x[rs] >= ps[[k]])
      }
      if (term@kind %in% c("seg", "jseg")) {
        out[rs] <- out[rs] +
          v[[paste0("gamma", k)]] * pmax(x[rs] - ps[[k]], 0)
      }
    }
    rows[[gi]] <- data.frame(group = names(bp$groups)[[gi]],
                             psi = seq_len(K), value = ps,
                             stringsAsFactors = FALSE)
  }
  list(eta = out, y = NULL, latent = do.call(rbind, rows))
}
