#' @include term_classes.R generics.R structural.R
NULL

#' @title S7 Class for Score-Driven Dynamics
#' @name GasTerm
#'
#' @description
#' The subclass of [structural_term()] holding a generalized autoregressive
#' score component: a time-varying level driven by the score of the
#' observation density, added to the predictor of one distribution parameter.
#' [gas()] constructs it. Its contribution is a state, so it implements
#' [term_filter()] and has no [term_matrix()] method at all.
#'
#' @details
#' # The seven properties of its own
#'
#' `p` is the number of score lags and `q` the number of autoregressive ones,
#' which together fix the parameter count at \eqn{1 + p + q} before any
#' subformula.
#'
#' `by` and `time` are the grouping and ordering expressions as written, kept
#' unevaluated; `NULL` means one series in row order.
#'
#' `links` holds whatever the caller overrode, empty where the defaults stand:
#' the identity on the level, the log on each loading, the rhobit on each
#' partial autocorrelation. `submodels` holds one right-hand side per
#' parameter developed over covariates.
#'
#' `blueprint` is filled by [term_build()] and carries the row order within
#' each group, the built sub-terms of each subformula and their designs. The
#' class overrides the branch's `blueprint` property because a structural term
#' has no design block to hang one on.
#'
#' # What the class is for
#'
#' The recursion, its exact Jacobian, the reverse pass, the curvature and the
#' contracted third derivative are all methods on it, so a fitting layer can
#' estimate the term's parameters beside the coefficients of every equation
#' and read the joint observed information.
#'
#' @inheritParams model_term
#' @param p The number of score lags, an integer of at least 0.
#' @param q The number of autoregressive lags, an integer of at least 0.
#' @param by An optional grouping expression; each group is filtered
#'   independently. `NULL` for one series.
#' @param time An optional ordering expression. `NULL` for row order.
#' @param links A named list of \pkg{linkfunctions7} links overriding the
#'   defaults, empty where none was given.
#' @param submodels A named list of one-sided formulas, one per parameter
#'   developed over covariates. Empty where none is.
#' @param blueprint A named list of the resolved ordering, grouping and
#'   sub-term designs, empty until [term_build()] fills it.
#'
#' @return An S7 object of class `GasTerm`, inheriting from
#'   [structural_term()] and [model_term()], with the seven properties above
#'   beside [model_term()]'s six.
#'
#' @seealso [gas()], the constructor; [term_filter()] for the recursion;
#'   [term_readable()] for the quantities a fitted one reports; [regime()] for
#'   the other dynamic term.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:60, y = c(rnorm(30), rnorm(30, 3)))
#'
#' tm <- gas(p = 1, q = 2, time = t)
#' S7::S7_inherits(tm, GasTerm)
#' c(p = tm@p, q = tm@q)
#'
#' # 1 + p + q parameters, and one chart each.
#' term_params(tm)
#' vapply(term_links(tm), function(l) l@link_name, character(1))
#'
#' # The build resolves the ordering; there is no block to read.
#' b <- term_build(tm, dd)
#' names(b@blueprint)
#' try(term_matrix(b))
#'
#' @export
GasTerm <- S7::new_class(
  name = "GasTerm",
  parent = structural_term,
  properties = list(
    p = S7::class_integer,
    q = S7::class_integer,
    by = S7::class_any,
    time = S7::class_any,
    links = S7::class_list,
    submodels = S7::class_list
  )
)

#' Score-Driven Dynamics
#'
#' @description
#' A generalized autoregressive score component (\cite{creal2013},
#' \cite{harvey2013}): a level \eqn{f_t} added to the linear predictor,
#' driven by the score of the observation density at the previous times,
#' \deqn{f_t = \omega + \sum_{i=1}^{p} a_i s_{t-i} +
#'   \sum_{j=1}^{q} b_j f_{t-j},}
#' with \eqn{s_t = \partial \ell_t / \partial \eta_t} the derivative of the
#' log-likelihood contribution with respect to the predictor it is
#' evaluated at.
#'
#' @details
#' The term adds no columns. The predictor at one time depends on the
#' data at the previous ones, so the contribution cannot be written as a
#' block, and [term_filter()] runs the recursion instead. That is what puts
#' it on the structural branch.
#'
#' What drives the recursion is the score of whatever distribution the
#' model carries, so the same term is a GARCH-like volatility model when
#' it enters the scale of a Gaussian, a dynamic count model when it enters
#' the mean of a Poisson, and a robust location filter when it enters a
#' Student t: a heavy-tailed score is bounded in the observation, so an
#' outlier moves the level by a bounded amount instead of in proportion to
#' its size.
#'
#' \subsection{The parameters and their chart}{
#' The parameters are the level \eqn{\omega}, the score loadings
#' \eqn{a_1, \dots, a_p}, and the persistence. Each is estimated on the
#' unconstrained scale of a link, and `links` overrides any of them;
#' the defaults are the following.
#'
#' The level carries the identity, being unconstrained. The loadings carry
#' the **log** link: a positive loading responds in the direction of the
#' score, which is the case the score-driven literature writes, and
#' positivity is then structural. A deviation or a submodel moves the loading
#' on the log scale, so no group and no observation can take a negative one.
#' A loading that must be free in sign is asked for with
#' `links = list(alpha1 = linkfunctions7::identity_link())`.
#'
#' The persistence is carried by **partial autocorrelations** rather
#' than by the coefficients \eqn{b_j}: the stationary region of an
#' autoregression is not a box, so no collection of scalar links covers
#' it, while the partial autocorrelations each range over \eqn{(-1, 1)}
#' independently and the Levinson-Durbin recursion carries them onto the
#' coefficients bijectively. At \eqn{q = 1} the two coincide. The
#' coordinate is named for the chart it lives on, `pacf1` and so on,
#' following the convention of \pkg{parameters7}.
#'
#' Whatever the links, a parameter modeled per group or per observation
#' stays in its own set: a departure acts on the unconstrained scale, so a
#' loading on the log link is positive and a persistence on the rhobit
#' link is stationary at every observation, whatever the departure is.
#' }
#'
#' \subsection{One parameter at a time}{
#' The level this term drives is a scalar, so it enters the predictor of
#' one distribution parameter. In the general formulation the level is a
#' vector with one entry per modeled parameter, \eqn{\omega} a vector and
#' the loadings \eqn{A_i} and \eqn{B_j} matrices, which lets the scale
#' respond to the score of the location and the other way round. The
#' recursion generalizes mechanically, and so does the derivative
#' propagated with it; what is missing is a way to say that one filter
#' spans several parameters, since a term written inside the formula of
#' one parameter has no place to declare it. That belongs to the model
#' layer. The persistence would also need a different chart: the
#' partial-autocorrelation construction below is a scalar one, and the
#' stationary region of a matrix autoregression is a bound on the
#' spectral radius of its companion matrix, which is not a box.
#'
#' The score driving the recursion is used unscaled. The general
#' formulation carries a scaling matrix, usually an inverse information,
#' which the curvature this term already receives would supply.
#' }
#'
#' \subsection{Groups and time}{
#' `by` filters each group independently, as a panel of short series needs,
#' and `time` gives the order within a group.
#' Without `time` the rows are taken in the order they appear.
#' }
#'
#' \subsection{A parameter developed with covariates}{
#' A two-sided formula in `...` whose left side names a parameter
#' develops it as \eqn{\psi_{j,t} = g_j^{-1}(z_t^\top\gamma_j)}, the
#' design \eqn{Z} built from the right-hand side through
#' [interpret_formula()], so it takes any additive term of the
#' package:
#' \preformatted{gas(p = 1, q = 1, omega ~ ridge(~g), alpha1 ~ s(x),
#'     pacf1 ~ random(~1 | id), by = id)}
#' The development acts on the unconstrained scale of the parameter's own
#' link. That is what keeps every per-observation value inside the
#' parameter's own set whatever the coefficients are: a loading on the
#' log link is positive at every observation, a persistence on the rhobit
#' chart is inside \eqn{(-1, 1)} at every observation, and at \eqn{q = 1}
#' that bounds the recursion's growth step by step. The coefficients
#' \eqn{\gamma_j} are the term's parameters, unconstrained and on the
#' identity link; the penalties the sub-terms carry are reported through
#' [term_penalties()] under the key `parameter::subterm`.
#' A parameter that varies by observation changes the recursion itself,
#' \deqn{f_t = \omega_t + \sum_i a_{i,t}\, s_{t-i}
#'   + \sum_j b_{j,t}\, f_{t-j},}
#' with \eqn{b_t} from the Levinson-Durbin map of that observation's
#' partial autocorrelations, and the filter, its derivative, the reverse
#' pass and the curvature all run the general recursion.
#'
#' `by = ~f` (a formula, where a grouping variable is a bare symbol)
#' is the shorthand giving the same subformula to every parameter; mixing
#' it with per-parameter formulas is an error. A structural term, and a
#' term whose block moves with its own coefficients, are rejected inside
#' a subformula.
#' }
#'
#' \subsection{A population value and a departure per group}{
#' The panel case is one subformula:
#' `gas(omega ~ random(~1 | id), by = id)` is a population value (the
#' intercept of the development) plus one unconstrained departure per
#' group, shrunk by the random intercept's own ridge, whose hyperparameter
#' a fitting layer estimates. `lasso(...)` in the subformula sets the
#' departures of the groups that do not need one exactly to zero.
#'
#' The departures act on the unconstrained scale of the parameter's chart,
#' so a persistence stays inside \eqn{(-1, 1)} whatever the departure is.
#' Their penalty is also what identifies them: a population value and
#' \eqn{m} departures are \eqn{m+1} numbers describing \eqn{m} group
#' values, so without the penalty the likelihood is flat along one
#' direction per developed parameter. This is the parametrization of a
#' random effect and it is identified the same way, there by a variance
#' component, here by the penalty. An unpenalized development over group
#' indicators is therefore for reading a filter at given parameters rather
#' than for fitting one.
#'
#' Earlier releases spelled this case as `deviations =` with a
#' `penalty =`; both arguments are gone, the subformula reproducing
#' them exactly (the same fit to the printed digit, hyperparameter
#' included) while covering what they could not.
#' }
#'
#' @param p The number of score lags. Defaults to 1.
#' @param q The number of autoregressive lags. Defaults to 1.
#' @param ... Two-sided formulas whose left side names a parameter, one
#'   per parameter to be developed with covariates, e.g.
#'   `alpha1 ~ s(x)`; see the section above.
#' @param by An optional grouping variable, evaluated in the data; each
#'   group is filtered independently, from its own starting level. A
#'   formula here is the shorthand giving the same subformula to every
#'   parameter, and then no grouping is implied.
#' @param time An optional ordering variable, evaluated in the data.
#' @param links An optional named list of \pkg{linkfunctions7} links over
#'   the parameters of [term_params()], overriding the defaults
#'   described above. A deviation cannot be named: it is unconstrained by
#'   construction, acting on the scale its parameter's own link defines.
#' @param label A single non-empty string naming the term.
#'
#' @return An object of class [GasTerm()] (a specification; see
#'   [term_build()]).
#'
#' @references
#' Creal, D., Koopman, S. J. and Lucas, A. (2013). Generalized
#' autoregressive score models with applications. *Journal of Applied
#' Econometrics*, 28(5), 777--795.
#'
#' Harvey, A. C. (2013). *Dynamic Models for Volatility and Heavy
#' Tails*. Cambridge University Press.
#'
#' @seealso [regime()] for the other dynamic term, [term_filter()] for the
#'   recursion, [term_params()] and [term_links()] for the parameters and
#'   their charts, [term_readable()] for what a fitted one reports.
#'
#' @examples
#' # A level, one loading and one persistence coordinate.
#' term_params(gas(p = 1, q = 1))
#' term_params(gas(p = 2, q = 2))
#'
#' # The loading rides a log chart so that it stays positive, and the
#' # persistence a rhobit so that the filter stays stationary.
#' vapply(term_links(gas(p = 1, q = 2)), function(l) l@link_name,
#'        character(1))
#'
#' # Running the filter: the term supplies the recursion and the caller
#' # supplies the density's score and curvature. Here a Gaussian mean.
#' set.seed(1)
#' dd <- data.frame(t = 1:60, y = c(rnorm(30), rnorm(30, 3)))
#' b <- term_build(gas(p = 1, q = 1, time = t), dd)
#' out <- term_filter(b, eta = rep(0, 60), y = dd$y,
#'                    score = function(e, i) dd$y[i] - e,
#'                    curvature = function(e, i) -1,
#'                    psi = list(omega = 0.1, alpha1 = 0.3, pacf1 = 0.8))
#'
#' # The level tracks the change in the mean.
#' round(out$eta[c(1, 15, 30, 32, 45, 60)], 3)
#'
#' # And the Jacobian it propagates is exact, not differenced: a central
#' # difference of the filter agrees with it to the step's own accuracy.
#' f <- function(v) sum(term_filter(b, rep(0, 60), dd$y,
#'                                  function(e, i) dd$y[i] - e,
#'                                  function(e, i) -1,
#'                                  list(omega = v[1], alpha1 = v[2],
#'                                       pacf1 = v[3]))$eta)
#' v0 <- c(0.1, 0.3, 0.8)
#' h  <- 1e-5
#' fd <- vapply(1:3, function(j) {
#'   e <- numeric(3); e[j] <- h
#'   (f(v0 + e) - f(v0 - e)) / (2 * h)
#' }, numeric(1))
#' max(abs(colSums(out$jacobian) - fd))
#'
#' # A panel: each group filtered from its own starting level.
#' dd$id <- rep(1:3, each = 20)
#' term_build(gas(p = 1, q = 1, by = id, time = t), dd)
#'
#' # A developed parameter expands in place, and brings its sub-term's
#' # penalty with it.
#' dd$g <- factor(rep(c("u", "v"), 30))
#' gb <- term_build(gas(p = 1, q = 1, omega ~ ridge(~ g), time = t), dd)
#' term_params(gb)
#' vapply(term_penalties(gb), function(e) e$name, character(1))
#'
#'
#' # Fitted. The data are simulated from a known truth, so the
#' # estimates below can be read against it.
#' if (requireNamespace("statmodels7", quietly = TRUE)) {
#'   # drawn from the recursion the term implements, with a gaussian mean
#'   set.seed(2)
#'   om <- 0.05; al <- 0.3; be <- 0.9
#'   lev <- om / (1 - be)
#'   yy <- numeric(600)
#'   for (i in seq_len(600)) {
#'     yy[i] <- lev + rnorm(1)
#'     lev <- om + al * (yy[i] - lev) + be * lev
#'   }
#'   fd <- data.frame(t = seq_len(600), y = yy)
#'   ft <- statmodels7::statmod(y ~ 0 + gas(p = 1, q = 1, time = t),
#'                              distributions7::gaussian1_distrib(), fd)
#'   # truth: omega 0.05, alpha 0.3, beta 0.9. beta is reported as the
#'   # autoregressive coefficient, which its coordinate is not.
#'   round(coef(ft)$mu, 3)
#' }
#' @export
gas <- function(p = 1, q = 1, ..., by = NULL, time = NULL,
                links = NULL, label = "gas") {
  chk <- function(v, nm, lo) {
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < lo ||
        v != round(v)) {
      stop(sprintf("'%s' must be a single integer of at least %d.", nm, lo),
           call. = FALSE)
    }
    as.integer(v)
  }
  p <- chk(p, "p", 1L)
  q <- chk(q, "q", 0L)
  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !nzchar(label)) {
    stop("'label' must be a single non-empty character string.",
         call. = FALSE)
  }
  by <- substitute(by)
  base <- .gas_base_params(p, q)
  submodels <- .gas_gather_submodels(list(...), base)
  # `by = ~f` is the shorthand giving every parameter the same subformula;
  # a grouping variable is a bare symbol, so the two are told apart by shape
  if (is.call(by) && identical(by[[1L]], as.name("~"))) {
    if (length(submodels)) {
      stop(paste("'by = ~f' gives every parameter the same subformula;",
                 "mixing it with per-parameter formulas is an error."),
           call. = FALSE)
    }
    if (length(by) != 2L) {
      stop("a formula 'by' must be one-sided, e.g. by = ~ridge(~g).",
           call. = FALSE)
    }
    f <- stats::as.formula(by, env = parent.frame())
    submodels <- stats::setNames(rep(list(f), length(base)), base)
    by <- NULL
  }
  if (!is.null(links)) {
    if (!is.list(links) || is.null(names(links)) || anyNA(names(links)) ||
        !all(nzchar(names(links)))) {
      stop("'links' must be a named list of linkfunctions7 links.",
           call. = FALSE)
    }
    bad <- setdiff(names(links), base)
    if (length(bad)) {
      stop(sprintf("'links' names '%s'; the parameters are %s.",
                   bad[1L], paste(base, collapse = ", ")), call. = FALSE)
    }
    for (lname in names(links)) {
      if (!S7::S7_inherits(links[[lname]], linkfunctions7::link)) {
        stop(sprintf("'links$%s' is not a linkfunctions7 link.", lname),
             call. = FALSE)
      }
    }
  }
  GasTerm(label = label, p = p, q = q,
          by = by, time = substitute(time),
          links = if (is.null(links)) list() else links,
          submodels = submodels,
          blueprint = list())
}

# the per-parameter subformulas of `...`: two-sided formulas whose left
# side names a base parameter, at most one per parameter
.gas_gather_submodels <- function(dots, base) {
  # a named entry is an argument the signature does not take -- the shape a
  # mistyped or removed argument arrives in -- and it is reported by name
  # rather than swallowed by the formula check below
  nms <- names(dots)
  if (!is.null(nms) && any(nzchar(nms))) {
    stop(sprintf(paste("unused argument '%s': '...' takes two-sided",
                       "formulas whose left side names a parameter."),
                 nms[nzchar(nms)][1L]), call. = FALSE)
  }
  out <- list()
  for (d in dots) {
    if (!inherits(d, "formula") || length(d) != 3L || !is.name(d[[2L]])) {
      stop(paste("arguments in '...' must be two-sided formulas whose left",
                 "side names a parameter, e.g. alpha1 ~ s(x)."),
           call. = FALSE)
    }
    pn <- as.character(d[[2L]])
    if (!(pn %in% base)) {
      stop(sprintf("the subformula names '%s'; the parameters are %s.",
                   pn, paste(base, collapse = ", ")), call. = FALSE)
    }
    if (!is.null(out[[pn]])) {
      stop(sprintf("parameter '%s' carries two subformulas.", pn),
           call. = FALSE)
    }
    out[[pn]] <- stats::as.formula(call("~", d[[3L]]), env = environment(d))
  }
  out
}

# the parameters of the filter itself, before any deviation
.gas_base_params <- function(p, q) {
  # The level and the score loadings carry the names the score-driven
  # literature gives them, spelled out: they are the quantities themselves,
  # each reported through its own link, so a name can promise what it
  # reports. The persistence cannot be named `beta`: it rides a partial
  # autocorrelation, the stationary region not being a box, and a free
  # coordinate named after the coefficient would promise the coefficient and
  # report the chart. The coefficient is what a fitted model REPORTS,
  # through `term_readable()`.
  c("omega",
    if (p > 0L) paste0("alpha", seq_len(p)),
    if (q > 0L) paste0("pacf", seq_len(q)))
}

#' @title The Parameters of a Score-Driven Term
#' @name term_params.GasTerm
#'
#' @description
#' `"omega"`, then `"alpha1"` ... `"alphap"`, then `"pacf1"` ... `"pacfq"`:
#' the level, one loading per score lag, and one partial autocorrelation per
#' autoregressive lag. A parameter carrying a subformula is **expanded in
#' place** into its coefficients, named `parameter.coefficient`.
#'
#' @details
#' The persistence coordinates are named for the chart they live on, never
#' for the quantity a reader reads. `pacf1` is a partial autocorrelation;
#' the autoregressive coefficient \eqn{\beta_1} the literature writes is a
#' function of the whole chart through Levinson-Durbin, and coincides with the
#' coordinate only at \eqn{q = 1}. [term_readable()] is what carries the
#' coordinates onto the coefficients.
#'
#' The order is the one [term_links()], [term_start()], [term_filter()]'s
#' Jacobian columns and the joint variance matrix are all indexed by.
#'
#' @param term A [GasTerm()]. Unbuilt, the subformulas have not been resolved,
#'   so a developed parameter still appears as itself.
#' @param ... Unused.
#'
#' @return A character vector of length [term_npar()].
#'
#' @seealso [term_links()] for the charts, [term_readable()] for the
#'   quantities, [gas()] for what they mean.
#'
#' @examples
#' term_params(gas(p = 1, q = 1))
#' term_params(gas(p = 2, q = 3))
#'
#' # A subformula expands its parameter in place.
#' set.seed(1)
#' dd <- data.frame(t = 1:40, y = rnorm(40), g = factor(rep(c("u", "v"), 20)))
#' term_params(term_build(gas(p = 1, q = 1, omega ~ g, time = t), dd))
#'
#' @keywords internal
S7::method(term_params, GasTerm) <- function(term, ...) {
  if (length(term@blueprint) && !is.null(term@blueprint$sub)) {
    return(.gas_sub_layout(term)$names)
  }
  .gas_base_params(term@p, term@q)
}

#' @title The Level of a Score-Driven Term
#' @name term_level_param.GasTerm
#' @description
#' `"omega"`, which adds a constant to the equation's predictor and is
#' therefore the direction an intercept there also spans. With the level
#' developed by a subformula, the coordinates whose design column is
#' constant, each of which shifts the predictor the same way.
#' @param term A [GasTerm()].
#' @param ... Unused.
#' @return A character vector, usually of length one.
#' @keywords internal
S7::method(term_level_param, GasTerm) <- function(term, ...) {
  sub <- if (length(term@blueprint)) term@blueprint$sub else NULL
  if (is.null(sub) || is.null(sub$omega)) return("omega")
  # with the level developed, the coordinates that shift the predictor by
  # a constant are the ones whose column is constant; any residual
  # constant hiding in a combination of non-constant columns is the
  # layer's subspace question, not the term's
  Z <- as.matrix(sub$omega$Z)
  const <- vapply(seq_len(ncol(Z)), function(k) {
    v <- Z[, k]
    all(v == v[1L])
  }, logical(1))
  term_params(term)[.gas_sub_layout(term)$idx$omega][const]
}

#' @title The Design of a Developed Level
#' @name term_level_design.GasTerm
#' @description
#' The design of `omega`'s development when it carries one, its
#' columns named after the coordinates, so a fitting layer can compare its
#' span with the equation's; `NULL` for a scalar level.
#' @param term A built [GasTerm()].
#' @param ... Unused.
#' @return A numeric matrix or `NULL`.
#' @keywords internal
S7::method(term_level_design, GasTerm) <- function(term, ...) {
  sub <- if (length(term@blueprint)) term@blueprint$sub else NULL
  if (is.null(sub) || is.null(sub$omega)) return(NULL)
  Z <- as.matrix(sub$omega$Z)
  colnames(Z) <- term_params(term)[.gas_sub_layout(term)$idx$omega]
  Z
}

#' @title What a Fitted Score-Driven Term Reports
#' @name term_readable.GasTerm
#' @description
#' The level, the score loadings and the autoregressive coefficients of the
#' literature, `omega`, `alpha1` and `beta1`, with the
#' Jacobian from the term's own parameters.
#' @details
#' The level and the loadings are reported through their own links, each a
#' function of its own coordinate alone, which the base method already
#' does. The persistence is not: it is carried on a partial
#' autocorrelation, and the coefficients come from the Levinson-Durbin
#' recursion, whose Jacobian the term already computes for the filter.
#' Chained onto the rhobit link of each coordinate, that Jacobian is what
#' a delta-method standard error for \eqn{\beta_j} needs. At \eqn{q = 1}
#' the two coincide and the chain factor is the link's alone; above it
#' they do not.
#'
#' A deviation is reported as it stands, being unconstrained and defined on
#' the scale of the parameter it departs from.
#' @param term A [GasTerm()].
#' @param zeta The parameters on the unconstrained scale.
#' @param ... Unused.
#' @return A list, as [term_readable()] documents.
#' @keywords internal
S7::method(term_readable, GasTerm) <- function(term, zeta, ...) {
  nm <- term_params(term)
  links <- term_links(term)
  z <- unlist(zeta[nm])
  q <- term@q
  base <- .gas_base_params(term@p, q)
  out <- S7::method(term_readable, model_term)(term, zeta)
  if (q == 0L) return(out)

  i_pa <- match(paste0("pacf", seq_len(q)), nm)
  # with any partial autocorrelation developed, the coefficients vary by
  # observation and there is no single beta to report: the coordinates are
  # reported as the base method gives them
  if (anyNA(i_pa)) return(out)
  lk <- links[[nm[i_pa[1L]]]]
  rho <- linkfunctions7::linkinv(lk, z[i_pa])
  k1 <- linkfunctions7::dlinkinv(lk, z[i_pa])
  ld <- gas_levinson(rho)
  out$name[i_pa] <- paste0("beta", seq_len(q))
  out$value[i_pa] <- ld$phi
  # the whole chart reaches every coefficient: row j of the recursion's
  # jacobian scaled by the link's derivative, and NOT the diagonal entry the
  # base method placed there
  out$jacobian[i_pa, ] <- 0
  out$jacobian[i_pa, i_pa] <- ld$jacobian * rep(k1, each = q)
  rownames(out$jacobian)[i_pa] <- out$name[i_pa]
  # a coefficient of a stationary autoregression is not confined to an
  # interval a scalar link expresses, the region not being a box, so its
  # interval is built on the identity scale
  out$scale[i_pa] <- rep(list(linkfunctions7::identity_link()), q)
  names(out$scale)[i_pa] <- out$name[i_pa]
  out
}

#' @title The Charts of a Score-Driven Term's Parameters
#' @name term_links.GasTerm
#'
#' @description
#' The identity on the level, the **log** on every loading and the **rhobit**
#' on every partial autocorrelation, unless the `links` argument of [gas()]
#' overrode one. Each carries its parameter's own admissible set onto the whole
#' real line, so an optimizer proposing anything at all gets an admissible
#' filter.
#'
#' @details
#' The log on a loading makes positivity structural: a subformula develops the
#' loading on that scale, so no group and no observation can take a negative
#' one. A loading that must be free in sign is asked for by name.
#'
#' The rhobit on a partial autocorrelation keeps it inside \eqn{(-1, 1)}, and
#' Levinson-Durbin then carries the whole chart onto a stationary
#' autoregression. That is the reason for the chart: the stationary region in
#' the coefficients is not a box, so no collection of scalar links covers it.
#'
#' Where a parameter carries a subformula the link is applied **inside** the
#' development, so the parameter is admissible at every observation and its
#' coefficients are unconstrained on the identity.
#'
#' @param term A [GasTerm()].
#' @param ... Unused.
#'
#' @return A named list of \pkg{linkfunctions7} links, one per entry of
#'   [term_params()].
#'
#' @seealso [term_params()], [gas()], [term_start()].
#'
#' @examples
#' vapply(term_links(gas(p = 2, q = 2)), function(l) l@link_name, character(1))
#'
#' # Each carries its own set onto the line, so nothing is out of range.
#' lk <- term_links(gas(p = 1, q = 1))
#' vapply(lk, function(l) paste(l@link_bounds, collapse = ", "), character(1))
#'
#' # Overridden: a loading free in sign.
#' vapply(term_links(gas(p = 1, q = 1,
#'                       links = list(alpha1 = linkfunctions7::identity_link()))),
#'        function(l) l@link_name, character(1))
#'
#' @keywords internal
S7::method(term_links, GasTerm) <- function(term, ...) {
  base <- .gas_base_params(term@p, term@q)
  nm <- term_params(term)
  stats::setNames(lapply(nm, function(p) {
    # a deviation, and a coordinate of a development, are unconstrained
    # already: they act on the scale the parameter's own link carries it
    # to, which the term applies inside
    if (!(p %in% base)) return(linkfunctions7::identity_link())
    .gas_param_link(term, p)
  }), nm)
}

#' @title Where a Score-Driven Term's Parameters Start
#' @name term_start.GasTerm
#' @description
#' Zero on the unconstrained scale of every parameter except the score
#' loadings, which start at \eqn{0.1} on the parameter scale, through
#' whatever link each one carries.
#' @details
#' Zero is the natural point of every other chart: a level of zero, no
#' persistence, no deviation. The term therefore starts as near the model
#' without it as its charts allow. The loadings are the exception because
#' zero on the log scale is a loading of one, a response strong enough to
#' destabilize the recursion at ordinary curvatures; \eqn{0.1} is a weak
#' response, and it is applied on the parameter scale so the start means
#' the same thing whatever chart a loading rides.
#' @param term A [GasTerm()].
#' @param ... Unused.
#' @return A named numeric vector on the unconstrained scale.
#' @keywords internal
S7::method(term_start, GasTerm) <- function(term, ...) {
  nm <- term_params(term)
  z <- stats::setNames(numeric(length(nm)), nm)
  base <- .gas_base_params(term@p, term@q)
  target <- function(j) {
    if (startsWith(j, "alpha")) {
      linkfunctions7::linkfun(.gas_param_link(term, j), 0.1)
    } else 0
  }
  sub <- if (length(term@blueprint)) term@blueprint$sub else NULL
  if (!is.null(sub)) {
    lay <- .gas_sub_layout(term)
    for (j in base) {
      z[lay$idx[[j]]] <- if (is.null(sub[[j]])) target(j) else {
        # the constant start projected onto the development's design
        .seg_proj_start(sub[[j]]$Z, target(j))
      }
    }
    return(z)
  }
  for (j in base) z[[j]] <- target(j)
  z
}

#' @title Penalties of a Score-Driven Term
#' @name term_penalties.GasTerm
#' @description
#' The penalties the subformulas' sub-terms declare, each under the key
#' `parameter::subterm` and covering that sub-term's coefficients in
#' the term's own numbering. The scalar parameters are unpenalized, and
#' the list is empty for a specification, whose developments do not exist
#' until the term is built.
#' @param term A built [GasTerm()].
#' @param ... Unused.
#' @return A list of entries, as [term_penalties()] documents.
#' @keywords internal
S7::method(term_penalties, GasTerm) <- function(term, ...) {
  sub <- if (length(term@blueprint)) term@blueprint$sub else NULL
  if (is.null(sub)) return(list())
  lay <- .gas_sub_layout(term)
  out <- list()
  for (j in names(sub)) {
    for (e in sub[[j]]$penalties) {
      out[[length(out) + 1L]] <- list(
        name = paste0(j, "::", e$name),
        index = lay$idx[[j]][e$index],
        penalty = e$penalty,
        # what the sub-term holds travels WITH the entry, so a structural
        # term propagates it by copying the entry and needs to know nothing
        # about hyperparameters
        fixed = e$fixed, n_values = e$n_values, values = e$values,
        min_ratio = e$min_ratio, search = e$search, ids = e$ids)
    }
  }
  out
}

#' The Autoregressive Coefficients Behind the Partial Autocorrelations
#'
#' @description
#' The Levinson-Durbin recursion carrying partial autocorrelations onto
#' the coefficients of a stationary autoregression, with the Jacobian of
#' that map propagated alongside.
#'
#' @param pacf A numeric vector of partial autocorrelations in
#'   \eqn{(-1, 1)}.
#'
#' @return A list with `phi`, the coefficients, and `jacobian`,
#'   the matrix of their derivatives with respect to `pacf`.
#'
#' @keywords internal
gas_levinson <- function(pacf) {
  q <- length(pacf)
  if (q == 0L) return(list(phi = numeric(0), jacobian = matrix(0, 0, 0)))
  phi <- numeric(0)
  jac <- matrix(0, 0, q)
  for (k in seq_len(q)) {
    new <- numeric(k)
    njac <- matrix(0, k, q)
    new[k] <- pacf[k]
    njac[k, k] <- 1
    if (k > 1L) {
      rev_idx <- rev(seq_len(k - 1L))
      new[seq_len(k - 1L)] <- phi - pacf[k] * phi[rev_idx]
      njac[seq_len(k - 1L), ] <- jac - pacf[k] * jac[rev_idx, , drop = FALSE]
      njac[seq_len(k - 1L), k] <- njac[seq_len(k - 1L), k] - phi[rev_idx]
    }
    phi <- new
    jac <- njac
  }
  list(phi = phi, jacobian = jac)
}

#' The Second Derivative of the Levinson-Durbin Map
#'
#' @description
#' [gas_levinson()] with the second derivatives of the
#' coefficients in the partial autocorrelations propagated as well.
#'
#' @details
#' The recursion is
#' \deqn{\phi^{(k)}_k = \rho_k, \qquad
#'   \phi^{(k)}_i = \phi^{(k-1)}_i - \rho_k\phi^{(k-1)}_{k-i},}
#' which is bilinear: \eqn{\rho_k} multiplies quantities that do not depend
#' on it. Differentiating twice therefore adds no new kind of term, only the
#' two places the product rule puts the first derivative,
#' \deqn{H^{(k)}_i = H^{(k-1)}_i - \rho_k H^{(k-1)}_{k-i}
#'   - e_k \left(J^{(k-1)}_{k-i}\right)^{\!\top}
#'   - J^{(k-1)}_{k-i} e_k^{\top},}
#' and the last coefficient's second derivative is zero at every order, it
#' being \eqn{\rho_k} itself.
#'
#' It is wanted because the observed information of a model carrying this
#' term needs the second derivative of the predictor in the term's own
#' parameters, and the persistence reaches the predictor through this map.
#'
#' @param pacf A numeric vector of partial autocorrelations in
#'   \eqn{(-1, 1)}.
#'
#' @return A list with `phi`, `jacobian` and `hessian`, the
#'   last a list of one symmetric matrix per coefficient.
#'
#' @seealso [gas_levinson()]
#'
#' @keywords internal
gas_levinson2 <- function(pacf) {
  q <- length(pacf)
  if (q == 0L) {
    return(list(phi = numeric(0), jacobian = matrix(0, 0, 0),
                hessian = list()))
  }
  phi <- numeric(0)
  jac <- matrix(0, 0, q)
  hes <- list()
  for (k in seq_len(q)) {
    new <- numeric(k)
    njac <- matrix(0, k, q)
    nhes <- replicate(k, matrix(0, q, q), simplify = FALSE)
    new[k] <- pacf[k]
    njac[k, k] <- 1
    if (k > 1L) {
      rev_idx <- rev(seq_len(k - 1L))
      new[seq_len(k - 1L)] <- phi - pacf[k] * phi[rev_idx]
      njac[seq_len(k - 1L), ] <- jac - pacf[k] * jac[rev_idx, , drop = FALSE]
      njac[seq_len(k - 1L), k] <- njac[seq_len(k - 1L), k] - phi[rev_idx]
      for (i in seq_len(k - 1L)) {
        r <- rev_idx[i]
        h <- hes[[i]] - pacf[k] * hes[[r]]
        h[k, ] <- h[k, ] - jac[r, ]
        h[, k] <- h[, k] - jac[r, ]
        nhes[[i]] <- h
      }
    }
    phi <- new
    jac <- njac
    hes <- nhes
  }
  list(phi = phi, jacobian = jac, hessian = hes)
}

#' The Third Derivative of the Levinson-Durbin Map, in One Direction
#'
#' @description
#' [gas_levinson2()]'s second derivatives differentiated once more
#' and contracted against a single direction, one matrix per coefficient.
#'
#' @details
#' The exact gradient of a marginal criterion over a penalty on this term's
#' own parameters needs the third derivative of the predictor, and the
#' persistence reaches the predictor through this map. It is needed only
#' contracted: the criterion asks for \eqn{\mathrm{tr}(M\,\partial K/\partial
#' u[v])}, a derivative along the single direction the penalized mode moves
#' in, so what is propagated is a matrix per coefficient and never a
#' three-index array.
#'
#' Differentiating the hessian recursion of [gas_levinson2()] once
#' more along \eqn{w} adds no new kind of term, the map being bilinear:
#' \deqn{T^{(k)}_i = T^{(k-1)}_i - \rho_k T^{(k-1)}_{k-i}
#'   - w_k H^{(k-1)}_{k-i}
#'   - e_k\left(H^{(k-1)}_{k-i}w\right)^{\!\top}
#'   - \left(H^{(k-1)}_{k-i}w\right)e_k^{\top},}
#' and the last coefficient's third derivative is zero at every order, it
#' being \eqn{\rho_k} itself.
#'
#' The map is multilinear of degree \eqn{k} in the first \eqn{k} partial
#' autocorrelations, so the result is identically zero for \eqn{q \le 2}:
#' at \eqn{q = 2} the only non-trivial coefficient is
#' \eqn{\phi_1 = \rho_1(1-\rho_2)}, which is bilinear. A check of this
#' function that stops at \eqn{q = 2} compares zero with zero and asserts
#' nothing.
#'
#' @param pacf A numeric vector of partial autocorrelations in
#'   \eqn{(-1, 1)}.
#' @param w The direction to contract against, as long as `pacf`.
#'
#' @return A list of one `q` by `q` matrix per coefficient.
#'
#' @seealso [gas_levinson2()]
#'
#' @keywords internal
gas_levinson3 <- function(pacf, w) {
  q <- length(pacf)
  if (q == 0L) return(list())
  phi <- numeric(0)
  jac <- matrix(0, 0, q)
  hes <- list()
  thi <- list()
  for (k in seq_len(q)) {
    new <- numeric(k)
    njac <- matrix(0, k, q)
    nhes <- replicate(k, matrix(0, q, q), simplify = FALSE)
    nthi <- replicate(k, matrix(0, q, q), simplify = FALSE)
    new[k] <- pacf[k]
    njac[k, k] <- 1
    if (k > 1L) {
      rev_idx <- rev(seq_len(k - 1L))
      new[seq_len(k - 1L)] <- phi - pacf[k] * phi[rev_idx]
      njac[seq_len(k - 1L), ] <- jac - pacf[k] * jac[rev_idx, , drop = FALSE]
      njac[seq_len(k - 1L), k] <- njac[seq_len(k - 1L), k] - phi[rev_idx]
      for (i in seq_len(k - 1L)) {
        r <- rev_idx[i]
        h <- hes[[i]] - pacf[k] * hes[[r]]
        h[k, ] <- h[k, ] - jac[r, ]
        h[, k] <- h[, k] - jac[r, ]
        nhes[[i]] <- h
        hw <- as.numeric(hes[[r]] %*% w)
        tt <- thi[[i]] - pacf[k] * thi[[r]] - w[[k]] * hes[[r]]
        tt[k, ] <- tt[k, ] - hw
        tt[, k] <- tt[, k] - hw
        nthi[[i]] <- tt
      }
    }
    phi <- new
    jac <- njac
    hes <- nhes
    thi <- nthi
  }
  thi
}

#' @title Build a Score-Driven Term
#' @name term_build.GasTerm
#'
#' @description
#' Resolves the grouping, the ordering and every subformula against the data,
#' builds the sub-terms of the developments and records all of it in the
#' blueprint. The recursion itself runs later, in [term_filter()].
#'
#' @details
#' `by` and `time` are evaluated in the data. Each must give one non-missing
#' value per row, and the blueprint records the row indices of each group in
#' time order, and the filter iterates over exactly that.
#'
#' Each subformula's right-hand side goes through [interpret_formula()] and its
#' terms are built, so their blueprints are recorded and a prediction reapplies
#' them. A structural sub-term, and one whose own block moves with its
#' coefficients, are rejected: a parameter's development must be a fixed
#' design. The penalties those sub-terms carry become the term's own
#' [term_penalties()] entries, keyed `parameter::subterm`.
#'
#' Without `time` the rows are taken as they come. That is a real choice: the
#' recursion is about order, so a data frame that is not already sorted gives a
#' different model.
#'
#' @param term A [GasTerm()].
#' @param data A data frame carrying whatever `by`, `time` and the subformulas
#'   name.
#' @param ... Unused.
#'
#' @return The term with `blueprint` filled. [term_is_built()] reads that
#'   property on this branch, so it is `TRUE` for the result.
#'
#' @seealso [gas()], [term_filter()], [term_penalties()].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:60, y = rnorm(60), id = rep(1:3, each = 20),
#'                  g = factor(rep(c("u", "v"), 30)))
#'
#' # One series, then a panel of three.
#' lengths(term_build(gas(p = 1, q = 1, time = t), dd)@blueprint$order)
#' lengths(term_build(gas(p = 1, q = 1, by = id, time = t), dd)@blueprint$order)
#'
#' # A development expands the parameter and brings its penalty.
#' gb <- term_build(gas(p = 1, q = 1, omega ~ ridge(~ g), time = t), dd)
#' term_params(gb)
#' vapply(term_penalties(gb), function(e) e$name, character(1))
#'
#' @keywords internal
S7::method(term_build, GasTerm) <- function(term, data, ...) {
  n <- nrow(data)
  gf <- if (is.null(term@by)) factor(rep(1L, n)) else {
    factor(eval(term@by, data, baseenv()))
  }
  grp <- as.integer(gf)
  tm <- if (is.null(term@time)) seq_len(n) else {
    eval(term@time, data, baseenv())
  }
  if (length(grp) != n || length(tm) != n) {
    stop("'by' and 'time' must evaluate to one value per row.", call. = FALSE)
  }
  if (anyNA(grp) || anyNA(tm)) {
    stop("'by' and 'time' must not contain missing values.", call. = FALSE)
  }
  # the rows of each group, in time order: the filter walks these and
  # scatters its result back to the original positions
  ord <- split(order(grp, tm), grp[order(grp, tm)])
  # the times are kept beside the ordering because a CONTINUATION of the
  # series has to know where the observed part ends: a new row belongs after
  # the last time of its own group, and nothing else records that
  term@blueprint <- list(order = ord, n = n, levels = levels(gf),
                         times = tm, group = grp,
                         by = term@by, time = term@time)
  # the developments, each through the interpreter exactly as a parameter
  # of nl(): the built sub-terms are kept, prediction reapplying each
  # one's own blueprint, and their penalties are collected at
  # term_penalties() under the key parameter::subterm
  if (length(term@submodels)) {
    sub <- list()
    for (pn in names(term@submodels)) {
      sub[[pn]] <- .nl_submodel(pn, term@submodels[[pn]], data)
    }
    term@blueprint$sub <- sub
  }
  # a penalty supplied as an object covers a fixed number of coefficients
  # and a development's width is counted here for the first time, so this
  # is the earliest point at which the two can be compared, and the latest
  # at which a caller can still read the mistake against what they wrote
  term_penalties(term)
  term
}

# The coefficients the recursion runs on, and the constants of the chart,
# from one set of parameter values: the score loadings as given, the
# autoregressive coefficients through Levinson-Durbin with its jacobian,
# and the level the recursion starts from.
.gas_coefs <- function(psi, p, q) {
  np <- length(psi)
  omega <- psi[[1L]]
  a <- if (p > 0L) psi[1L + seq_len(p)] else numeric(0)
  ld <- gas_levinson(if (q > 0L) psi[1L + p + seq_len(q)] else numeric(0))
  b <- ld$phi
  i_a <- if (p > 0L) 1L + seq_len(p) else integer(0)
  i_pa <- if (q > 0L) 1L + p + seq_len(q) else integer(0)
  db <- matrix(0, max(q, 1L), np)
  if (q > 0L) db[seq_len(q), i_pa] <- ld$jacobian

  sb <- if (q > 0L) sum(b) else 0
  if (abs(1 - sb) < 1e-10) {
    stop("the autoregressive polynomial is at the unit root; the filter has no starting level.",
         call. = FALSE)
  }
  f0 <- omega / (1 - sb)
  df0 <- numeric(np)
  df0[1L] <- 1 / (1 - sb)
  if (q > 0L) {
    df0 <- df0 + (omega / (1 - sb)^2) *
      colSums(db[seq_len(q), , drop = FALSE])
  }
  list(omega = omega, a = a, b = b, db = db, f0 = f0, df0 = df0,
       i_a = i_a, np = np)
}

#' @title Filter a Score-Driven Term
#' @name term_filter.GasTerm
#' @description
#' Runs the score-driven recursion over each group in time order and
#' returns the predictor with its dynamic level added, together with the
#' exact derivative of that predictor with respect to the term's
#' parameters, propagated alongside the state.
#' @param term A built `GasTerm`.
#' @param eta The static part of the predictor.
#' @param y The response, unused directly: it reaches the filter through
#'   `score` and `curvature`.
#' @param score A function of the predictor returning
#'   \eqn{\partial\ell/\partial\eta} per observation.
#' @param curvature A function of the predictor returning
#'   \eqn{\partial^2\ell/\partial\eta^2} per observation.
#' @param psi The parameters, named as [term_params()].
#' @param ... Unused.
#' @param fast The fast context of the caller, or `NULL`: a list with
#'   `family` (the distribution's S7 class name), `link` (the
#'   parameter's link name), `k` (the parameter's 1-based index),
#'   `bounds`, `y` and `theta` (the per-observation
#'   parameters). Where the C registries of \pkg{distributions7} and
#'   \pkg{linkfunctions7} cover the pair, the recursion reads the score and
#'   the curvature through their scalar entry points instead of the R
#'   callbacks, bit-identically; where they do not, the context is inert
#'   and the callbacks run as before.
#' @param threads How many threads the recursion may use, over groups and
#'   only on the fast route: a group's filter is independent of the others
#'   and its writes land on its own rows, so no reduction is split and the
#'   result does not depend on the count, bit for bit.
#' @return A list with `eta`, `jacobian` and `curv`, the
#'   curvature read at each predictor.
#' @keywords internal
S7::method(term_filter, GasTerm) <- function(term, eta, y, score, curvature,
                                             psi, ..., fast = NULL,
                                             threads = 1L) {
  bp <- term@blueprint
  if (!length(bp)) {
    stop("the term has not been built; call term_build(term, data) first.",
         call. = FALSE)
  }
  nm <- term_params(term)
  psi <- unlist(psi[nm])
  if (length(psi) != length(nm) || anyNA(psi)) {
    stop(sprintf("'psi' must supply %s.", paste(nm, collapse = ", ")),
         call. = FALSE)
  }
  if (!is.null(bp$sub)) {
    return(.gas_filter_sub(term, eta, y, score, curvature, psi,
                           fast = fast, threads = threads))
  }
  p <- term@p
  q <- term@q
  cf <- .gas_coefs(psi[.gas_base_params(p, q)], p, q)
  out <- gas_filter_cpp(eta, bp$order, p, q, cf$omega, cf$a, cf$b, cf$db,
                        cf$f0, cf$df0, cf$i_a, cf$np, score, curvature,
                        fast = fast, threads = as.integer(threads))
  jac <- out$jacobian
  colnames(jac) <- nm
  list(eta = out$eta, jacobian = jac, curv = out$curv)
}

#' @title Filter a Score-Driven Term Backwards
#' @name term_adjoint.GasTerm
#' @description
#' Runs the recursion of [term_filter()] in reverse, returning the
#' derivative of a caller's objective with respect to the static predictor
#' it supplied and with respect to the sequence of scores it returned.
#' @param term A built `GasTerm`.
#' @param eta The static part of the predictor.
#' @param y The response, unused directly.
#' @param score,curvature The callbacks of [term_filter()].
#' @param psi The parameters, named as [term_params()].
#' @param g The direct derivative of the objective in the predictor the
#'   filter produced, one value per observation.
#' @param ... Unused.
#' @param fast,threads The fast context and the thread count of
#'   [term_filter()], passed to the forward pass the adjoint
#'   re-runs. The reverse pass reads the curvature sequence that pass
#'   returns, so with a covered context the adjoint evaluates no R
#'   callback at all; without one the callbacks run in the forward pass
#'   only, once per observation instead of twice.
#' @return A list with `deta` and `dscore`.
#' @keywords internal
S7::method(term_adjoint, GasTerm) <- function(term, eta, y, score, curvature,
                                              psi, g, ..., fast = NULL,
                                              threads = 1L) {
  bp <- term@blueprint
  if (!length(bp)) {
    stop("the term has not been built; call term_build(term, data) first.",
         call. = FALSE)
  }
  nm <- term_params(term)
  v <- unlist(psi[nm])
  if (length(v) != length(nm) || anyNA(v)) {
    stop(sprintf("'psi' must supply %s.", paste(nm, collapse = ", ")),
         call. = FALSE)
  }
  g <- as.numeric(g)
  if (length(g) != bp$n) {
    stop(sprintf("'g' must give one value per observation (%d).", bp$n),
         call. = FALSE)
  }
  if (!is.null(bp$sub)) {
    return(.gas_adjoint_sub(term, eta, y, score, curvature, v, g,
                            fast = fast, threads = threads))
  }
  # the forward pass already read the curvature at every predictor it
  # produced, so the reverse pass looks the sequence up rather than
  # evaluating the callback again at the same points
  cvv <- term_filter(term, eta, y, score, curvature, psi,
                     fast = fast, threads = threads)$curv

  p <- term@p
  q <- term@q
  base <- .gas_base_params(p, q)

  cf <- .gas_coefs(v[base], p, q)
  deta <- numeric(bp$n)
  dscore <- numeric(bp$n)
  for (l in seq_along(bp$order)) {
    rows <- bp$order[[l]]
    a <- cf$a
    b <- cf$b
    k <- length(rows)
    fb <- numeric(k)
    sb <- numeric(k)
    for (t in rev(seq_len(k))) {
      row <- rows[t]
      # e_t reaches the objective directly and through every later score
      eb <- g[row] + sb[t] * cvv[row]
      fb[t] <- fb[t] + eb
      deta[row] <- eb
      dscore[row] <- sb[t]
      if (q > 0L) {
        for (j in seq_len(q)) {
          if (t - j >= 1L) fb[t - j] <- fb[t - j] + b[[j]] * fb[t]
        }
      }
      if (p > 0L) {
        for (i in seq_len(p)) {
          if (t - i >= 1L) sb[t - i] <- sb[t - i] + a[[i]] * fb[t]
        }
      }
    }
  }
  list(deta = deta, dscore = dscore)
}

#' The Chart's Derivatives in the Unconstrained Parameters
#'
#' @description
#' The level, the score loadings and the autoregressive coefficients as
#' functions of the term's unconstrained parameters, with their first and
#' second derivatives in those.
#'
#' @details
#' The level and the loadings each reach the recursion through their own
#' link, so their first derivative is the link's and their second, on the
#' diagonal, is [linkfunctions7::d2linkinv()]; on the identity
#' both collapse to one and zero. The persistence reaches the coefficients
#' through two maps, the link onto the partial autocorrelations and
#' Levinson-Durbin onto the coefficients, so its second derivative
#' carries both a term in the map's own curvature and one in the link's.
#'
#' @param zeta The term's parameters on the unconstrained scale.
#' @param p,q The score and autoregressive orders.
#' @param links The links, as [term_links()] gives them.
#'
#' @return A list with the values and the derivative arrays.
#'
#' @keywords internal
.gas_chart_derivs <- function(zeta, p, q, links) {
  np <- length(zeta)
  nm <- names(zeta)
  i_om <- 1L
  i_a <- if (p > 0L) 1L + seq_len(p) else integer(0)
  i_pa <- if (q > 0L) 1L + p + seq_len(q) else integer(0)

  lk_om <- links[[nm[i_om]]]
  omega <- linkfunctions7::linkinv(lk_om, zeta[[i_om]])
  a <- if (p > 0L) vapply(seq_len(p), function(i)
    linkfunctions7::linkinv(links[[nm[i_a[i]]]], zeta[[i_a[i]]]),
    numeric(1)) else numeric(0)
  d_omega <- numeric(np)
  d_omega[i_om] <- linkfunctions7::dlinkinv(lk_om, zeta[[i_om]])
  h_omega <- matrix(0, np, np)
  h_omega[i_om, i_om] <- linkfunctions7::d2linkinv(lk_om, zeta[[i_om]])
  d_a <- lapply(seq_len(max(p, 1L)), function(i) {
    v <- numeric(np)
    if (p > 0L) {
      v[i_a[i]] <- linkfunctions7::dlinkinv(links[[nm[i_a[i]]]],
                                            zeta[[i_a[i]]])
    }
    v
  })
  h_a <- lapply(seq_len(max(p, 1L)), function(i) {
    h <- matrix(0, np, np)
    if (p > 0L) {
      h[i_a[i], i_a[i]] <- linkfunctions7::d2linkinv(links[[nm[i_a[i]]]],
                                                     zeta[[i_a[i]]])
    }
    h
  })

  b <- numeric(0)
  d_b <- list()
  h_b <- list()
  if (q > 0L) {
    lk <- links[[nm[i_pa[1L]]]]
    z <- zeta[i_pa]
    rho <- linkfunctions7::linkinv(lk, z)
    k1 <- linkfunctions7::dlinkinv(lk, z)
    k2 <- linkfunctions7::d2linkinv(lk, z)
    ld <- gas_levinson2(rho)
    b <- ld$phi
    for (j in seq_len(q)) {
      g1 <- numeric(np)
      g1[i_pa] <- ld$jacobian[j, ] * k1
      d_b[[j]] <- g1
      h <- matrix(0, np, np)
      # the map's curvature scaled by the chart, plus the chart's own on
      # the diagonal, which is where a non-identity link contributes. The
      # sub-block is built whole and then placed: at q = 1 an index pair of
      # length one collapses to a scalar, and `diag<-` on a scalar is an
      # error rather than a one-by-one assignment
      sub <- ld$hessian[[j]] * outer(k1, k1) +
        diag(ld$jacobian[j, ] * k2, nrow = q)
      h[i_pa, i_pa] <- sub
      h_b[[j]] <- h
    }
  }
  list(omega = omega, a = a, b = b, d_omega = d_omega, h_omega = h_omega,
       d_a = d_a, h_a = h_a, d_b = d_b, h_b = h_b, np = np)
}

#' The Chart's Third Derivatives, in One Direction
#'
#' @description
#' [.gas_chart_derivs()] differentiated once more and contracted
#' against a single direction in the term's own coordinates: one matrix for
#' the level, one per score loading and one per autoregressive coefficient.
#'
#' @details
#' The level and the loadings ride scalar links, so each of their third
#' derivatives is a single diagonal entry, the link's own \eqn{h'''} times
#' the direction's component there. The persistence is a composition, the
#' Levinson-Durbin map read at \eqn{\rho = h^{-1}(z)}, and differentiating
#' \eqn{B(z) = \phi(\rho(z))} three times and contracting the last slot gives
#' \deqn{T_{kl}h'_kh'_l + P_{kl}(h''_kv_kh'_l + h'_kh''_lv_l)
#'   + \delta_{kl}\big((Hw)_k h''_k + P_kh'''_kv_k\big),}
#' with \eqn{w_m = h'_mv_m} the direction pushed onto the partial
#' autocorrelations, \eqn{T} the contracted third derivative of the map and
#' \eqn{P}, \eqn{H} its first two.
#'
#' @param zeta The term's base parameters on the unconstrained scale.
#' @param p,q The score and autoregressive orders.
#' @param links The links, as [term_links()] gives them.
#' @param vz The direction, in the same coordinates as `zeta`.
#'
#' @return A list with `t_omega`, `t_a` and `t_b`, each a
#'   matrix or a list of matrices over the term's base coordinates.
#'
#' @seealso [.gas_chart_derivs()], [gas_levinson3()]
#'
#' @keywords internal
.gas_chart_derivs3 <- function(zeta, p, q, links, vz) {
  np <- length(zeta)
  nm <- names(zeta)
  i_om <- 1L
  i_a <- if (p > 0L) 1L + seq_len(p) else integer(0)
  i_pa <- if (q > 0L) 1L + p + seq_len(q) else integer(0)

  t_omega <- matrix(0, np, np)
  t_omega[i_om, i_om] <- linkfunctions7::d3linkinv(links[[nm[i_om]]],
                                                   zeta[[i_om]]) * vz[[i_om]]
  t_a <- lapply(seq_len(max(p, 1L)), function(i) {
    h <- matrix(0, np, np)
    if (p > 0L) {
      h[i_a[i], i_a[i]] <- linkfunctions7::d3linkinv(links[[nm[i_a[i]]]],
                                                     zeta[[i_a[i]]]) *
        vz[[i_a[i]]]
    }
    h
  })

  t_b <- list()
  if (q > 0L) {
    lk <- links[[nm[i_pa[1L]]]]
    z <- zeta[i_pa]
    v <- vz[i_pa]
    rho <- linkfunctions7::linkinv(lk, z)
    k1 <- linkfunctions7::dlinkinv(lk, z)
    k2 <- linkfunctions7::d2linkinv(lk, z)
    k3 <- linkfunctions7::d3linkinv(lk, z)
    w <- k1 * v
    ld <- gas_levinson2(rho)
    t3 <- gas_levinson3(rho, w)
    for (j in seq_len(q)) {
      Hj <- ld$hessian[[j]]
      Hw <- as.numeric(Hj %*% w)
      # the sub-block is built whole and then placed: at q = 1 an index pair
      # of length one collapses to a scalar, which is the trap the second
      # order already records
      sub <- t3[[j]] * outer(k1, k1) +
        Hj * outer(k2 * v, k1) + Hj * outer(k1, k2 * v) +
        diag(Hw * k2 + ld$jacobian[j, ] * k3 * v, nrow = q)
      h <- matrix(0, np, np)
      h[i_pa, i_pa] <- sub
      t_b[[j]] <- h
    }
  }
  list(t_omega = t_omega, t_a = t_a, t_b = t_b)
}

#' @title Second Derivatives of a Score-Driven Predictor
#' @name term_curvature.GasTerm
#' @description
#' The forward Jacobian of the filter's predictor in a caller's unknowns and
#' the second derivative contracted against the caller's weights, both
#' propagated through the recursion beside the state.
#' @details
#' With subformulas the general per-observation route runs instead, and
#' the second derivative is accumulated on each group's active set rather
#' than as a square over all the unknowns: a development's coordinate
#' reaches only the groups where its column is not identically zero, so
#' with grouping indicators the active set has the same size whether the
#' panel has ten groups or a thousand. Measured, the full square cost
#' 0.39 s at 124 unknowns over 1600 rows and would have reached about
#' twelve minutes at five hundred groups.
#'
#' `blocks` is called with the row of the jacobian restricted to the
#' active set and with that set, and returns its pieces in the same
#' coordinates. A callback of the earlier three-argument shape is still
#' accepted and given the full row, its result being subset here; it costs
#' the quadratic allocation the restriction exists to avoid.
#' @param term A built `GasTerm`.
#' @param eta The static part of the predictor.
#' @param y The response, unused directly.
#' @param score,curvature The callbacks of [term_filter()].
#' @param psi The parameters on the parameter scale, named as
#'   [term_params()].
#' @param g The weights the second derivative is contracted against.
#' @param seed The derivative of the static predictor in the unknowns.
#' @param blocks The model's own derivative pieces; see
#'   [term_curvature()].
#' @param ... Unused.
#' @param score_values,curvature_values The score and curvature of the
#'   model's log-density evaluated at the current predictors, one value per
#'   observation, on the parameter scale the callbacks read. Supplying both,
#'   together with `blocks_data`, routes the second-order recursion of
#'   the subformula route through the compiled kernel; either `NULL`
#'   keeps the R route.
#' @param blocks_data The model's derivative pieces as data instead of as a
#'   callback: a list with `H` (the mixed second derivatives, one
#'   column per distribution parameter), `D3` (the third derivatives,
#'   one column per parameter pair, pair `(r, r2)` at column
#'   `(r - 1) * np + r2`), `Vs` (the per-parameter jacobian rows
#'   of the other equations) and `ap` (the filter's own parameter
#'   index). Read only by the compiled route.
#' @param threads Threads for the compiled route's group loop, as
#'   `numericals7::n_threads()` counts them; 1 is sequential.
#' @return A list with `jacobian` and `curvature`.
#' @keywords internal
S7::method(term_curvature, GasTerm) <- function(term, eta, y, score,
                                                curvature, psi, g, seed,
                                                blocks, ...,
                                                score_values = NULL,
                                                curvature_values = NULL,
                                                blocks_data = NULL,
                                                threads = 1L) {
  .gas_curvature_core(term, eta, y, score, curvature, psi, g, seed, blocks,
                      direction = NULL, score_values = score_values,
                      curvature_values = curvature_values,
                      blocks_data = blocks_data, threads = threads)
}

#' @title Third Derivatives of a Score-Driven Predictor
#' @name term_third.GasTerm
#' @description
#' The second derivative of the filter's predictor differentiated once more
#' along one direction, propagated beside the state exactly as the first two
#' orders are.
#' @details
#' The recursion gains one state, \eqn{\Psi_t = \partial^3f_t/\partial u^3[v]},
#' seeded by the third derivative of the score in the same way \eqn{\Phi} is
#' seeded by its second. Everything else it needs, meaning the directional
#' derivatives of \eqn{F}, \eqn{\Phi}, \eqn{\dot S} and \eqn{\ddot S}, is a
#' contraction of a quantity the second-order recursion already carries, so
#' no second recursion is run and no three-index array is formed.
#'
#' The chart contributes its own third derivatives: the level and the
#' loadings through their scalar links, the persistence through
#' [.gas_chart_derivs3()], which composes
#' [gas_levinson3()] with them.
#' @param term A built `GasTerm`.
#' @param eta The static part of the predictor.
#' @param y The response, unused directly.
#' @param score,curvature The callbacks of [term_filter()].
#' @param psi The parameters on the parameter scale.
#' @param g The weights the third derivative is contracted against.
#' @param seed The derivative of the static predictor in the unknowns.
#' @param blocks The model's derivative pieces; see [term_third()].
#' @param direction The direction to contract against.
#' @param ... Unused.
#' @return A list with `jacobian`, `dphi` and `curvature`.
#' @keywords internal
S7::method(term_third, GasTerm) <- function(term, eta, y, score, curvature,
                                            psi, g, seed, blocks, direction,
                                            ...) {
  .gas_curvature_core(term, eta, y, score, curvature, psi, g, seed, blocks,
                      direction = direction)
}

#' The Score-Driven Recursion's Second and Third Derivatives
#'
#' @description
#' The body [term_curvature()] and [term_third()] share.
#' With `direction` `NULL` it propagates the first two derivatives
#' of the predictor; with a direction it propagates the third as well,
#' contracted against it.
#'
#' @details
#' The two orders are written here once, in place of a method each. The
#' third order's recursion reads \eqn{F}, \eqn{\Phi}, \eqn{\dot S} and
#' \eqn{\ddot S} at every lag, so a separate implementation would carry a
#' second copy of the first two orders, and the two would drift.
#'
#' @param term A built `GasTerm`.
#' @param eta The static part of the predictor.
#' @param y The response.
#' @param score,curvature The callbacks of [term_filter()].
#' @param psi The parameters on the parameter scale.
#' @param g The weights the contraction is taken against.
#' @param seed The derivative of the static predictor in the unknowns.
#' @param blocks The model's derivative pieces.
#' @param direction The direction, or `NULL` for the second order alone.
#'
#' @return A list with `jacobian` and `curvature`, and `dphi`
#'   where a direction was given.
#'
#' @keywords internal
.gas_curvature_core <- function(term, eta, y, score, curvature, psi, g, seed,
                                blocks, direction = NULL,
                                score_values = NULL,
                                curvature_values = NULL,
                                blocks_data = NULL, threads = 1L) {
  bp <- term@blueprint
  if (!length(bp)) {
    stop("the term has not been built; call term_build(term, data) first.",
         call. = FALSE)
  }
  nm <- term_params(term)
  links <- term_links(term)
  psiv <- unlist(psi[nm])
  third <- !is.null(direction)
  if (!is.null(bp$sub)) {
    return(.gas_curvature_sub(term, eta, y, score, curvature, psiv, g,
                              seed, blocks, direction,
                              score_values = score_values,
                              curvature_values = curvature_values,
                              blocks_data = blocks_data, threads = threads))
  }
  base <- .gas_base_params(term@p, term@q)
  nb <- length(base)
  ng <- length(bp$order)
  # the recursion is driven by the parameters, and the caller's unknowns
  # reach them through the links, so the chart is differentiated on the
  # UNCONSTRAINED scale, which is what the model estimates
  zeta <- vapply(base, function(j)
    linkfunctions7::linkfun(links[[j]], psiv[[j]]), numeric(1))

  seed <- as.matrix(seed)
  m <- ncol(seed)
  np <- length(nm)
  if (nrow(seed) != bp$n) {
    stop(sprintf("'seed' must have one row per observation (%d).", bp$n),
         call. = FALSE)
  }
  g <- as.numeric(g)
  if (length(g) != bp$n) {
    stop(sprintf("'g' must give one value per observation (%d).", bp$n),
         call. = FALSE)
  }
  # the term's own parameters are the LAST np columns of the unknowns, by
  # the convention the caller and this method share
  if (m < np) {
    stop(sprintf("'seed' has %d columns and the term has %d parameters.",
                 m, np), call. = FALSE)
  }
  zcol <- m - np + seq_len(np)
  p <- term@p
  q <- term@q
  if (third) {
    direction <- as.numeric(direction)
    if (length(direction) != m) {
      stop(sprintf("'direction' must have one value per unknown (%d).", m),
           call. = FALSE)
    }
  }

  # every group shares the scalar parameters, so the chart, the lifted
  # derivative arrays and the starting level are built once
  active_of <- function(l) seq_len(m)

  lift_of <- function(l, act) {
    A <- matrix(0, nb, length(act))
    A[cbind(seq_len(nb), match(zcol[seq_len(nb)], act))] <- 1
    A
  }

  prep <- function(l, act) {
    mk <- length(act)
    ch <- .gas_chart_derivs(zeta, p, q, links)
    A <- lift_of(l, act)
    lift <- function(v) as.numeric(crossprod(A, v))
    lift2 <- function(h) crossprod(A, h %*% A)
    om_u <- lift(ch$d_omega)
    om_uu <- lift2(ch$h_omega)
    a_u <- lapply(ch$d_a, lift)
    a_uu <- lapply(ch$h_a, lift2)
    b_u <- lapply(ch$d_b, lift)
    b_uu <- lapply(ch$h_b, lift2)
    sb <- if (q > 0L) sum(ch$b) else 0
    if (abs(1 - sb) < 1e-10) {
      stop("the autoregressive polynomial is at the unit root; the filter has no starting level.",
           call. = FALSE)
    }
    # the starting level and its two derivatives: f0 = omega/(1 - sum b),
    # with omega's own chart curvature entering the first term
    db_sum <- if (q > 0L) Reduce(`+`, b_u) else numeric(mk)
    f0 <- ch$omega / (1 - sb)
    f0_u <- om_u / (1 - sb) + ch$omega * db_sum / (1 - sb)^2
    f0_uu <- om_uu / (1 - sb) +
      (outer(om_u, db_sum) + outer(db_sum, om_u)) / (1 - sb)^2 +
      2 * ch$omega * outer(db_sum, db_sum) / (1 - sb)^3
    if (q > 0L) {
      f0_uu <- f0_uu + ch$omega * Reduce(`+`, b_uu) / (1 - sb)^2
    }
    out <- list(omega = ch$omega, a = ch$a, b = ch$b, om_u = om_u,
                om_uu = om_uu, a_u = a_u, a_uu = a_uu, b_u = b_u,
                b_uu = b_uu, f0 = f0, f0_u = f0_u, f0_uu = f0_uu)
    if (!third) return(out)

    # the third order, contracted against the direction. Every quantity here
    # is either a contraction of one the second order already carries or the
    # chart's own third derivative; nothing of three indices is formed
    vk <- direction[act]
    v_base <- direction[zcol[seq_len(nb)]]
    ch3 <- .gas_chart_derivs3(zeta, p, q, links, v_base)
    out$t_om <- lift2(ch3$t_omega)
    out$t_a <- lapply(ch3$t_a, lift2)
    out$t_b <- lapply(ch3$t_b, lift2)
    out$dom <- sum(om_u * vk)
    out$dom2 <- as.numeric(om_uu %*% vk)
    out$da <- lapply(a_u, function(x) sum(x * vk))
    out$da2 <- lapply(a_uu, function(x) as.numeric(x %*% vk))
    out$db <- lapply(b_u, function(x) sum(x * vk))
    out$db2 <- lapply(b_uu, function(x) as.numeric(x %*% vk))
    out$vk <- vk

    # the starting level's third derivative. f0 = omega/(1 - S) with
    # S the sum of the autoregressive coefficients, so writing c = 1/(1-S)
    # the second order is c*om_uu + c^2(om_u S_u + S_u om_u)
    # + 2 omega c^3 S_u S_u + omega c^2 S_uu, and this differentiates it
    cst <- 1 / (1 - sb)
    S_u <- db_sum
    S_uu <- if (q > 0L) Reduce(`+`, b_uu) else matrix(0, mk, mk)
    S_3 <- if (q > 0L) Reduce(`+`, out$t_b) else matrix(0, mk, mk)
    dS <- sum(S_u * vk)
    dS2 <- as.numeric(S_uu %*% vk)
    om <- ch$omega
    out$f0_3 <- out$t_om * cst + om_uu * (cst^2 * dS) +
      2 * cst^3 * dS * (outer(om_u, S_u) + outer(S_u, om_u)) +
      cst^2 * (outer(out$dom2, S_u) + outer(om_u, dS2) +
               outer(dS2, om_u) + outer(S_u, out$dom2)) +
      (2 * out$dom * cst^3 + 6 * om * cst^4 * dS) * outer(S_u, S_u) +
      2 * om * cst^3 * (outer(dS2, S_u) + outer(S_u, dS2)) +
      out$dom * cst^2 * S_uu + 2 * om * cst^3 * dS * S_uu + om * cst^2 * S_3
    out$df0 <- sum(f0_u * vk)
    out$dphi0 <- as.numeric(f0_uu %*% vk)
    out
  }
  shared <- prep(1L, seq_len(m))

  # does the caller take the active set? A three-argument callback is the
  # earlier contract and is given the full row instead
  wants_act <- length(formals(blocks)) >= 4L

  D <- matrix(0, bp$n, m)
  W <- matrix(0, m, m)
  dP <- if (third) matrix(0, bp$n, m) else NULL
  for (l in seq_len(ng)) {
    rows <- bp$order[[l]]
    act <- active_of(l)
    mk <- length(act)
    gp <- if (is.null(shared)) prep(l, act) else shared
    Wl <- matrix(0, mk, mk)
    a <- gp$a
    b <- gp$b
    om_u <- gp$om_u
    om_uu <- gp$om_uu
    a_u <- gp$a_u
    a_uu <- gp$a_uu
    b_u <- gp$b_u
    b_uu <- gp$b_uu
    f0 <- gp$f0
    f0_u <- gp$f0_u
    f0_uu <- gp$f0_uu
    vk <- gp$vk
    k <- length(rows)
    f <- numeric(k)
    s <- numeric(k)
    F_ <- vector("list", k)
    Phi <- vector("list", k)
    Sd <- vector("list", k)
    Sdd <- vector("list", k)
    Psi <- if (third) vector("list", k) else NULL
    Sddd <- if (third) vector("list", k) else NULL
    for (t in seq_len(k)) {
      row <- rows[t]
      ft <- gp$omega
      Ft <- om_u
      # the level's own chart curvature, zero on the identity
      Pt <- om_uu
      Tt <- if (third) gp$t_om else NULL
      if (p > 0L) {
        for (i in seq_len(p)) {
          lag <- t - i
          s_l <- if (lag >= 1L) s[lag] else 0
          Sd_l <- if (lag >= 1L) Sd[[lag]] else numeric(mk)
          Sdd_l <- if (lag >= 1L) Sdd[[lag]] else matrix(0, mk, mk)
          ft <- ft + a[[i]] * s_l
          Ft <- Ft + a[[i]] * Sd_l + s_l * a_u[[i]]
          Pt <- Pt + a[[i]] * Sdd_l +
            outer(a_u[[i]], Sd_l) + outer(Sd_l, a_u[[i]]) +
            s_l * a_uu[[i]]
          if (third) {
            Sddd_l <- if (lag >= 1L) Sddd[[lag]] else matrix(0, mk, mk)
            dSd_l <- if (lag >= 1L) sum(Sd_l * vk) else 0
            dSdd_l <- if (lag >= 1L) as.numeric(Sdd_l %*% vk) else numeric(mk)
            Tt <- Tt + gp$da[[i]] * Sdd_l + a[[i]] * Sddd_l +
              outer(gp$da2[[i]], Sd_l) + outer(Sd_l, gp$da2[[i]]) +
              outer(a_u[[i]], dSdd_l) + outer(dSdd_l, a_u[[i]]) +
              dSd_l * a_uu[[i]] + s_l * gp$t_a[[i]]
          }
        }
      }
      if (q > 0L) {
        for (j in seq_len(q)) {
          lag <- t - j
          f_l <- if (lag >= 1L) f[lag] else f0
          F_l <- if (lag >= 1L) F_[[lag]] else f0_u
          Phi_l <- if (lag >= 1L) Phi[[lag]] else f0_uu
          ft <- ft + b[[j]] * f_l
          Ft <- Ft + b[[j]] * F_l + f_l * b_u[[j]]
          Pt <- Pt + b[[j]] * Phi_l +
            outer(b_u[[j]], F_l) + outer(F_l, b_u[[j]]) + f_l * b_uu[[j]]
          if (third) {
            Psi_l <- if (lag >= 1L) Psi[[lag]] else gp$f0_3
            dF_l <- if (lag >= 1L) sum(F_l * vk) else gp$df0
            dPhi_l <- if (lag >= 1L) as.numeric(Phi_l %*% vk) else gp$dphi0
            Tt <- Tt + gp$db[[j]] * Phi_l + b[[j]] * Psi_l +
              outer(gp$db2[[j]], F_l) + outer(F_l, gp$db2[[j]]) +
              outer(b_u[[j]], dPhi_l) + outer(dPhi_l, b_u[[j]]) +
              dF_l * b_uu[[j]] + f_l * gp$t_b[[j]]
          }
        }
      }
      f[t] <- ft
      F_[[t]] <- Ft
      Phi[[t]] <- Pt
      if (third) Psi[[t]] <- Tt

      e_t <- eta[row] + ft
      # the row of the jacobian on the active set; every other column of it
      # is zero, the seed carrying no other group's parameters
      Dt <- seed[row, act] + Ft
      D[row, act] <- Dt
      Wl <- Wl + g[row] * (if (third) Tt else Pt)

      s[t] <- score(e_t, row)
      cv <- curvature(e_t, row)
      bl <- if (wants_act) blocks(e_t, row, Dt, act) else {
        full <- numeric(m)
        full[act] <- Dt
        b3 <- blocks(e_t, row, full)
        list(cross = b3$cross[act], M = b3$M[act, act, drop = FALSE],
             dcurv = b3$dcurv[act],
             N = if (third) b3$N[act, act, drop = FALSE] else NULL)
      }
      Sd[[t]] <- cv * Dt + bl$cross
      Sdd[[t]] <- cv * Pt + bl$M
      if (third) {
        # the second derivative of the predictor contracted along the
        # direction is the derivative of the jacobian ROW, which is what
        # differentiating M's own V_p asks for; the derivative of the
        # PREDICTOR along it is the scalar Dt . v, which the caller forms
        dPhi_t <- as.numeric(Pt %*% vk)
        dP[row, act] <- dPhi_t
        Sddd[[t]] <- sum(bl$dcurv * vk) * Pt + cv * Tt + bl$N +
          outer(dPhi_t, bl$dcurv) + outer(bl$dcurv, dPhi_t)
      }
    }
    W[act, act] <- W[act, act] + Wl
  }
  # The matrix is symmetric and the accumulation is not: an entry and its
  # transpose collect the same terms in a different ORDER, and (x + p) + q
  # is not (x + q) + p in floating point. The gap is of the order of the
  # rounding, and a caller about to factor this wants an exactly symmetric
  # matrix rather than one that is nearly so.
  W <- (W + t(W)) / 2
  if (third) return(list(jacobian = D, dphi = dP, curvature = W))
  list(jacobian = D, curvature = W)
}

#' @title Print a Score-Driven Term
#' @name print.GasTerm
#'
#' @description
#' Prints the label, the two orders and, for a built term, over how many groups
#' the filter runs. A second line names any developed parameters, and a third
#' lists the parameters themselves.
#'
#' @details
#' The form is
#'
#' ```
#' <GasTerm> 'gas': score-driven, p = 1, q = 1; 3 group(s)
#'   developed: omega
#'   parameters: omega.(Intercept), omega.g2, alpha1, pacf1
#' ```
#'
#' The `developed` line appears only where a subformula was given. The
#' parameter list is [term_params()], so a developed parameter shows there as
#' its coefficients.
#'
#' A built term is described by its group count rather than by the word
#' "built", the count being the more useful line.
#'
#' @param x A [GasTerm()], built or not.
#' @param ... Unused, and accepted so that the signature matches [print()]'s.
#'
#' @return `x`, invisibly. Called for the lines it writes.
#'
#' @seealso [gas()], [term_params()].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:60, y = rnorm(60), id = rep(1:3, each = 20),
#'                  g = factor(rep(c("u", "v"), 30)))
#'
#' # A specification, and the same term built over three groups.
#' gas(p = 1, q = 2)
#' term_build(gas(p = 1, q = 1, by = id, time = t), dd)
#'
#' # A developed parameter is named on its own line.
#' term_build(gas(p = 1, q = 1, omega ~ g, by = id, time = t), dd)
#'
#' @keywords internal
S7::method(print, GasTerm) <- function(x, ...) {
  built <- length(x@blueprint) > 0L
  cat(sprintf("<GasTerm> '%s': score-driven, p = %d, q = %d%s\n",
              x@label, x@p, x@q,
              if (built) sprintf("; %d group(s)", length(x@blueprint$order))
              else " (specification)"))
  if (length(x@submodels)) {
    cat("  developed: ", paste(names(x@submodels), collapse = ", "),
        "\n", sep = "")
  }
  cat("  parameters: ", paste(term_params(x), collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' The Score-Driven Recursion in R
#'
#' @description
#' The loop `gas_filter_cpp()` replaces, kept so the compiled
#' route has something to be compared against that shares none of its code.
#'
#' @param eta The static predictor.
#' @param order A list of row indices, one entry per group, in time order.
#' @param p,q The score and autoregressive orders.
#' @param omega The level.
#' @param a,b The score loadings and the autoregressive coefficients.
#' @param db The derivative of the coefficients in the parameters.
#' @param f0,df0 The starting level and its derivative.
#' @param i_a The positions of the score loadings among the parameters.
#' @param np The number of parameters.
#' @param score,curvature The callbacks of [term_filter()].
#'
#' @return A list with `eta`, `jacobian` and `curv`.
#'
#' @keywords internal
gas_filter_r <- function(eta, order, p, q, omega, a, b, db, f0, df0,
                         i_a, np, score, curvature) {
  n <- length(eta)
  eta_out <- numeric(n)
  jac <- matrix(0, n, np)
  cv <- numeric(n)

  for (rows in order) {
    m <- length(rows)
    f <- numeric(m)
    df <- matrix(0, m, np)
    s <- numeric(m)
    ds <- matrix(0, m, np)

    for (t in seq_len(m)) {
      f[t] <- omega
      df[t, ] <- 0
      df[t, 1L] <- 1
      if (p > 0L) {
        for (i in seq_len(p)) {
          s_lag <- if (t - i >= 1L) s[t - i] else 0
          ds_lag <- if (t - i >= 1L) ds[t - i, ] else numeric(np)
          f[t] <- f[t] + a[[i]] * s_lag
          df[t, ] <- df[t, ] + a[[i]] * ds_lag
          df[t, i_a[i]] <- df[t, i_a[i]] + s_lag
        }
      }
      if (q > 0L) {
        for (j in seq_len(q)) {
          f_lag <- if (t - j >= 1L) f[t - j] else f0
          df_lag <- if (t - j >= 1L) df[t - j, ] else df0
          f[t] <- f[t] + b[[j]] * f_lag
          df[t, ] <- df[t, ] + b[[j]] * df_lag + db[j, ] * f_lag
        }
      }
      e_t <- eta[rows[t]] + f[t]
      s[t] <- score(e_t, rows[t])
      cv[rows[t]] <- curvature(e_t, rows[t])
      ds[t, ] <- cv[rows[t]] * df[t, ]
    }

    eta_out[rows] <- eta[rows] + f
    jac[rows, ] <- df
  }
  list(eta = eta_out, jacobian = jac, curv = cv)
}


#' @name term_static_deriv.GasTerm
#'
#' @title The Filter's Sensitivity to Its Own Equation's Coefficients
#'
#' @description
#' The forward recursion of [term_static_deriv()] for a
#' score-driven term.
#'
#' @details
#' The recursion is the filter's own with the seed replaced: the starting
#' level depends on the term's parameters, never on the static predictor,
#' so the propagated derivative starts at zero and everything after it comes
#' from the scores. The autoregressive and loading coefficients are read at
#' each row, which covers a term whose parameters carry submodels of their
#' own as well as one whose parameters are constant.
#'
#' It is written in R. It evaluates no callback, runs once at a variance
#' instead of at every iteration of a fit, and costs one
#' pass over the data per column.
#'
#' @param term A built score-driven term.
#' @param curv The curvature at each predictor.
#' @param X The directions to propagate.
#' @param psi The term's parameters, on the parameter scale.
#' @param ... Ignored.
#'
#' @return A matrix of `X`'s dimensions.
#'
#' @seealso [term_static_deriv()], [term_filter()]
#'
#' @keywords internal
S7::method(term_static_deriv, GasTerm) <- function(term, curv, X, psi, ...) {
  bp <- term@blueprint
  if (!length(bp)) {
    stop("the term has not been built; call term_build(term, data) first.",
         call. = FALSE)
  }
  X <- as.matrix(X)
  n <- nrow(X)
  k <- ncol(X)
  if (!k) return(X)
  p <- term@p
  q <- term@q
  nm <- term_params(term)
  psi <- unlist(psi[nm])
  if (is.null(bp$sub)) {
    cf <- .gas_coefs(psi[.gas_base_params(p, q)], p, q)
    A <- matrix(rep(cf$a, each = n), n, p)
    B <- matrix(rep(cf$b, each = n), n, q)
  } else {
    vals <- .gas_sub_values(term, psi, "parameter")
    bd <- .gas_sub_b(term, vals)
    A <- matrix(0, n, p)
    for (i in seq_len(p)) A[, i] <- vals[[paste0("alpha", i)]]$v
    B <- matrix(0, n, q)
    for (j in seq_len(q)) B[, j] <- bd$B[, j]
  }
  U <- matrix(0, n, k)
  for (rows in bp$order) {
    m <- length(rows)
    u <- matrix(0, m, k)
    dv <- matrix(0, m, k)
    for (t in seq_len(m)) {
      r <- rows[[t]]
      ut <- numeric(k)
      if (p > 0L) {
        for (i in seq_len(p)) {
          if (t - i >= 1L) ut <- ut + A[r, i] * dv[t - i, ]
        }
      }
      if (q > 0L) {
        for (j in seq_len(q)) {
          if (t - j >= 1L) ut <- ut + B[r, j] * u[t - j, ]
        }
      }
      u[t, ] <- ut
      dv[t, ] <- curv[[r]] * (X[r, ] + ut)
    }
    U[rows, ] <- u
  }
  X + U
}


#' @name term_continue.GasTerm
#'
#' @title Continuing a Score-Driven Filter Past the Series
#'
#' @description
#' The level the recursion reaches at rows that come after the observed
#' ones.
#'
#' @details
#' The recursion is the filter's own, started from the level and score the
#' observed part ended at, never from the term's own seed. What changes
#' past the data is the score: it has zero conditional mean by construction,
#' the model's own defining property, so at a row whose response is not
#' observed the driving term is its expectation and the continuation is the
#' deterministic recursion
#' \deqn{f_{n+h} = \omega + \sum_i \alpha_i s_{n+h-i} +
#'   \sum_j \beta_j f_{n+h-j},}
#' the loadings contributing only while \eqn{n+h-i} is still an observed
#' time. No simulation and no integration is involved.
#'
#' A new row is placed by its own time within its own group, and must come
#' after every observed time of that group: a row falling inside the series
#' is not a continuation but a re-reading of it, where the response is known
#' and the filter must be run instead of continued, so it is rejected with
#' the rows named. A group the fit never saw is rejected for the same
#' reason: there is no state to continue.
#'
#' @param term A built score-driven term.
#' @param psi The term's parameters, named as [term_params()].
#' @param f_past The level at each observed row.
#' @param s_past The score at each observed row.
#' @param newdata The rows to continue onto.
#' @param ... Ignored.
#'
#' @return A numeric vector of `nrow(newdata)` levels.
#'
#' @seealso [term_continue()], [term_filter()]
#'
#' @keywords internal
S7::method(term_continue, GasTerm) <- function(term, psi, f_past, s_past,
                                               newdata, ...) {
  bp <- term@blueprint
  if (!length(bp)) {
    stop("the term has not been built; call term_build(term, data) first.",
         call. = FALSE)
  }
  nn <- nrow(newdata)
  if (!nn) return(numeric(0))
  p <- term@p
  q <- term@q
  nm <- term_params(term)
  u <- unlist(psi[nm])

  gf <- if (is.null(term@by)) factor(rep(1L, nn), levels = "1") else
    factor(as.character(eval(term@by, newdata, baseenv())),
           levels = bp$levels)
  tv <- if (is.null(term@time)) rep(NA_real_, nn) else
    eval(term@time, newdata, baseenv())
  if (is.null(term@time)) {
    stop("continuing a score-driven term past the series needs its `time`: ",
         "without one a row's
  place in the series is its position in the ",
         "data frame, which says nothing about new rows.", call. = FALSE)
  }
  if (anyNA(gf)) {
    stop(sprintf(paste0("rows %s name a group the fit never saw, so there ",
                        "is no state to continue."),
                 paste(utils::head(which(is.na(gf)), 5L), collapse = ", ")),
         call. = FALSE)
  }
  if (anyNA(tv)) stop("`time` must not be missing.", call. = FALSE)

  # the parameters at each row: constant, or the development read at the
  # NEW rows through each sub-term's own blueprint
  if (is.null(bp$sub)) {
    cf <- .gas_coefs(u[.gas_base_params(p, q)], p, q)
    om_new <- rep(cf$omega, nn)
    A_new <- matrix(rep(cf$a, each = nn), nn, p)
    B_new <- matrix(rep(cf$b, each = nn), nn, q)
  } else {
    v <- .gas_sub_new(term, u, newdata)
    om_new <- v$om
    A_new <- v$A
    B_new <- v$B
  }

  out <- numeric(nn)
  gi <- as.integer(gf)
  for (lv in sort(unique(gi))) {
    rows_old <- bp$order[[as.character(lv)]]
    if (is.null(rows_old)) {
      stop("a group of the new data has no observed rows to continue from.",
           call. = FALSE)
    }
    new <- which(gi == lv)
    new <- new[order(tv[new])]
    last <- max(bp$times[rows_old])
    if (any(tv[new] <= last)) {
      bad <- new[tv[new] <= last]
      stop(sprintf(paste0("rows %s fall inside the observed series (their ",
                          "group ends at time %s).\n  A row whose response ",
                          "is known is read by running the filter, not by ",
                          "continuing it."),
                   paste(utils::head(bad, 5L), collapse = ", "),
                   format(last)), call. = FALSE)
    }
    # the history the continuation starts from, oldest first
    f <- c(f_past[rows_old], numeric(length(new)))
    sc <- c(s_past[rows_old], numeric(length(new)))
    n0 <- length(rows_old)
    for (h in seq_along(new)) {
      t <- n0 + h
      r <- new[[h]]
      ft <- om_new[[r]]
      if (p > 0L) {
        for (i in seq_len(p)) {
          if (t - i >= 1L) ft <- ft + A_new[r, i] * sc[[t - i]]
        }
      }
      if (q > 0L) {
        for (j in seq_len(q)) {
          if (t - j >= 1L) ft <- ft + B_new[r, j] * f[[t - j]]
        }
      }
      f[[t]] <- ft
      # the score past the data is at its conditional mean, which is zero
      sc[[t]] <- 0
      out[[r]] <- ft
    }
  }
  out
}


#' @name term_simulate.GasTerm
#'
#' @title Generating From a Score-Driven Filter
#'
#' @description
#' Runs the recursion forward drawing the response at each step.
#'
#' @details
#' No separate recursion is written. The filter's own recursion is the
#' generative one. The level at one time is a function of the scores before
#' it, and a score is a function of a response and a predictor, so
#' the only difference is where the response comes from.
#' [term_filter()] calls its `score` callback exactly once
#' per observation, in time order within each group, at the predictor the
#' recursion has just produced; a callback that draws the response there,
#' keeps it, and returns the score of what it drew turns the filter into a
#' generator.
#'
#' The curvature is read at the same point, so it reads the response the
#' score drew, drawing no second one.
#'
#' The fast route is not taken: it reads the response through a registered
#' C entry point, and here the response does not exist until the step that
#' needs it.
#'
#' @param term A built score-driven term.
#' @param psi The term's parameters, on the parameter scale.
#' @param eta The static part of the predictor.
#' @param draw A function `(e, i)` returning one response value.
#' @param score,curvature Functions of `(y, e, i)` giving the
#'   derivatives of the log-density in the predictor at observation `i`.
#'   They take the response as an argument, unlike the filter's own, because
#'   here it is being made.
#' @param ... Ignored.
#'
#' @return A list with `eta`, `y` and `latent`, the level.
#'
#' @seealso [term_simulate()], [term_filter()]
#'
#' @keywords internal
S7::method(term_simulate, GasTerm) <- function(term, psi, eta, draw,
                                               score = NULL,
                                               curvature = NULL, ...) {
  n <- length(eta)
  yv <- rep(NA_real_, n)
  if (is.null(score)) score <- function(y, e, i) y - e
  if (is.null(curvature)) curvature <- function(y, e, i) -1
  get <- function(e, i) {
    if (is.na(yv[[i]])) yv[[i]] <<- as.numeric(draw(e, i))
    yv[[i]]
  }
  out <- term_filter(term, eta, yv,
                     function(e, i) as.numeric(score(get(e, i), e, i)),
                     function(e, i) as.numeric(curvature(get(e, i), e, i)),
                     psi, fast = NULL, threads = 1L)
  list(eta = out$eta, y = yv, latent = as.numeric(out$eta) - as.numeric(eta))
}
