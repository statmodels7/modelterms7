#' @include term_classes.R generics.R structural.R
NULL

#' @title S7 Class for Score-Driven Dynamics
#' @name GasTerm
#'
#' @description
#' A subclass of \code{\link{structural_term}} for a generalized
#' autoregressive score component: a time-varying level driven by the
#' score of the observation density, added to the linear predictor.
#' Constructed by \code{\link{gas}}.
#'
#' @inheritParams model_term
#' @param p The number of score lags.
#' @param q The number of autoregressive lags.
#' @param by An optional grouping expression, filtered independently.
#' @param time An optional ordering expression.
#' @param links The links overriding the defaults, if any.
#' @param submodels One optional subformula per parameter.
#' @param blueprint The resolved ordering and grouping.
#'
#' @return An object of class \code{GasTerm}.
#'
#' @seealso \code{\link{gas}}
#' @examples
#' S7::S7_inherits(gas(), GasTerm)
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
    submodels = S7::class_list,
    blueprint = S7::class_list
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
#' block, and \code{\link{term_filter}} runs the recursion instead. That
#' is what makes it a \code{\link{structural_term}}.
#'
#' What drives the recursion is the score of whatever distribution the
#' model carries, so the same term is a GARCH-like volatility model when
#' it enters the scale of a Gaussian, a dynamic count model when it enters
#' the mean of a Poisson, and a robust location filter when it enters a
#' Student t: a heavy-tailed score is bounded in the observation, so an
#' outlier moves the level by a bounded amount rather than in proportion
#' to its size.
#'
#' \subsection{The parameters and their chart}{
#' The parameters are the level \eqn{\omega}, the score loadings
#' \eqn{a_1, \dots, a_p}, and the persistence. Each is estimated on the
#' unconstrained scale of a link, and \code{links} overrides any of them;
#' the defaults are the following.
#'
#' The level carries the identity, being unconstrained. The loadings carry
#' the \strong{log} link: a positive loading responds in the direction of
#' the score, which is the case the score-driven literature writes, and
#' positivity is then structural -- a deviation or a submodel moves the
#' loading on the log scale and no group or observation can take a
#' negative one. A loading that must be free in sign is asked for with
#' \code{links = list(alpha1 = linkfunctions7::identity_link())}.
#'
#' The persistence is carried by \strong{partial autocorrelations} rather
#' than by the coefficients \eqn{b_j}: the stationary region of an
#' autoregression is not a box, so no collection of scalar links covers
#' it, while the partial autocorrelations each range over \eqn{(-1, 1)}
#' independently and the Levinson-Durbin recursion carries them onto the
#' coefficients bijectively. At \eqn{q = 1} the two coincide. The
#' coordinate is named for the chart it lives on, \code{pacf1} and so on,
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
#' spectral radius of its companion matrix rather than a box.
#'
#' The score driving the recursion is used unscaled. The general
#' formulation carries a scaling matrix, usually an inverse information,
#' which the curvature this term already receives would supply.
#' }
#'
#' \subsection{Groups and time}{
#' \code{by} filters each group independently, which is what a panel of
#' short series needs, and \code{time} gives the order within a group.
#' Without \code{time} the rows are taken in the order they appear.
#' }
#'
#' \subsection{A parameter developed with covariates}{
#' A two-sided formula in \code{...} whose left side names a parameter
#' develops it as \eqn{\psi_{j,t} = g_j^{-1}(z_t^\top\gamma_j)}, the
#' design \eqn{Z} built from the right-hand side through
#' \code{\link{interpret_formula}}, so it takes any additive term of the
#' package:
#' \preformatted{gas(p = 1, q = 1, omega ~ ridge(~g), alpha1 ~ s(x),
#'     pacf1 ~ random(~1 | id), by = id)}
#' The development acts on the unconstrained scale of the parameter's own
#' link, which is what keeps every per-observation value in the
#' parameter's own set whatever the coefficients are: a loading on the
#' log link is positive at every observation, a persistence on the rhobit
#' chart is inside \eqn{(-1, 1)} at every observation, and at \eqn{q = 1}
#' that bounds the recursion's growth step by step. The coefficients
#' \eqn{\gamma_j} are the term's parameters, unconstrained and on the
#' identity link; the penalties the sub-terms carry are reported through
#' \code{\link{term_penalties}} under the key \code{parameter::subterm}.
#' A parameter that varies by observation changes the recursion itself,
#' \deqn{f_t = \omega_t + \sum_i a_{i,t}\, s_{t-i}
#'   + \sum_j b_{j,t}\, f_{t-j},}
#' with \eqn{b_t} from the Levinson-Durbin map of that observation's
#' partial autocorrelations, and the filter, its derivative, the reverse
#' pass and the curvature all run the general recursion.
#'
#' \code{by = ~f} (a formula, where a grouping variable is a bare symbol)
#' is the shorthand giving the same subformula to every parameter; mixing
#' it with per-parameter formulas is an error. A structural term, and a
#' term whose block moves with its own coefficients, are rejected inside
#' a subformula.
#' }
#'
#' \subsection{A population value and a departure per group}{
#' The panel case is one subformula:
#' \code{gas(omega ~ random(~1 | id), by = id)} is a population value (the
#' intercept of the development) plus one unconstrained departure per
#' group, shrunk by the random intercept's own ridge, whose hyperparameter
#' a fitting layer estimates. \code{lasso(...)} in the subformula sets the
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
#' Earlier releases spelled this case as \code{deviations =} with a
#' \code{penalty =}; both arguments are gone, the subformula reproducing
#' them exactly (the same fit to the printed digit, hyperparameter
#' included) while covering what they could not.
#' }
#'
#' @param p The number of score lags. Defaults to 1.
#' @param q The number of autoregressive lags. Defaults to 1.
#' @param ... Two-sided formulas whose left side names a parameter, one
#'   per parameter to be developed with covariates, e.g.
#'   \code{alpha1 ~ s(x)}; see the section above.
#' @param by An optional grouping variable, evaluated in the data; each
#'   group is filtered independently, from its own starting level. A
#'   FORMULA here is the shorthand giving the same subformula to every
#'   parameter, and then no grouping is implied.
#' @param time An optional ordering variable, evaluated in the data.
#' @param links An optional named list of \pkg{linkfunctions7} links over
#'   the parameters of \code{\link{term_params}}, overriding the defaults
#'   described above. A deviation cannot be named: it is unconstrained by
#'   construction, acting on the scale its parameter's own link defines.
#' @param label A single non-empty string naming the term.
#'
#' @return An object of class \code{\link{GasTerm}} (a specification; see
#'   \code{\link{term_build}}).
#'
#' @references
#' Creal, D., Koopman, S. J. and Lucas, A. (2013). Generalized
#' autoregressive score models with applications. \emph{Journal of Applied
#' Econometrics}, 28(5), 777--795.
#'
#' Harvey, A. C. (2013). \emph{Dynamic Models for Volatility and Heavy
#' Tails}. Cambridge University Press.
#'
#' @examples
#' term_params(gas(p = 1, q = 2))
#'
#' @seealso \code{\link{regime}}
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

S7::method(term_params, GasTerm) <- function(term, ...) {
  if (length(term@blueprint) && !is.null(term@blueprint$sub)) {
    return(.gas_sub_layout(term)$names)
  }
  .gas_base_params(term@p, term@q)
}

#' @title The Level of a Score-Driven Term
#' @name term_level_param.GasTerm
#' @description
#' \code{"omega"}, which adds a constant to the equation's predictor and is
#' therefore the direction an intercept there also spans. With the level
#' developed by a subformula, the coordinates whose design column is
#' constant, each of which shifts the predictor the same way.
#' @param term A \code{\link{GasTerm}}.
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
#' The design of \code{omega}'s development when it carries one, its
#' columns named after the coordinates, so a fitting layer can compare its
#' span with the equation's; \code{NULL} for a scalar level.
#' @param term A built \code{\link{GasTerm}}.
#' @param ... Unused.
#' @return A numeric matrix or \code{NULL}.
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
#' The level, the score loadings and the AUTOREGRESSIVE COEFFICIENTS of the
#' literature -- \code{omega}, \code{alpha1}, \code{beta1} -- with the
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
#' @param term A \code{\link{GasTerm}}.
#' @param zeta The parameters on the unconstrained scale.
#' @param ... Unused.
#' @return A list, as \code{\link{term_readable}} documents.
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
#' Zero is the natural point of every other chart -- a level of zero, no
#' persistence, no deviation -- so the term starts as near the model
#' without it as its charts allow. The loadings are the exception because
#' zero on the log scale is a loading of ONE, a response strong enough to
#' destabilize the recursion at ordinary curvatures; \eqn{0.1} is a weak
#' response, and it is applied on the parameter scale so the start means
#' the same thing whatever chart a loading rides.
#' @param term A \code{\link{GasTerm}}.
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
#' \code{parameter::subterm} and covering that sub-term's coefficients in
#' the term's own numbering. The scalar parameters are unpenalized, and
#' the list is empty for a specification, whose developments do not exist
#' until the term is built.
#' @param term A built \code{\link{GasTerm}}.
#' @param ... Unused.
#' @return A list of entries, as \code{\link{term_penalties}} documents.
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
        min_ratio = e$min_ratio)
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
#' @return A list with \code{phi}, the coefficients, and \code{jacobian},
#'   the matrix of their derivatives with respect to \code{pacf}.
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
#' \code{\link{gas_levinson}} with the second derivatives of the
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
#' @return A list with \code{phi}, \code{jacobian} and \code{hessian}, the
#'   last a list of one symmetric matrix per coefficient.
#'
#' @seealso \code{\link{gas_levinson}}
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
#' \code{\link{gas_levinson2}}'s second derivatives differentiated once more
#' and contracted against a single direction, one matrix per coefficient.
#'
#' @details
#' The exact gradient of a marginal criterion over a penalty on this term's
#' own parameters needs the third derivative of the predictor, and the
#' persistence reaches the predictor through this map. It is needed only
#' CONTRACTED: the criterion asks for \eqn{\mathrm{tr}(M\,\partial K/\partial
#' u[v])}, a derivative along the single direction the penalized mode moves
#' in, so what is propagated is a matrix per coefficient and never a
#' three-index array.
#'
#' Differentiating the hessian recursion of \code{\link{gas_levinson2}} once
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
#' @param w The direction to contract against, as long as \code{pacf}.
#'
#' @return A list of one \code{q} by \code{q} matrix per coefficient.
#'
#' @seealso \code{\link{gas_levinson2}}
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
  term@blueprint <- list(order = ord, n = n, levels = levels(gf),
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
#' @param term A built \code{GasTerm}.
#' @param eta The static part of the predictor.
#' @param y The response, unused directly: it reaches the filter through
#'   \code{score} and \code{curvature}.
#' @param score A function of the predictor returning
#'   \eqn{\partial\ell/\partial\eta} per observation.
#' @param curvature A function of the predictor returning
#'   \eqn{\partial^2\ell/\partial\eta^2} per observation.
#' @param psi The parameters, named as \code{\link{term_params}}.
#' @param ... Unused.
#' @return A list with \code{eta} and \code{jacobian}.
#' @keywords internal
S7::method(term_filter, GasTerm) <- function(term, eta, y, score, curvature,
                                             psi, ...) {
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
    return(.gas_filter_sub(term, eta, y, score, curvature, psi))
  }
  p <- term@p
  q <- term@q
  cf <- .gas_coefs(psi[.gas_base_params(p, q)], p, q)
  out <- gas_filter_cpp(eta, bp$order, p, q, cf$omega, cf$a, cf$b, cf$db,
                        cf$f0, cf$df0, cf$i_a, cf$np, score, curvature)
  jac <- out$jacobian
  colnames(jac) <- nm
  list(eta = out$eta, jacobian = jac)
}

#' @title Filter a Score-Driven Term Backwards
#' @name term_adjoint.GasTerm
#' @description
#' Runs the recursion of \code{\link{term_filter}} in reverse, returning the
#' derivative of a caller's objective with respect to the static predictor
#' it supplied and with respect to the sequence of scores it returned.
#' @param term A built \code{GasTerm}.
#' @param eta The static part of the predictor.
#' @param y The response, unused directly.
#' @param score,curvature The callbacks of \code{\link{term_filter}}.
#' @param psi The parameters, named as \code{\link{term_params}}.
#' @param g The direct derivative of the objective in the predictor the
#'   filter produced, one value per observation.
#' @param ... Unused.
#' @return A list with \code{deta} and \code{dscore}.
#' @keywords internal
S7::method(term_adjoint, GasTerm) <- function(term, eta, y, score, curvature,
                                              psi, g, ...) {
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
    return(.gas_adjoint_sub(term, eta, y, score, curvature, v, g))
  }
  # the predictor the recursion produced, which is where the callbacks are read
  e <- term_filter(term, eta, y, score, curvature, psi)$eta

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
      eb <- g[row] + sb[t] * curvature(e[row], row)
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
#' diagonal, is \code{\link[linkfunctions7]{d2linkinv}}; on the identity
#' both collapse to one and zero. The persistence reaches the coefficients
#' through two maps -- the link onto the partial autocorrelations and
#' Levinson-Durbin onto the coefficients -- so its second derivative
#' carries both a term in the map's own curvature and one in the link's.
#'
#' @param zeta The term's parameters on the unconstrained scale.
#' @param p,q The score and autoregressive orders.
#' @param links The links, as \code{\link{term_links}} gives them.
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
#' \code{\link{.gas_chart_derivs}} differentiated once more and contracted
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
#' @param links The links, as \code{\link{term_links}} gives them.
#' @param vz The direction, in the same coordinates as \code{zeta}.
#'
#' @return A list with \code{t_omega}, \code{t_a} and \code{t_b}, each a
#'   matrix or a list of matrices over the term's base coordinates.
#'
#' @seealso \code{\link{.gas_chart_derivs}}, \code{\link{gas_levinson3}}
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
#' the second derivative is accumulated on each group's ACTIVE SET rather
#' than as a square over all the unknowns: a development's coordinate
#' reaches only the groups where its column is not identically zero, so
#' with grouping indicators the active set has the same size whether the
#' panel has ten groups or a thousand. Measured, the full square cost
#' 0.39 s at 124 unknowns over 1600 rows and would have reached about
#' twelve minutes at five hundred groups.
#'
#' \code{blocks} is called with the row of the jacobian RESTRICTED to the
#' active set and with that set, and returns its pieces in the same
#' coordinates. A callback of the earlier three-argument shape is still
#' accepted and given the full row, its result being subset here; it costs
#' the quadratic allocation the restriction exists to avoid.
#' @param term A built \code{GasTerm}.
#' @param eta The static part of the predictor.
#' @param y The response, unused directly.
#' @param score,curvature The callbacks of \code{\link{term_filter}}.
#' @param psi The parameters on the PARAMETER scale, named as
#'   \code{\link{term_params}}.
#' @param g The weights the second derivative is contracted against.
#' @param seed The derivative of the static predictor in the unknowns.
#' @param blocks The model's own derivative pieces; see
#'   \code{\link{term_curvature}}.
#' @param ... Unused.
#' @return A list with \code{jacobian} and \code{curvature}.
#' @keywords internal
S7::method(term_curvature, GasTerm) <- function(term, eta, y, score,
                                                curvature, psi, g, seed,
                                                blocks, ...) {
  .gas_curvature_core(term, eta, y, score, curvature, psi, g, seed, blocks,
                      direction = NULL)
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
#' seeded by its second. Everything else it needs -- the directional
#' derivatives of \eqn{F}, \eqn{\Phi}, \eqn{\dot S} and \eqn{\ddot S} -- is a
#' contraction of a quantity the second-order recursion already carries, so
#' no second recursion is run and no three-index array is formed.
#'
#' The chart contributes its own third derivatives: the level and the
#' loadings through their scalar links, the persistence through
#' \code{\link{.gas_chart_derivs3}}, which composes
#' \code{\link{gas_levinson3}} with them.
#' @param term A built \code{GasTerm}.
#' @param eta The static part of the predictor.
#' @param y The response, unused directly.
#' @param score,curvature The callbacks of \code{\link{term_filter}}.
#' @param psi The parameters on the PARAMETER scale.
#' @param g The weights the third derivative is contracted against.
#' @param seed The derivative of the static predictor in the unknowns.
#' @param blocks The model's derivative pieces; see \code{\link{term_third}}.
#' @param direction The direction to contract against.
#' @param ... Unused.
#' @return A list with \code{jacobian}, \code{dphi} and \code{curvature}.
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
#' The body \code{\link{term_curvature}} and \code{\link{term_third}} share.
#' With \code{direction} \code{NULL} it propagates the first two derivatives
#' of the predictor; with a direction it propagates the third as well,
#' contracted against it.
#'
#' @details
#' The two orders are written here once rather than in a method each. The
#' third order's recursion reads \eqn{F}, \eqn{\Phi}, \eqn{\dot S} and
#' \eqn{\ddot S} at every lag, so a separate implementation would carry a
#' second copy of the first two orders, and the two would drift.
#'
#' @param term A built \code{GasTerm}.
#' @param eta The static part of the predictor.
#' @param y The response.
#' @param score,curvature The callbacks of \code{\link{term_filter}}.
#' @param psi The parameters on the parameter scale.
#' @param g The weights the contraction is taken against.
#' @param seed The derivative of the static predictor in the unknowns.
#' @param blocks The model's derivative pieces.
#' @param direction The direction, or \code{NULL} for the second order alone.
#'
#' @return A list with \code{jacobian} and \code{curvature}, and \code{dphi}
#'   where a direction was given.
#'
#' @keywords internal
.gas_curvature_core <- function(term, eta, y, score, curvature, psi, g, seed,
                                blocks, direction = NULL) {
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
                              seed, blocks, direction))
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
#' The loop \code{gas_filter_cpp()} replaces, kept so the compiled
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
#' @param score,curvature The callbacks of \code{\link{term_filter}}.
#'
#' @return A list with \code{eta} and \code{jacobian}.
#'
#' @keywords internal
gas_filter_r <- function(eta, order, p, q, omega, a, b, db, f0, df0,
                         i_a, np, score, curvature) {
  n <- length(eta)
  eta_out <- numeric(n)
  jac <- matrix(0, n, np)

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
      ds[t, ] <- curvature(e_t, rows[t]) * df[t, ]
    }

    eta_out[rows] <- eta[rows] + f
    jac[rows, ] <- df
  }
  list(eta = eta_out, jacobian = jac)
}
