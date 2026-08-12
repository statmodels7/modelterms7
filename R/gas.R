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
#' @param deviations Which parameters carry a deviation per group.
#' @param penalty_kind The penalty on the deviations, if any.
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
    deviations = S7::class_any,
    penalty_kind = S7::class_character,
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
#' \eqn{a_1, \dots, a_p}, and the persistence. The persistence is carried
#' by \strong{partial autocorrelations} rather than by the coefficients
#' \eqn{b_j}: the stationary region of an autoregression is not a box, so
#' no collection of scalar links covers it, while the partial
#' autocorrelations each range over \eqn{(-1, 1)} independently and the
#' Levinson-Durbin recursion carries them onto the coefficients
#' bijectively. At \eqn{q = 1} the two coincide. The coordinate is named
#' for the chart it lives on, \code{pacf1} and so on, following the
#' convention of \pkg{parameters7}.
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
#' \subsection{A population value and a deviation per group}{
#' By default every group of a panel is filtered with the same parameters.
#' \code{deviations} gives each group its own, written as a population
#' value and a departure from it,
#' \deqn{\psi_{j,i} = g_j^{-1}\!\left(g_j(\psi_j) + \delta_{j,i}\right),}
#' the deviation acting on the unconstrained scale of the chart the
#' parameter lives on, so that a persistence stays inside \eqn{(-1, 1)}
#' whatever the deviation is. The deviations are parameters of the term,
#' named \code{omega.dev.1} and so on after the parameter and the level,
#' and they carry the identity link, being unconstrained already.
#'
#' They are parameters and not a penalty on the per-group values through a
#' difference matrix, which is what the same model looks like written the
#' other way. The difference decides what can be fitted: a penalty over a
#' general map is the generalized-lasso problem, whose proximal operator
#' does not split by coordinate, whereas a deviation named as a coordinate
#' is reached by a soft threshold and by a coordinate descent unchanged.
#' \code{penalty} shrinks them towards zero, which is towards a panel that
#' is homogeneous in that parameter, and \code{"lasso"} sets the
#' deviations of the groups that do not need one exactly to it.
#'
#' The penalty is also what identifies them. A parameter and its \eqn{m}
#' deviations are \eqn{m+1} numbers describing \eqn{m} group values, so
#' adding a constant to \eqn{g_j(\psi_j)} and subtracting it from every
#' \eqn{\delta_{j,i}} leaves the filter and its likelihood exactly
#' unchanged: without a penalty on the deviations the likelihood is flat
#' along one direction per parameter carrying them. This is the
#' parametrization of a random effect, and it is identified in the same way
#' -- there by a variance component, here by the penalty, which selects the
#' deviations of smallest size among those that describe the same panel.
#' \code{penalty = "none"} is therefore for reading a filter at given
#' parameters rather than for fitting one.
#'
#' The parameters a specification reports are the population ones alone:
#' how many groups there are is a property of the data, so the deviations
#' appear once the term is built.
#' }
#'
#' @param p The number of score lags. Defaults to 1.
#' @param q The number of autoregressive lags. Defaults to 1.
#' @param by An optional grouping variable, evaluated in the data; each
#'   group is filtered independently, from its own starting level.
#' @param time An optional ordering variable, evaluated in the data.
#' @param deviations Whether each group carries a deviation from the
#'   population parameters: \code{FALSE} (default), \code{TRUE} for every
#'   parameter, or a character vector naming the parameters that carry
#'   one. Requires \code{by}.
#' @param penalty One of \code{"none"} (default), \code{"lasso"} or
#'   \code{"ridge"}, applied to the deviations. Requires them.
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
gas <- function(p = 1, q = 1, by = NULL, time = NULL, deviations = FALSE,
                penalty = c("none", "lasso", "ridge"), label = "gas") {
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
  penalty <- match.arg(penalty)
  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !nzchar(label)) {
    stop("'label' must be a single non-empty character string.",
         call. = FALSE)
  }
  by <- substitute(by)
  base <- .gas_base_params(p, q)
  if (!isFALSE(deviations)) {
    if (!isTRUE(deviations) &&
        (!is.character(deviations) || !length(deviations) ||
         anyNA(deviations))) {
      stop(paste("'deviations' must be FALSE, TRUE, or a character vector",
                 "naming the parameters that carry one."), call. = FALSE)
    }
    if (is.null(by)) {
      stop("'deviations' needs 'by': a deviation is a departure per group.",
           call. = FALSE)
    }
    bad <- setdiff(if (isTRUE(deviations)) character(0) else deviations, base)
    if (length(bad)) {
      stop(sprintf("'deviations' names '%s'; the parameters are %s.",
                   bad[1L], paste(base, collapse = ", ")), call. = FALSE)
    }
  } else if (penalty != "none") {
    stop(paste("'penalty' reaches the deviations of a panel, so it needs",
               "'deviations'; the population parameters of a filter are not",
               "shrunk towards zero."), call. = FALSE)
  }
  GasTerm(label = label, p = p, q = q,
          by = by, time = substitute(time),
          deviations = deviations, penalty_kind = penalty,
          blueprint = list())
}

# the parameters of the filter itself, before any deviation
.gas_base_params <- function(p, q) {
  c("omega",
    if (p > 0L) paste0("a", seq_len(p)),
    if (q > 0L) paste0("pacf", seq_len(q)))
}

# which of them carry a deviation per group
.gas_dev_params <- function(term) {
  d <- term@deviations
  if (isFALSE(d)) return(character(0))
  if (isTRUE(d)) return(.gas_base_params(term@p, term@q))
  d
}

S7::method(term_params, GasTerm) <- function(term, ...) {
  base <- .gas_base_params(term@p, term@q)
  dv <- .gas_dev_params(term)
  levs <- term@blueprint$levels
  if (!length(dv) || is.null(levs)) return(base)
  c(base, unlist(lapply(dv, function(p) paste(p, "dev", levs, sep = "."))))
}

#' @title The Level of a Score-Driven Term
#' @name term_level_param.GasTerm
#' @description
#' \code{"omega"}, which adds a constant to the equation's predictor and is
#' therefore the direction an intercept there also spans.
#' @param term A \code{\link{GasTerm}}.
#' @param ... Unused.
#' @return A single string.
#' @keywords internal
S7::method(term_level_param, GasTerm) <- function(term, ...) "omega"

S7::method(term_links, GasTerm) <- function(term, ...) {
  base <- .gas_base_params(term@p, term@q)
  nm <- term_params(term)
  stats::setNames(lapply(nm, function(p) {
    # a deviation is already unconstrained: it acts on the scale the
    # population parameter's own link carries it to
    if (!(p %in% base)) linkfunctions7::identity_link()
    else if (startsWith(p, "pacf")) linkfunctions7::rhobit_link()
    else linkfunctions7::identity_link()
  }), nm)
}

#' @title Penalties of a Score-Driven Term
#' @name term_penalties.GasTerm
#' @description
#' One entry per parameter carrying deviations, named after it and covering
#' its deviations across the groups. The population parameters are
#' unpenalized, and the list is empty when \code{penalty = "none"}, and for
#' a specification, whose deviations do not exist until the data say how
#' many groups there are.
#' @param term A built \code{\link{GasTerm}}.
#' @param ... Unused.
#' @return A list of entries, as \code{\link{term_penalties}} documents.
#' @keywords internal
S7::method(term_penalties, GasTerm) <- function(term, ...) {
  if (identical(term@penalty_kind, "none")) return(list())
  levs <- term@blueprint$levels
  if (is.null(levs)) return(list())
  base <- .gas_base_params(term@p, term@q)
  dv <- .gas_dev_params(term)
  m <- length(levs)
  factory <- .penalty_factory(term@penalty_kind)
  lapply(seq_along(dv), function(i) {
    list(name = dv[i], index = length(base) + (i - 1L) * m + seq_len(m),
         penalty = factory(m))
  })
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

# One group's parameter values and the two chain factors the jacobian needs:
# with psi_i = g^-1(g(psi) + delta), the derivative in the population value
# is g^-1'(g(psi) + delta) g'(psi) and the one in the deviation is the same
# without the second factor. At delta = 0 the first is 1, by the inverse
# function theorem, which is what makes an unused deviation cost nothing.
.gas_group_values <- function(base_psi, links, dv, delta) {
  v <- base_psi
  dpop <- rep(1, length(base_psi))
  ddev <- stats::setNames(numeric(length(dv)), dv)
  for (p in dv) {
    lk <- links[[p]]
    e <- linkfunctions7::linkfun(lk, base_psi[[p]]) + delta[[p]]
    v[[p]] <- linkfunctions7::linkinv(lk, e)
    di <- linkfunctions7::dlinkinv(lk, e)
    ddev[[p]] <- di
    dpop[[which(names(base_psi) == p)]] <-
      di * linkfunctions7::dlinkfun(lk, base_psi[[p]])
  }
  list(value = v, dpop = dpop, ddev = ddev)
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
  p <- term@p
  q <- term@q
  base <- .gas_base_params(p, q)
  dv <- .gas_dev_params(term)
  base_psi <- psi[base]

  if (!length(dv)) {
    cf <- .gas_coefs(base_psi, p, q)
    out <- gas_filter_cpp(eta, bp$order, p, q, cf$omega, cf$a, cf$b, cf$db,
                          cf$f0, cf$df0, cf$i_a, cf$np, score, curvature)
    jac <- out$jacobian
    colnames(jac) <- nm
    return(list(eta = out$eta, jacobian = jac))
  }

  # With deviations each group runs on parameters of its own, so the filter
  # is called once per group and the columns of its jacobian are chained
  # onto the population values and scattered into that group's deviations.
  # A group's rows reach no other group's deviation, the groups being
  # filtered independently.
  levs <- bp$levels
  m <- length(levs)
  links <- term_links(term)
  nb <- length(base)
  eta_out <- numeric(bp$n)
  jac <- matrix(0, bp$n, length(nm))
  for (l in seq_len(m)) {
    delta <- stats::setNames(
      as.list(psi[nb + (seq_along(dv) - 1L) * m + l]), dv)
    gv <- .gas_group_values(base_psi, links, dv, delta)
    cf <- .gas_coefs(gv$value, p, q)
    rows <- bp$order[[l]]
    out <- gas_filter_cpp(eta, list(rows), p, q, cf$omega, cf$a, cf$b,
                          cf$db, cf$f0, cf$df0, cf$i_a, cf$np, score,
                          curvature)
    eta_out[rows] <- out$eta[rows]
    g <- out$jacobian[rows, , drop = FALSE]
    jac[rows, seq_len(nb)] <- g * rep(gv$dpop, each = length(rows))
    for (i in seq_along(dv)) {
      j <- match(dv[i], base)
      jac[rows, nb + (i - 1L) * m + l] <- g[, j] * gv$ddev[[dv[i]]]
    }
  }
  colnames(jac) <- nm
  list(eta = eta_out, jacobian = jac)
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
  # the predictor the recursion produced, which is where the callbacks are read
  e <- term_filter(term, eta, y, score, curvature, psi)$eta

  p <- term@p
  q <- term@q
  base <- .gas_base_params(p, q)
  dv <- .gas_dev_params(term)
  base_psi <- v[base]
  links <- if (length(dv)) term_links(term) else NULL
  levs <- bp$levels
  nb <- length(base)
  m <- length(levs)

  deta <- numeric(bp$n)
  dscore <- numeric(bp$n)
  for (l in seq_along(bp$order)) {
    rows <- bp$order[[l]]
    vals <- if (!length(dv)) base_psi else {
      delta <- stats::setNames(
        as.list(v[nb + (seq_along(dv) - 1L) * m + l]), dv)
      .gas_group_values(base_psi, links, dv, delta)$value
    }
    cf <- .gas_coefs(vals, p, q)
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
#' The level and the loadings carry the identity link, so their first
#' derivative is one and their second is zero. The persistence reaches the
#' coefficients through two maps -- the link onto the partial
#' autocorrelations and Levinson-Durbin onto the coefficients -- so its
#' second derivative carries both a term in the map's own curvature and one
#' in the link's, which is where \code{\link[linkfunctions7]{d2linkinv}}
#' enters.
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

  omega <- zeta[[i_om]]
  a <- if (p > 0L) zeta[i_a] else numeric(0)
  d_omega <- numeric(np)
  d_omega[i_om] <- 1
  d_a <- lapply(seq_len(max(p, 1L)), function(i) {
    v <- numeric(np)
    if (p > 0L) v[i_a[i]] <- 1
    v
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
  list(omega = omega, a = a, b = b, d_omega = d_omega, d_a = d_a,
       d_b = d_b, h_b = h_b, np = np)
}

#' @title Second Derivatives of a Score-Driven Predictor
#' @name term_curvature.GasTerm
#' @description
#' The forward Jacobian of the filter's predictor in a caller's unknowns and
#' the second derivative contracted against the caller's weights, both
#' propagated through the recursion beside the state.
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
  bp <- term@blueprint
  if (!length(bp)) {
    stop("the term has not been built; call term_build(term, data) first.",
         call. = FALSE)
  }
  if (length(.gas_dev_params(term))) {
    stop(paste("term_curvature() does not carry deviations yet: the",
               "per-group chain adds a factor to every derivative and is",
               "not written."), call. = FALSE)
  }
  nm <- term_params(term)
  links <- term_links(term)
  psiv <- unlist(psi[nm])
  # the recursion is driven by the parameters, and the caller's unknowns
  # reach them through the links, so the chart is differentiated on the
  # UNCONSTRAINED scale, which is what the model estimates
  zeta <- vapply(nm, function(j)
    linkfunctions7::linkfun(links[[j]], psiv[[j]]), numeric(1))
  ch <- .gas_chart_derivs(zeta, term@p, term@q, links)

  seed <- as.matrix(seed)
  m <- ncol(seed)
  np <- ch$np
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
  lift <- function(v) {
    out <- numeric(m)
    out[zcol] <- v
    out
  }
  lift2 <- function(h) {
    out <- matrix(0, m, m)
    out[zcol, zcol] <- h
    out
  }
  om_u <- lift(ch$d_omega)
  a_u <- lapply(ch$d_a, lift)
  b_u <- lapply(ch$d_b, lift)
  b_uu <- lapply(ch$h_b, lift2)

  p <- term@p
  q <- term@q
  a <- ch$a
  b <- ch$b
  sb <- if (q > 0L) sum(b) else 0
  if (abs(1 - sb) < 1e-10) {
    stop("the autoregressive polynomial is at the unit root; the filter has no starting level.",
         call. = FALSE)
  }
  # the starting level and its two derivatives: f0 = omega/(1 - sum b)
  db_sum <- if (q > 0L) Reduce(`+`, b_u) else numeric(m)
  f0 <- ch$omega / (1 - sb)
  f0_u <- om_u / (1 - sb) + ch$omega * db_sum / (1 - sb)^2
  f0_uu <- (outer(om_u, db_sum) + outer(db_sum, om_u)) / (1 - sb)^2 +
    2 * ch$omega * outer(db_sum, db_sum) / (1 - sb)^3
  if (q > 0L) {
    hb_sum <- Reduce(`+`, b_uu)
    f0_uu <- f0_uu + ch$omega * hb_sum / (1 - sb)^2
  }

  D <- matrix(0, bp$n, m)
  W <- matrix(0, m, m)
  for (rows in bp$order) {
    k <- length(rows)
    f <- numeric(k)
    s <- numeric(k)
    F_ <- vector("list", k)
    Phi <- vector("list", k)
    Sd <- vector("list", k)
    Sdd <- vector("list", k)
    for (t in seq_len(k)) {
      row <- rows[t]
      ft <- ch$omega
      Ft <- om_u
      Pt <- matrix(0, m, m)
      if (p > 0L) {
        for (i in seq_len(p)) {
          lag <- t - i
          s_l <- if (lag >= 1L) s[lag] else 0
          Sd_l <- if (lag >= 1L) Sd[[lag]] else numeric(m)
          Sdd_l <- if (lag >= 1L) Sdd[[lag]] else matrix(0, m, m)
          ft <- ft + a[[i]] * s_l
          Ft <- Ft + a[[i]] * Sd_l + s_l * a_u[[i]]
          Pt <- Pt + a[[i]] * Sdd_l +
            outer(a_u[[i]], Sd_l) + outer(Sd_l, a_u[[i]])
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
        }
      }
      f[t] <- ft
      F_[[t]] <- Ft
      Phi[[t]] <- Pt

      e_t <- eta[row] + ft
      Dt <- seed[row, ] + Ft
      D[row, ] <- Dt
      W <- W + g[row] * Pt

      s[t] <- score(e_t, row)
      cv <- curvature(e_t, row)
      bl <- blocks(e_t, row, Dt)
      Sd[[t]] <- cv * Dt + bl$cross
      Sdd[[t]] <- cv * Pt + bl$M
    }
  }
  # The matrix is symmetric and the accumulation is not: an entry and its
  # transpose collect the same terms in a different ORDER, and (x + p) + q
  # is not (x + q) + p in floating point. The gap is of the order of the
  # rounding, and a caller about to factor this wants an exactly symmetric
  # matrix rather than one that is nearly so.
  W <- (W + t(W)) / 2
  list(jacobian = D, curvature = W)
}

S7::method(print, GasTerm) <- function(x, ...) {
  built <- length(x@blueprint) > 0L
  cat(sprintf("<GasTerm> '%s': score-driven, p = %d, q = %d%s\n",
              x@label, x@p, x@q,
              if (built) sprintf("; %d group(s)", length(x@blueprint$order))
              else " (specification)"))
  dv <- .gas_dev_params(x)
  if (length(dv)) {
    cat(sprintf("  deviations on: %s%s\n", paste(dv, collapse = ", "),
                if (x@penalty_kind != "none")
                  sprintf("; %s", x@penalty_kind) else ""))
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
