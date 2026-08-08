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
#' @param p The number of score lags. Defaults to 1.
#' @param q The number of autoregressive lags. Defaults to 1.
#' @param by An optional grouping variable, evaluated in the data; each
#'   group is filtered independently, from its own starting level.
#' @param time An optional ordering variable, evaluated in the data.
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
#' @export
gas <- function(p = 1, q = 1, by = NULL, time = NULL, label = "gas") {
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
  GasTerm(label = label, p = p, q = q,
          by = substitute(by), time = substitute(time),
          blueprint = list())
}

S7::method(term_params, GasTerm) <- function(term, ...) {
  c("omega",
    if (term@p > 0L) paste0("a", seq_len(term@p)),
    if (term@q > 0L) paste0("pacf", seq_len(term@q)))
}

S7::method(term_links, GasTerm) <- function(term, ...) {
  nm <- term_params(term)
  stats::setNames(lapply(nm, function(p) {
    if (startsWith(p, "pacf")) linkfunctions7::rhobit_link()
    else linkfunctions7::identity_link()
  }), nm)
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

S7::method(term_build, GasTerm) <- function(term, data, ...) {
  n <- nrow(data)
  grp <- if (is.null(term@by)) rep(1L, n) else {
    as.integer(factor(eval(term@by, data, baseenv())))
  }
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
  term@blueprint <- list(order = ord, n = n,
                         by = term@by, time = term@time)
  term
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
  n <- bp$n
  p <- term@p
  q <- term@q
  np <- length(psi)

  omega <- psi[["omega"]]
  a <- if (p > 0L) psi[paste0("a", seq_len(p))] else numeric(0)
  ld <- gas_levinson(if (q > 0L) psi[paste0("pacf", seq_len(q))] else numeric(0))
  b <- ld$phi

  # d omega, d a_i and d b_j with respect to psi: constants of the chart
  i_om <- 1L
  i_a <- if (p > 0L) 1L + seq_len(p) else integer(0)
  i_pa <- if (q > 0L) 1L + p + seq_len(q) else integer(0)
  db <- matrix(0, max(q, 1L), np)
  if (q > 0L) db[seq_len(q), i_pa] <- ld$jacobian

  # the level the recursion starts from, and its derivative
  sb <- if (q > 0L) sum(b) else 0
  if (abs(1 - sb) < 1e-10) {
    stop("the autoregressive polynomial is at the unit root; the filter has no starting level.",
         call. = FALSE)
  }
  f0 <- omega / (1 - sb)
  df0 <- numeric(np)
  df0[i_om] <- 1 / (1 - sb)
  if (q > 0L) df0 <- df0 + (omega / (1 - sb)^2) * colSums(db[seq_len(q), , drop = FALSE])

  out <- gas_filter_cpp(eta, bp$order, p, q, omega, a, b, db, f0, df0,
                        i_a, np, score, curvature)
  eta_out <- out$eta
  jac <- out$jacobian

  colnames(jac) <- nm
  list(eta = eta_out, jacobian = jac)
}

S7::method(print, GasTerm) <- function(x, ...) {
  built <- length(x@blueprint) > 0L
  cat(sprintf("<GasTerm> '%s': score-driven, p = %d, q = %d%s\n",
              x@label, x@p, x@q,
              if (built) sprintf("; %d group(s)", length(x@blueprint$order))
              else " (specification)"))
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
