#' @include term_classes.R generics.R structural.R
NULL

#' @title Log-Likelihood Contributions of a Structural Term
#'
#' @description
#' The per-observation log-likelihood contributions a structural term
#' produces, with their derivatives in the term's own parameters. This is
#' the second shape the structural branch takes, beside
#' \code{\link{term_filter}}: a term that shifts the predictor implements
#' the filter, and one that rewrites the likelihood itself -- a mixture
#' over latent states, say -- implements this, because its contribution is
#' not a predictor and cannot be reported as one.
#'
#' @details
#' The contribution of observation \eqn{t} is the logarithm of its one-step
#' predictive density given everything before it,
#'
#' \deqn{\ell_t(\psi) = \log f(y_t \mid y_1, \dots, y_{t-1}; \psi),
#'   \qquad \sum_{t=1}^{n} \ell_t(\psi)
#'     = \log f(y_1, \dots, y_n; \psi),}
#'
#' so the vector returned sums to the term's log-likelihood by the chain
#' rule of probability whatever the dependence between observations. For
#' \code{\link{regime}} it is the normalizing constant of the forward
#' recursion, \eqn{\ell_t = \log \sum_{k} \pi_{t \mid t-1, k}
#' f(y_t \mid S_t = k)}, and the Jacobian
#' \eqn{\partial \ell_t / \partial \psi_j} is propagated beside the
#' filtered distribution rather than differenced.
#'
#' @param term A built structural term.
#' @param eta The static part of the linear predictor.
#' @param y The response.
#' @param logdens A function of a predictor value and a row index,
#'   returning the log-density of that observation at that predictor.
#' @param score A function of the same two arguments returning the
#'   derivative of that log-density with respect to the predictor.
#' @param psi The term's parameters, named as \code{\link{term_params}}.
#' @param ... Passed to methods.
#'
#' @return A list with \code{loglik}, one contribution per observation
#'   summing to the term's log-likelihood, and \code{jacobian}, an
#'   \code{n} by \code{length(psi)} matrix of its derivatives.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:40, y = c(rnorm(20), rnorm(20, 3)))
#' term <- term_build(regime(2, time = t), dd)
#' out <- term_loglik(term, rep(0, 40), dd$y,
#'                    logdens = function(e, i) dnorm(dd$y[i], e, log = TRUE),
#'                    score = function(e, i) dd$y[i] - e,
#'                    psi = list(level1 = 0, gap2 = 3,
#'                               alr1.1 = 2, alr2.1 = -2))
#' sum(out$loglik)
#'
#' @seealso \code{\link{regime}}, \code{\link{term_filter}}
#' @export
term_loglik <- S7::new_generic("term_loglik", "term",
  function(term, eta, y, logdens, score, psi, ...) S7::S7_dispatch())

S7::method(term_loglik, structural_term) <- function(term, eta, y, logdens,
                                                     score, psi, ...) {
  stop(sprintf("the term class '%s' does not implement term_loglik().",
               attr(S7::S7_class(term), "name")), call. = FALSE)
}

#' @title S7 Class for Markov Regime Terms
#' @name RegimeTerm
#'
#' @description
#' A subclass of \code{\link{structural_term}} for a latent Markov chain
#' of regimes, each shifting the linear predictor by a level of its own.
#' Constructed by \code{\link{regime}}.
#'
#' @inheritParams model_term
#' @param k The number of regimes.
#' @param by An optional grouping expression, run independently.
#' @param time An optional ordering expression.
#' @param chain The \pkg{parameters7} transition matrix.
#' @param blueprint The resolved ordering and grouping.
#'
#' @return An object of class \code{RegimeTerm}.
#'
#' @seealso \code{\link{regime}}
#' @examples
#' S7::S7_inherits(regime(2), RegimeTerm)
#' @export
RegimeTerm <- S7::new_class(
  name = "RegimeTerm",
  parent = structural_term,
  properties = list(
    k = S7::class_integer,
    by = S7::class_any,
    time = S7::class_any,
    chain = S7::class_any,
    blueprint = S7::class_list
  )
)

#' Markov Regime Switching
#'
#' @description
#' A latent Markov chain of \eqn{K} regimes, each shifting the linear
#' predictor by a level of its own (\cite{hamilton1989}). The likelihood
#' is the mixture over the unobserved state path, evaluated by the forward
#' recursion, and it is built from whatever density the model carries: the
#' term supplies the chain and the levels, the distribution supplies
#' everything else.
#'
#' @details
#' Writing \eqn{f_j(t)} for the density of observation \eqn{t} at the
#' predictor shifted by the level of regime \eqn{j}, the forward recursion
#' is
#' \deqn{\tilde\alpha_t(j) = f_j(t) \sum_i \alpha_{t-1}(i) P_{ij},
#'   \qquad c_t = \sum_j \tilde\alpha_t(j), \qquad
#'   \alpha_t = \tilde\alpha_t / c_t,}
#' started at the chain's stationary distribution, and the log-likelihood
#' is \eqn{\sum_t \log c_t}. Normalizing at every step is what keeps the
#' recursion representable: the unnormalized quantities are products of
#' \eqn{t} densities, so they decay geometrically and reach zero in
#' double precision on a series of a few hundred observations. The
#' contributions \eqn{\log c_t} are what \code{\link{term_loglik}}
#' returns, one per observation, together with their exact derivatives,
#' propagated through the recursion beside the state.
#'
#' This is the second dynamic model of the package and it is the
#' complement of the first: \code{\link{gas}} is driven by the score of
#' the density and moves continuously, while a regime chain moves in
#' jumps between a finite number of states. Both are built from the
#' density rather than from an error structure, so both apply to any
#' family the model carries.
#'
#' \subsection{The parameters and their charts}{
#' The levels are \strong{ordered by construction}: the first is free and
#' each of the others is the previous one plus a positive gap, carried on
#' a log link. Without an ordering the regimes are exchangeable and the
#' likelihood has \eqn{K!} identical maxima, which is not a hard problem
#' to fit but is one whose answer cannot be reported. The transition
#' matrix is \code{\link[parameters7]{transition_matrix}}, whose free
#' values are the additive log-ratios of each row, so every row is a
#' probability vector by construction.
#'
#' The initial distribution is the chain's stationary one, which costs no
#' parameters and whose derivative is obtained from the linear system it
#' solves.
#' }
#'
#' @param k The number of regimes, at least 2.
#' @param by An optional grouping variable; each group runs its own
#'   recursion from the stationary distribution.
#' @param time An optional ordering variable.
#' @param label A single non-empty string naming the term.
#'
#' @return An object of class \code{\link{RegimeTerm}} (a specification;
#'   see \code{\link{term_build}}).
#'
#' @references
#' Hamilton, J. D. (1989). A new approach to the economic analysis of
#' nonstationary time series and the business cycle. \emph{Econometrica},
#' 57(2), 357--384.
#'
#' @examples
#' term_params(regime(2))
#'
#' @export
regime <- function(k = 2, by = NULL, time = NULL, label = "regime") {
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 2 ||
      k != round(k)) {
    stop("'k' must be a single integer of at least 2.", call. = FALSE)
  }
  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !nzchar(label)) {
    stop("'label' must be a single non-empty character string.",
         call. = FALSE)
  }
  k <- as.integer(k)
  RegimeTerm(label = label, k = k,
             by = substitute(by), time = substitute(time),
             chain = parameters7::transition_matrix(k),
             blueprint = list())
}

S7::method(term_params, RegimeTerm) <- function(term, ...) {
  c("level1",
    if (term@k > 1L) paste0("gap", seq.int(2L, term@k)),
    term@chain@free_names)
}

S7::method(term_links, RegimeTerm) <- function(term, ...) {
  nm <- term_params(term)
  stats::setNames(lapply(nm, function(p) {
    if (startsWith(p, "gap")) linkfunctions7::log_link()
    else linkfunctions7::identity_link()
  }), nm)
}

S7::method(term_build, RegimeTerm) <- function(term, data, ...) {
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
  ord <- split(order(grp, tm), grp[order(grp, tm)])
  term@blueprint <- list(order = ord, n = n)
  term
}

#' The Stationary Distribution of a Chain, and Its Derivative
#'
#' @description
#' The row vector \eqn{\delta} solving \eqn{\delta P = \delta} with
#' \eqn{\sum_j \delta_j = 1}, together with its derivative in whatever the
#' transition matrix depends on.
#'
#' @details
#' Differentiating \eqn{\delta(I - P) = 0} under the normalization gives
#' \eqn{d\delta\,(I - P) = \delta\,dP} with \eqn{\sum_j d\delta_j = 0}, so
#' both the value and the derivative come from the same linear system with
#' one column replaced by the normalization.
#'
#' @param P A row-stochastic matrix.
#' @param dP A list of derivative matrices of \code{P}.
#'
#' @return A list with \code{delta} and \code{ddelta}, the latter one row
#'   per element of \code{dP}.
#'
#' @keywords internal
regime_stationary <- function(P, dP) {
  k <- nrow(P)
  A <- diag(k) - P
  A[, k] <- 1
  delta <- as.numeric(solve(t(A), c(rep(0, k - 1L), 1)))
  ddelta <- lapply(dP, function(D) {
    rhs <- as.numeric(delta %*% D)
    rhs[k] <- 0
    as.numeric(solve(t(A), rhs))
  })
  list(delta = delta, ddelta = ddelta)
}

#' @title Log-Likelihood of a Regime Term
#' @name term_loglik.RegimeTerm
#' @description
#' Runs the forward recursion over each group in time order, normalizing
#' at every step, and returns the per-observation contributions with their
#' exact derivatives.
#' @param term A built \code{RegimeTerm}.
#' @param eta The static predictor.
#' @param y The response, reaching the recursion through the two callbacks.
#' @param logdens,score The log-density and its derivative in the predictor.
#' @param psi The parameters, named as \code{\link{term_params}}.
#' @param ... Unused.
#' @return A list with \code{loglik} and \code{jacobian}.
#' @keywords internal
S7::method(term_loglik, RegimeTerm) <- function(term, eta, y, logdens, score,
                                                psi, ...) {
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
  k <- term@k
  np <- length(nm)
  n <- bp$n

  # the levels, ordered by construction, and their derivatives
  gaps <- if (k > 1L) v[paste0("gap", seq.int(2L, k))] else numeric(0)
  if (any(gaps <= 0)) {
    stop("every gap must be positive: the levels are ordered by construction.",
         call. = FALSE)
  }
  mu <- cumsum(c(v[["level1"]], gaps))
  dmu <- matrix(0, k, np)
  dmu[, 1L] <- 1
  for (j in seq_len(k)) {
    if (j > 1L) dmu[j, 1L + seq_len(j - 1L)] <- 1
  }

  # the chain and its derivatives, from parameters7
  n_lev <- k
  i_tr <- (1L + (k - 1L)) + seq_len(k * (k - 1L))
  eta_tr <- v[term@chain@free_names]
  P <- parameters7::param_value(term@chain, eta_tr)
  dP_list <- parameters7::param_d1(term@chain, eta_tr)
  # widened to the term's own parameter vector: zero outside the chain's slots
  dP <- vector("list", np)
  for (i in seq_len(np)) dP[[i]] <- matrix(0, k, k)
  for (i in seq_along(term@chain@free_names)) {
    dP[[i_tr[i]]] <- dP_list[[term@chain@free_names[i]]]
  }

  st <- regime_stationary(P, dP)

  # A regime shifts the predictor by a level of its own, so the density
  # and the score of every observation under every regime are known
  # before the recursion starts: k vectorized calls replace the 2nk
  # scalar ones the loop used to make. The closures must therefore
  # accept the whole index vector, which the generic's contract -- one
  # value per observation -- already asks of them.
  idx <- seq_len(n)
  LF <- matrix(0, n, k)
  SC <- matrix(0, n, k)
  for (jj in seq_len(k)) {
    e <- eta + mu[jj]
    a <- as.numeric(logdens(e, idx))
    b <- as.numeric(score(e, idx))
    if (length(a) != n || length(b) != n) {
      stop(paste0("'logdens' and 'score' must return one value per ",
                  "observation when given the whole index vector."),
           call. = FALSE)
    }
    LF[, jj] <- a
    SC[, jj] <- b
  }

  ddm <- do.call(rbind, st$ddelta)
  out <- regime_forward_cpp(bp$order, LF, SC, dmu, unclass(P), dP,
                            st$delta, ddm)
  loglik <- out$loglik
  jac <- out$jacobian

  colnames(jac) <- nm
  list(loglik = loglik, jacobian = jac)
}

S7::method(print, RegimeTerm) <- function(x, ...) {
  built <- length(x@blueprint) > 0L
  cat(sprintf("<RegimeTerm> '%s': %d regimes%s\n", x@label, x@k,
              if (built) sprintf("; %d group(s)", length(x@blueprint$order))
              else " (specification)"))
  cat("  parameters: ", paste(term_params(x), collapse = ", "), "\n", sep = "")
  invisible(x)
}

# The recursion the compiled kernel replaces, kept as the twin the tests
# hold it to. The arithmetic is written the way the derivation reads:
# the state is normalized at every step and the density factored by its
# largest value, without either of which the unnormalized quantities
# underflow within a few hundred observations.
.regime_forward_r <- function(order, LF, SC, dmu, P, dP, delta, ddelta) {
  n <- nrow(LF); k <- ncol(LF); np <- ncol(dmu)
  loglik <- numeric(n)
  jac <- matrix(0, n, np)
  for (rows in order) {
    a <- delta
    da <- ddelta
    for (t in seq_along(rows)) {
      row <- rows[t]
      lf <- LF[row, ]
      mx <- max(lf)
      w <- exp(lf - mx)
      sc <- SC[row, ]
      # element (i, j) is w_j * s_j * dmu[j, i]: the recycling has to run
      # DOWN the columns of the np x k result, which rep(..., each = np)
      # does and a bare vector product does not
      dw <- t(dmu) * rep(w * sc, each = np)

      if (t == 1L) {
        pred <- a
        dpred <- da
      } else {
        pred <- as.numeric(a %*% P)
        dpred <- da %*% P
        for (i in seq_len(np)) {
          dpred[i, ] <- dpred[i, ] + as.numeric(a %*% dP[[i]])
        }
      }

      atil <- w * pred
      datil <- dw * rep(pred, each = np) + dpred * rep(w, each = np)
      ct <- sum(atil)
      dct <- rowSums(datil)

      loglik[row] <- log(ct) + mx
      jac[row, ] <- dct / ct

      a <- atil / ct
      da <- (datil - outer(dct, a)) / ct
    }
  }
  list(loglik = loglik, jacobian = jac)
}
