#' @include term_classes.R generics.R structural.R
NULL

#' @title Log-Likelihood Contributions of a Structural Term
#'
#' @description
#' The per-observation log-likelihood contributions a structural term
#' produces, with their derivatives in the term's own parameters. This is the
#' second shape the structural branch takes, beside [term_filter()]: a term
#' that shifts the predictor implements the filter, and one that rewrites the
#' likelihood itself implements this, its contribution being no predictor at
#' all.
#'
#' @details
#' # What the contributions are
#'
#' The contribution of observation \eqn{t} is the logarithm of its one-step
#' predictive density given everything before it,
#'
#' \deqn{\ell_t(\psi) = \log f(y_t \mid y_1, \dots, y_{t-1}; \psi),
#'   \qquad \sum_{t=1}^{n} \ell_t(\psi)
#'     = \log f(y_1, \dots, y_n; \psi),}
#'
#' so the vector returned sums to the term's log-likelihood by the chain rule
#' of probability, whatever the dependence between observations.
#'
#' For [regime()] it is the normalizing constant of the forward recursion,
#' \eqn{\ell_t = \log \sum_{k} \pi_{t \mid t-1, k} f(y_t \mid S_t = k)}, and
#' the Jacobian \eqn{\partial \ell_t / \partial \psi_j} is propagated beside
#' the filtered distribution, never differenced.
#'
#' # How the model reaches it
#'
#' The term knows the chain and the levels and nothing about the family. The
#' two callbacks are how the model's own density enters: `logdens` returns the
#' log-density of one observation at a given predictor and `score` its
#' derivative in that predictor. Both are called with the predictor shifted by
#' each regime's level, so a term of \eqn{K} regimes evaluates them \eqn{K}
#' times per observation.
#'
#' Unlike [term_filter()]'s callbacks these do **not** depend on the state the
#' recursion has reached: a regime shifts a predictor known in advance. That
#' is why the density and the score of every observation under every regime
#' can be computed once, vectorized, before the recursion starts.
#'
#' The method on [structural_term()] throws, naming the class: a term of the
#' filter shape implements [term_filter()] instead, and a fitting layer tells
#' the two apart by which one answers.
#'
#' @param term A built structural term.
#' @param eta The static part of the linear predictor, one value per
#'   observation.
#' @param y The response. It reaches the recursion through the callbacks; the
#'   argument is passed for methods that need it directly.
#' @param logdens A function `(e, i)` returning the log-density of
#'   observation `i` at predictor `e`.
#' @param score A function `(e, i)` returning the derivative of that
#'   log-density with respect to the predictor.
#' @param psi The term's parameters on the **parameter** scale, named as
#'   [term_params()].
#' @param ... Passed to methods.
#'
#' @return A list of two: `loglik`, a numeric vector of `n` contributions
#'   summing to the term's log-likelihood, and `jacobian`, an `n` by
#'   `length(psi)` matrix of their derivatives in the term's parameters,
#'   exact and propagated through the recursion.
#'
#' @seealso [regime()] and the marginal break-point terms for the two
#'   implementations, [term_filter()] for the other structural shape,
#'   [term_posterior()] for the smoothed states, [term_hessian()] for the
#'   observed information.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:40, y = c(rnorm(20), rnorm(20, 3)))
#' term <- term_build(regime(2, time = t), dd)
#' psi <- list(level1 = 0, gap2 = 3, alr1.1 = 2, alr2.1 = -2)
#'
#' out <- term_loglik(term, rep(0, 40), dd$y,
#'                    logdens = function(e, i) dnorm(dd$y[i], e, log = TRUE),
#'                    score = function(e, i) dd$y[i] - e,
#'                    psi = psi)
#' sum(out$loglik)
#' dim(out$jacobian)
#'
#' # The forward recursion is the sum over every state path. On eight
#' # observations that sum can be taken in full, and it agrees exactly.
#' d2 <- data.frame(t = 1:8, y = c(rnorm(4), rnorm(4, 3)))
#' t2 <- term_build(regime(2, time = t), d2)
#' fwd <- sum(term_loglik(t2, rep(0, 8), d2$y,
#'                        logdens = function(e, i) dnorm(d2$y[i], e, log = TRUE),
#'                        score = function(e, i) d2$y[i] - e,
#'                        psi = psi)$loglik)
#'
#' P <- as.matrix(parameters7::param_value(
#'   parameters7::transition_matrix(2), c(alr1.1 = 2, alr2.1 = -2)))
#' p0 <- Re(eigen(t(P))$vectors[, 1]); p0 <- p0 / sum(p0)
#' lev <- c(0, 3)                       # psi is on the parameter scale
#' paths <- as.matrix(expand.grid(rep(list(1:2), 8)))
#' tot <- sum(apply(paths, 1, function(s) {
#'   pr <- p0[s[1]]
#'   for (j in 2:8) pr <- pr * P[s[j - 1], s[j]]
#'   pr * prod(dnorm(d2$y, lev[s]))
#' }))
#' c(forward = fwd, all_256_paths = log(tot))
#'
#' @export
#' @aliases term_loglik.structural_term
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
#' The subclass of [structural_term()] holding a latent Markov chain of
#' regimes, each shifting the linear predictor by a level of its own.
#' [regime()] constructs it. Its contribution is a likelihood mixed over the
#' unobserved state path, so it implements [term_loglik()] and has no
#' predictor to report.
#'
#' @details
#' # The five properties of its own
#'
#' `k` is the number of regimes. `by` and `time` are the grouping and ordering
#' expressions as written, kept unevaluated; `NULL` means one group in row
#' order.
#'
#' `chain` is the [parameters7::transition_matrix()] object, whose free values
#' are the additive log-ratios of each row, so every row is a probability
#' vector at any coordinate. Its `free_names` are the tail of
#' [term_params()].
#'
#' `blueprint` is filled by [term_build()] and holds the row order within each
#' group and the observation count. The class overrides the branch's
#' `blueprint` property because a structural term carries no design block to
#' hang one on.
#'
#' @inheritParams model_term
#' @param k The number of regimes, an integer of at least 2.
#' @param by An optional grouping expression; each group runs its own
#'   recursion. `NULL` for one group.
#' @param time An optional ordering expression. `NULL` for row order.
#' @param chain A [parameters7::transition_matrix()] of side `k`.
#' @param blueprint A named list of the resolved ordering and grouping, empty
#'   until [term_build()] fills it.
#'
#' @return An S7 object of class `RegimeTerm`, inheriting from
#'   [structural_term()] and [model_term()], with the five properties above
#'   beside [model_term()]'s six.
#'
#' @seealso [regime()], the constructor; [term_loglik()] for what it computes;
#'   [term_posterior()] for the smoothed states; [gas()] for the other
#'   dynamic term.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:40, y = c(rnorm(20), rnorm(20, 3)))
#'
#' tm <- regime(2, time = t)
#' S7::S7_inherits(tm, RegimeTerm)
#' c(k = tm@k, chain = class(tm@chain)[1])
#'
#' # The build resolves the order and fills the blueprint.
#' b <- term_build(tm, dd)
#' names(b@blueprint)
#'
#' # It is on the structural branch, so there is no block.
#' try(term_matrix(b))
#'
#' @export
RegimeTerm <- S7::new_class(
  name = "RegimeTerm",
  parent = structural_term,
  properties = list(
    k = S7::class_integer,
    by = S7::class_any,
    time = S7::class_any,
    chain = S7::class_any
  )
)

#' Markov Regime Switching
#'
#' @description
#' A latent Markov chain of \eqn{K} regimes, each shifting the linear
#' predictor by a level of its own (Hamilton, 1989). The likelihood is the
#' mixture over the unobserved state path, evaluated by the forward recursion,
#' and it is built from whatever density the model carries: the term supplies
#' the chain and the levels, the distribution supplies everything else.
#'
#' @details
#' # The forward recursion
#'
#' Writing \eqn{f_j(t)} for the density of observation \eqn{t} at the
#' predictor shifted by the level of regime \eqn{j},
#'
#' \deqn{\tilde\alpha_t(j) = f_j(t) \sum_i \alpha_{t-1}(i) P_{ij},
#'   \qquad c_t = \sum_j \tilde\alpha_t(j), \qquad
#'   \alpha_t = \tilde\alpha_t / c_t,}
#'
#' started at the chain's stationary distribution, with the log-likelihood
#' \eqn{\sum_t \log c_t}. Those contributions are what [term_loglik()]
#' returns, one per observation, with their exact derivatives propagated
#' beside the state.
#'
#' Normalizing at every step is what keeps the recursion representable. The
#' unnormalized quantities are products of \eqn{t} densities, so they decay
#' geometrically and reach zero in double precision on a series of a few
#' hundred observations.
#'
#' # It is the complement of the other dynamic term
#'
#' [gas()] is driven by the score of the density and moves continuously; a
#' regime chain moves in jumps between a finite number of states. Both are
#' built from the density itself, so both apply to any family the model
#' carries, and both propagate their exact derivative beside the state.
#'
#' @section The parameters and their charts:
#' The levels are **ordered by construction**: `level1` is free on the
#' identity, and each of the others is the previous one plus a positive gap,
#' `gap2` ... `gapK`, each carried on a log link. Without an ordering the
#' regimes are exchangeable and the likelihood has \eqn{K!} identical maxima,
#' which is not a hard problem to fit but is one whose answer cannot be
#' reported.
#'
#' The transition matrix is [parameters7::transition_matrix()], whose free
#' values are the additive log-ratios of each row, named `alr1.1`, `alr2.1`
#' and so on, so every row is a probability vector at any coordinate. A chain
#' of \eqn{K} regimes therefore has \eqn{K(K-1)} of them, and
#' [term_params()] returns \eqn{K + K(K-1)} names in all: one level, \eqn{K-1}
#' gaps and the rest.
#'
#' The initial distribution is the chain's stationary one, which costs no
#' parameters; its derivative comes from the linear system it solves.
#'
#' [term_level_param()] answers `"level1"`, since a constant added to it
#' shifts every regime and is the direction an intercept in the same equation
#' also spans. The gaps are unaffected: what a constant cannot express is a
#' difference between regimes.
#'
#' @param k The number of regimes, a single whole number of at least 2.
#'   Anything else throws.
#' @param by An optional grouping variable, given as a bare expression. Each
#'   group runs its own recursion from the stationary distribution, so a panel
#'   of independent series is `by = id`. `NULL`, the default, is one group.
#' @param time An optional ordering variable, a bare expression. `NULL`, the
#'   default, takes the rows in the order they are given. Both `by` and `time`
#'   must evaluate to one non-missing value per row at build time.
#' @param label A single non-empty character string naming the term,
#'   `"regime"` by default.
#'
#' @return An unbuilt [RegimeTerm()]: a specification, whose `blueprint` is
#'   empty until [term_build()] resolves the ordering.
#'
#' @references
#' Hamilton, J. D. (1989). A new approach to the economic analysis of
#' nonstationary time series and the business cycle. *Econometrica*, 57(2),
#' 357--384.
#'
#' @seealso [gas()] for the continuous dynamic term, [term_loglik()] for the
#'   likelihood it computes, [term_posterior()] for the smoothed state
#'   probabilities, [term_hessian()] for the observed information.
#'
#' @examples
#' # One level, one gap, and two log-ratios of a two-state chain.
#' term_params(regime(2))
#' vapply(term_links(regime(2)), function(l) l@link_name, character(1))
#'
#' # K + K(K-1) parameters in all.
#' c(k3 = length(term_params(regime(3))), expected = 3 + 3 * 2)
#'
#' # The levels are ordered, so level1 is the one an intercept collides with.
#' term_level_param(regime(2))
#'
#' # Fitted through a model, the term is what carries the state; on its own
#' # it computes the mixed likelihood from a density the caller supplies.
#' set.seed(1)
#' dd <- data.frame(t = 1:40, y = c(rnorm(20), rnorm(20, 3)))
#' b <- term_build(regime(2, time = t), dd)
#' out <- term_loglik(b, rep(0, 40), dd$y,
#'                    logdens = function(e, i) dnorm(dd$y[i], e, log = TRUE),
#'                    score = function(e, i) dd$y[i] - e,
#'                    psi = list(level1 = 0, gap2 = 3, alr1.1 = 2, alr2.1 = -2))
#' sum(out$loglik)
#'
#' # And the smoothed probability of the second regime finds the change.
#' pp <- term_posterior(b, rep(0, 40), dd$y,
#'                      logdens = function(e, i) dnorm(dd$y[i], e, log = TRUE),
#'                      psi = list(level1 = 0, gap2 = 3, alr1.1 = 2,
#'                                 alr2.1 = -2))
#' round(pp[c(1, 19, 20, 21, 22, 40), 2], 3)
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

#' @title The Parameters of a Regime Term
#' @name term_params.RegimeTerm
#'
#' @description
#' `"level1"`, then `"gap2"` ... `"gapK"`, then the free names of the
#' transition matrix, `"alr1.1"` and the rest. A chain of \eqn{K} regimes has
#' \eqn{K + K(K-1)} parameters in all: one free level, \eqn{K-1} positive gaps
#' that order the others, and \eqn{K(K-1)} additive log-ratios, \eqn{K-1} per
#' row of the matrix.
#'
#' @details
#' The order matters: it is the order [term_links()], [term_start()],
#' [term_loglik()]'s Jacobian columns and the joint variance matrix are all
#' indexed by, and `psi` must supply the names whatever order they come in.
#'
#' The log-ratio names come from [parameters7::transition_matrix()]'s own
#' `free_names`, so the spelling is that package's.
#'
#' @param term A [RegimeTerm()], built or not: the names depend on `k` alone.
#' @param ... Unused.
#'
#' @return A character vector of length \eqn{K + K(K-1)}.
#'
#' @seealso [term_links()] for the chart each rides, [regime()] for what they
#'   mean.
#'
#' @examples
#' term_params(regime(2))
#' term_params(regime(3))
#' vapply(2:4, function(k) length(term_params(regime(k))), integer(1))
#'
#' @keywords internal
S7::method(term_params, RegimeTerm) <- function(term, ...) {
  c("level1",
    if (term@k > 1L) paste0("gap", seq.int(2L, term@k)),
    term@chain@free_names)
}

#' @title The Level of a Regime Term
#' @name term_level_param.RegimeTerm
#' @description
#' `"level1"`. The levels are ordered by construction, each the
#' previous plus a positive gap, so the first one shifts every regime and is
#' the direction an intercept in the same equation also spans. The gaps are
#' unaffected: what a constant cannot express is a difference between
#' regimes.
#' @param term A [RegimeTerm()].
#' @param ... Unused.
#' @return A single string.
#' @keywords internal
S7::method(term_level_param, RegimeTerm) <- function(term, ...) "level1"

#' @title The Charts of a Regime Term's Parameters
#' @name term_links.RegimeTerm
#'
#' @description
#' The **log** link on every gap and the identity on everything else. A gap
#' must be positive, that being what orders the levels; a level and an
#' additive log-ratio are already unconstrained.
#'
#' @details
#' The log-ratios need no chart because
#' [parameters7::transition_matrix()] has already provided one: they are the
#' additive log-ratios of each row, so any real values at all give a row of
#' probabilities summing to one.
#'
#' The consequence for a caller is that `psi` and `zeta` differ only in the
#' gaps. [term_loglik()] takes the parameter scale, where `gap2 = 3` is a gap
#' of three; [term_readable()] takes the unconstrained scale, where the same
#' number is a gap of \eqn{e^3}.
#'
#' @param term A [RegimeTerm()].
#' @param ... Unused.
#'
#' @return A named list of \pkg{linkfunctions7} links, one per entry of
#'   [term_params()].
#'
#' @seealso [term_params()], [regime()].
#'
#' @examples
#' vapply(term_links(regime(3)), function(l) l@link_name, character(1))
#'
#' @keywords internal
S7::method(term_links, RegimeTerm) <- function(term, ...) {
  nm <- term_params(term)
  stats::setNames(lapply(nm, function(p) {
    if (startsWith(p, "gap")) linkfunctions7::log_link()
    else linkfunctions7::identity_link()
  }), nm)
}

#' @title Build a Regime Term
#' @name term_build.RegimeTerm
#'
#' @description
#' Resolves the grouping and the ordering against the data and records them in
#' the blueprint. Nothing else is computed: the term has no design block, and
#' its chain and levels are already fixed by the constructor.
#'
#' @details
#' `by` and `time` are evaluated in `data` with [baseenv()] as the enclosure.
#' Each must give one value per row, and neither may be missing; both are
#' checked here, since this is the first point at which they meet data.
#'
#' The blueprint holds `order`, the row indices of each group sorted by time,
#' and `n`, the observation count. [term_loglik()] runs one forward recursion
#' per group in that order, each from the chain's stationary distribution, so
#' a panel of independent series is fitted by giving `by`.
#'
#' Without `time` the rows are taken in the order they appear. That is a real
#' choice worth making deliberately: the recursion is about order, so a data
#' frame that is not already sorted gives a different model.
#'
#' @param term A [RegimeTerm()].
#' @param data A data frame carrying whatever `by` and `time` name.
#' @param ... Unused.
#'
#' @return The term with `blueprint` filled. [term_is_built()] reads that
#'   property on this branch, so it is `TRUE` for the result.
#'
#' @seealso [regime()], [term_loglik()].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:40, id = rep(1:2, each = 20),
#'                  y = c(rnorm(20), rnorm(20, 3)))
#'
#' # One group in row order.
#' b <- term_build(regime(2, time = t), dd)
#' names(b@blueprint)
#' lengths(b@blueprint$order)
#'
#' # Two independent series, each with its own recursion.
#' b2 <- term_build(regime(2, by = id, time = t), dd)
#' lengths(b2@blueprint$order)
#'
#' # Both must give one value per row.
#' try(term_build(regime(2, time = c(1, 2)), dd))
#'
#' @keywords internal
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
#' Differentiating a second time gives
#' \eqn{d^2\delta\,(I - P) = d_i\delta\,d_jP + d_j\delta\,d_iP
#' + \delta\,d^2P} under \eqn{\sum_j d^2\delta_j = 0}, so the second
#' derivative comes from the same system again.
#'
#' @param P A row-stochastic matrix.
#' @param dP A list of derivative matrices of `P`.
#' @param d2P Optionally, the second derivatives: one list per row of
#'   `P` whose elements are the matrices of second derivatives of that
#'   row's entries, indexed as `dP` is.
#'
#' @return A list with `delta` and `ddelta`, the latter one row
#'   per element of `dP`, and, when `d2P` is given,
#'   `d2delta`, one square matrix per state.
#'
#' @keywords internal
regime_stationary <- function(P, dP, d2P = NULL) {
  k <- nrow(P)
  A <- diag(k) - P
  A[, k] <- 1
  delta <- as.numeric(solve(t(A), c(rep(0, k - 1L), 1)))
  ddelta <- lapply(dP, function(D) {
    rhs <- as.numeric(delta %*% D)
    rhs[k] <- 0
    as.numeric(solve(t(A), rhs))
  })
  out <- list(delta = delta, ddelta = ddelta)
  if (!is.null(d2P)) {
    np <- length(dP)
    d2 <- rep(list(matrix(0, np, np)), k)
    for (i in seq_len(np)) {
      for (j in seq_len(i)) {
        # the second derivative of P's entries in this pair of directions
        D2 <- matrix(0, k, k)
        for (a in seq_len(k)) for (b in seq_len(k)) {
          D2[a, b] <- d2P[[a]][[b]][i, j]
        }
        rhs <- as.numeric(ddelta[[i]] %*% dP[[j]] + ddelta[[j]] %*% dP[[i]] +
                            delta %*% D2)
        rhs[k] <- 0
        sol <- as.numeric(solve(t(A), rhs))
        for (s in seq_len(k)) {
          d2[[s]][i, j] <- sol[s]
          d2[[s]][j, i] <- sol[s]
        }
      }
    }
    out$d2delta <- d2
  }
  out
}

#' @title Log-Likelihood of a Regime Term
#' @name term_loglik.RegimeTerm
#' @description
#' Runs the forward recursion over each group in time order, normalizing
#' at every step, and returns the per-observation contributions with their
#' exact derivatives.
#' @param term A built `RegimeTerm`.
#' @param eta The static predictor.
#' @param y The response, reaching the recursion through the two callbacks.
#' @param logdens,score The log-density and its derivative in the predictor.
#' @param psi The parameters, named as [term_params()].
#' @param ... Unused.
#' @return A list with `loglik` and `jacobian`.
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

#' @title The Smoothed State Probabilities of a Latent Markov Term
#'
#' @description
#' The probability of each regime at each observation given the whole
#' series, which is everything a model layer needs to differentiate a
#' likelihood mixed over states.
#'
#' @details
#' [term_loglik()] returns the derivative of the mixed likelihood
#' in the term's own parameters, which is the piece estimating those needs.
#' Estimating the coefficients needs something else, and for this term that
#' something is not a second recursion carrying derivatives: it is one
#' quantity, by Fisher's identity. Writing \eqn{\gamma_t(k)} for the probability
#' returned here and \eqn{\theta_t(k)} for the parameters the model has at
#' observation \eqn{t} under regime \eqn{k},
#'
#' \deqn{\frac{\partial L}{\partial \eta_{q,t}}
#'   = \sum_k \gamma_t(k)\,
#'     \frac{\partial \ell(y_t; \theta_t(k))}{\partial \eta_q},}
#'
#' for every predictor the model carries, not only the one the regimes
#' shift. A caller therefore differentiates its own likelihood \eqn{K} times
#' vectorized and weights the results, and needs no callback per
#' observation: the regimes shift a predictor that is known before the
#' recursion starts, which is the property that made the forward pass
#' compilable and makes this cheap.
#'
#' The probabilities come from the forward pass this term already runs and a
#' backward pass beside it, both normalized at every step. Without the
#' normalization the quantities are products of \eqn{t} densities and reach
#' zero in double precision within a few hundred observations.
#'
#' @param term A built [RegimeTerm()].
#' @param eta The static predictor.
#' @param y The response.
#' @param logdens The log-density as a function of the predictor and the
#'   observation index, as [term_loglik()] takes it.
#' @param psi The term's parameters, named as [term_params()].
#' @param ... Passed to methods.
#'
#' @return A numeric matrix of `n` rows and \eqn{K} columns, every row
#'   summing to one, giving \eqn{P(S_t = k \mid y_1, \dots, y_n)}.
#'
#' @seealso [term_loglik()] for the likelihood the same recursion computes,
#'   [term_hessian()] for the observed information, [regime()] for the model.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:40, y = c(rnorm(20), rnorm(20, 3)))
#' term <- term_build(regime(2, time = t), dd)
#' psi <- list(level1 = 0, gap2 = 3, alr1.1 = 2, alr2.1 = -2)
#' g <- term_posterior(term, rep(0, 40), dd$y,
#'                     logdens = function(e, i) dnorm(dd$y[i], e, log = TRUE),
#'                     psi = psi)
#'
#' # Every row is a distribution over the regimes.
#' dim(g)
#' range(rowSums(g))
#'
#' # The data switch level at observation 20, and the smoothed
#' # probability of the second regime finds it.
#' round(g[c(1, 19, 20, 21, 22, 40), 2], 3)
#'
#' @export
#' @aliases term_posterior.structural_term
term_posterior <- S7::new_generic("term_posterior", "term",
  function(term, eta, y, logdens, psi, ...) S7::S7_dispatch())

S7::method(term_posterior, structural_term) <- function(term, eta, y, logdens,
                                                        psi, ...) {
  stop(sprintf("the term class '%s' does not implement term_posterior().",
               attr(S7::S7_class(term), "name")), call. = FALSE)
}

#' @title Smoothed State Probabilities of a Regime Term
#' @name term_posterior.RegimeTerm
#' @description
#' The forward pass of [term_loglik()] with a backward pass beside
#' it, both normalized, giving the probability of each regime at each
#' observation given the whole series.
#' @param term A built [RegimeTerm()].
#' @param eta The static predictor.
#' @param y The response.
#' @param logdens The log-density.
#' @param psi The term's parameters.
#' @param ... Unused.
#' @return A numeric matrix, one row per observation and one column per
#'   regime.
#' @keywords internal
S7::method(term_posterior, RegimeTerm) <- function(term, eta, y, logdens,
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
  n <- bp$n

  gaps <- if (k > 1L) v[paste0("gap", seq.int(2L, k))] else numeric(0)
  if (any(gaps <= 0)) {
    stop("every gap must be positive: the levels are ordered by construction.",
         call. = FALSE)
  }
  mu <- cumsum(c(v[["level1"]], gaps))
  P <- parameters7::param_value(term@chain, v[term@chain@free_names])
  st <- regime_stationary(P, list())
  delta <- st$delta

  idx <- seq_len(n)
  LF <- matrix(0, n, k)
  for (jj in seq_len(k)) {
    LF[, jj] <- as.numeric(logdens(eta + mu[jj], idx))
  }

  out <- matrix(0, n, k)
  for (rows in bp$order) {
    m <- length(rows)
    lf <- LF[rows, , drop = FALSE]
    mx <- apply(lf, 1L, max)
    W <- exp(lf - mx)

    # forward, normalized at every step
    al <- matrix(0, m, k)
    a <- delta
    for (t in seq_len(m)) {
      pred <- if (t == 1L) a else as.numeric(a %*% P)
      atil <- W[t, ] * pred
      a <- atil / sum(atil)
      al[t, ] <- a
    }
    # backward, normalized the same way: the scale cancels in the product
    be <- matrix(1, m, k)
    b <- rep(1, k)
    if (m > 1L) {
      for (t in seq.int(m - 1L, 1L)) {
        b <- as.numeric(P %*% (W[t + 1L, ] * b))
        s <- sum(b)
        if (s > 0) b <- b / s
        be[t, ] <- b
      }
    }
    g <- al * be
    out[rows, ] <- g / rowSums(g)
  }
  out
}

#' @title The Observed Hessian of a Likelihood Mixed Over States
#'
#' @description
#' The log-likelihood of a term that mixes over latent states, with its
#' exact gradient and its exact Hessian in the whole of a caller's unknown
#' vector: the coefficients of every equation together with the term's own
#' parameters.
#'
#' @details
#' [term_posterior()] gives the gradient by Fisher's identity, and
#' the matrix a caller can assemble from the same smoothed probabilities is
#' the complete-data information, the ordinary one averaged over the states.
#' That is the matrix an EM step inverts. It is not the observed information
#' of the mixture, which is smaller by the information the unobserved states
#' cost, and a standard error read off it is too small.
#'
#' Louis's identity expresses the difference as the conditional variance of
#' the complete-data score, which needs the pairwise smoothed probabilities
#' and a recursion carrying a second moment. The route taken here is
#' shorter. The scaled forward recursion computes the observed
#' log-likelihood exactly, as \eqn{\log L = \sum_t \log c_t} with \eqn{c_t}
#' the normalizing constant of one step, so differentiating that arithmetic
#' twice gives the observed Hessian with no identity involved. Carrying
#' \eqn{a_t(k)} together with its first and second derivatives in the
#' unknown vector \eqn{u},
#'
#' \deqn{\frac{\partial \log c_t}{\partial u} = \frac{\dot c_t}{c_t},
#'   \qquad
#'   \frac{\partial^2 \log c_t}{\partial u \partial u^\top}
#'     = \frac{\ddot c_t}{c_t}
#'       - \frac{\dot c_t}{c_t}\frac{\dot c_t^\top}{c_t},}
#'
#' and the state is renormalized by the quotient rule,
#'
#' \deqn{\ddot a = \bigl(\ddot{\tilde a} - \dot a \otimes \dot c
#'   - \dot c \otimes \dot a - a\,\ddot c\bigr)/c.}
#'
#' The emission enters through \eqn{\ddot w = w(gg^\top + G)}, with \eqn{g}
#' and \eqn{G} the score and the Hessian of the log-density in \eqn{u} at
#' that state, which is where the caller's derivatives are used. Groups are
#' independent series and their contributions add.
#'
#' The cost is \eqn{O(nK^2m^2)} in time and \eqn{O(Km^2)} in storage, and
#' the computation is meant to be run once, at a fitted point, never per
#' iteration.
#'
#' @param term A built [RegimeTerm()].
#' @param eta The static part of the predictor the regimes shift.
#' @param y The response.
#' @param logdens A function of a predictor value and a row index returning
#'   the log-density there, as [term_loglik()] takes it.
#' @param grad A function of the same two arguments returning a matrix with
#'   one row per observation and one column per distribution parameter, the
#'   derivative of the log-density in each predictor.
#' @param hess A function of the same two arguments returning an array of
#'   second derivatives in the predictors, one slice per observation.
#' @param psi The term's parameters, named as [term_params()].
#' @param seed A list with one matrix per distribution parameter, each with
#'   one row per observation and one column per unknown, giving the
#'   derivative of that predictor in the caller's unknowns with zeros in the
#'   columns the term's own parameters occupy.
#' @param cols The columns of the unknown vector the term's own parameters
#'   occupy, in the order of [term_params()].
#' @param level The index, among the distribution parameters, of the one the
#'   regimes shift.
#' @param weights Optional observation weights.
#' @param ... Passed to methods.
#'
#' @return A list of three: `loglik`, the `n` per-observation contributions;
#'   `gradient`, the derivative of their weighted sum in the caller's whole
#'   unknown vector; and `hessian`, the observed Hessian of the mixture in
#'   that same vector, an `m` by `m` matrix with `m` its length. The Hessian
#'   is the mixture's own, smaller in the Loewner order than the
#'   complete-data information by the information the unobserved states cost.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:30, y = c(rnorm(15), rnorm(15, 3)))
#' term <- term_build(regime(2, time = t), dd)
#' # a gaussian mean of unit variance: one predictor, one unknown besides
#' # the term's own four
#' m <- 5L
#' out <- term_hessian(
#'   term, eta = rep(0, 30), y = dd$y,
#'   logdens = function(e, i) dnorm(dd$y[i], e, log = TRUE),
#'   grad = function(e, i) matrix(dd$y[i] - e, ncol = 1L),
#'   hess = function(e, i) array(-1, c(length(i), 1L, 1L)),
#'   psi = list(level1 = 0, gap2 = 3, alr1.1 = 2, alr2.1 = -2),
#'   seed = list(matrix(c(rep(1, 30), rep(0, 30 * 4)), 30, m)),
#'   cols = 2:5, level = 1L)
#' dim(out$hessian)
#'
#' @seealso [term_posterior()], [term_loglik()]
#' @export
#' @aliases term_hessian.structural_term
term_hessian <- S7::new_generic("term_hessian", "term",
  function(term, eta, y, logdens, grad, hess, psi, seed, cols, level,
           weights = NULL, ...) S7::S7_dispatch())

S7::method(term_hessian, structural_term) <- function(term, eta, y, logdens,
                                                      grad, hess, psi, seed,
                                                      cols, level,
                                                      weights = NULL, ...) {
  stop(sprintf("the term class '%s' does not implement term_hessian().",
               attr(S7::S7_class(term), "name")), call. = FALSE)
}

# The chart of a regime's own parameters, on the unconstrained scale the
# model layer optimizes on: the levels are cumulative sums of positive
# gaps, so a gap under its log link contributes itself to both the first
# and the second derivative, and the transition parameters are identity
# linked, so the chain's own derivative arrays are already in zeta.
.regime_chart <- function(v, k, nm, chain) {
  np <- length(nm)
  gaps <- if (k > 1L) v[paste0("gap", seq.int(2L, k))] else numeric(0)
  if (any(gaps <= 0)) {
    stop("every gap must be positive: the levels are ordered by construction.",
         call. = FALSE)
  }
  mu <- cumsum(c(v[["level1"]], gaps))

  dmu <- matrix(0, k, np)
  dmu[, 1L] <- 1
  d2mu <- rep(list(matrix(0, np, np)), k)
  for (j in seq_len(k)) {
    if (j > 1L) {
      slot <- 1L + seq_len(j - 1L)
      dmu[j, slot] <- gaps[seq_len(j - 1L)]
      d2mu[[j]][cbind(slot, slot)] <- gaps[seq_len(j - 1L)]
    }
  }

  i_tr <- k + seq_len(k * (k - 1L))
  fn <- chain@free_names
  eta_tr <- v[fn]
  P <- parameters7::param_value(chain, eta_tr)
  d1 <- parameters7::param_d1(chain, eta_tr)
  d2 <- parameters7::param_d2(chain, eta_tr)

  # dPz[j, l, i] = dP[j, l]/dzeta_i, and D2P[[j]][[l]] the np x np second
  # derivative of that entry; both are zero outside the chain's own slots
  dPz <- array(0, c(k, k, np))
  for (i in seq_along(fn)) dPz[, , i_tr[i]] <- d1[[fn[i]]]
  D2P <- lapply(seq_len(k), function(j) rep(list(matrix(0, np, np)), k))
  for (i in seq_along(fn)) {
    for (jj in seq_along(fn)) {
      key <- if (paste0(fn[i], ":", fn[jj]) %in% names(d2)) {
        paste0(fn[i], ":", fn[jj])
      } else {
        paste0(fn[jj], ":", fn[i])
      }
      M <- d2[[key]]
      if (is.null(M)) next
      for (a in seq_len(k)) for (b in seq_len(k)) {
        D2P[[a]][[b]][i_tr[i], i_tr[jj]] <- M[a, b]
      }
    }
  }

  st <- regime_stationary(P, lapply(seq_len(np), function(i) dPz[, , i]),
                          D2P)
  list(mu = mu, dmu = dmu, d2mu = d2mu, P = P, dPz = dPz, D2P = D2P,
       delta = st$delta, ddelta = st$ddelta, d2delta = st$d2delta, np = np)
}

#' @title Observed Hessian of a Regime Term
#' @name term_hessian.RegimeTerm
#' @description
#' Forward-mode second-order propagation through the scaled forward
#' recursion, giving the exact Hessian of the mixed log-likelihood.
#' @param term A built [RegimeTerm()].
#' @param eta The static predictor.
#' @param y The response.
#' @param logdens,grad,hess The log-density and its first two derivatives in
#'   the predictors.
#' @param psi The term's parameters.
#' @param seed The derivative of each predictor in the caller's unknowns.
#' @param cols The columns the term's own parameters occupy.
#' @param level The distribution parameter the regimes shift.
#' @param weights Observation weights.
#' @param ... Unused.
#' @return A list with `loglik`, `gradient` and `hessian`.
#' @keywords internal
S7::method(term_hessian, RegimeTerm) <- function(term, eta, y, logdens, grad,
                                                 hess, psi, seed, cols, level,
                                                 weights = NULL, ...) {
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
  n <- bp$n
  if (!is.list(seed) || !length(seed)) {
    stop("'seed' must be a list with one matrix per distribution parameter.",
         call. = FALSE)
  }
  seed <- lapply(seed, as.matrix)
  m <- ncol(seed[[1L]])
  npar <- length(seed)
  if (any(vapply(seed, nrow, 1L) != n) ||
      any(vapply(seed, ncol, 1L) != m)) {
    stop(sprintf("every element of 'seed' must be %d by %d.", n, m),
         call. = FALSE)
  }
  cols <- as.integer(cols)
  level <- as.integer(level)
  if (length(cols) != length(nm) || any(cols < 1L) || any(cols > m)) {
    stop(sprintf("'cols' must give the %d columns the term's parameters occupy.",
                 length(nm)), call. = FALSE)
  }
  if (length(level) != 1L || level < 1L || level > npar) {
    stop("'level' must index one of the distribution parameters.",
         call. = FALSE)
  }
  w <- if (is.null(weights)) rep(1, n) else rep_len(as.numeric(weights), n)

  ch <- .regime_chart(v, k, nm, term@chain)
  np <- ch$np

  # every regime shifts one predictor by a level of its own, so the
  # density and its derivatives under every regime are known before the
  # recursion starts, as they are for the first-order pass
  idx <- seq_len(n)
  LF <- matrix(0, n, k)
  GP <- vector("list", k)
  HP <- vector("list", k)
  GM <- vector("list", k)
  for (j in seq_len(k)) {
    e <- eta + ch$mu[j]
    LF[, j] <- as.numeric(logdens(e, idx))
    g <- as.matrix(grad(e, idx))
    h <- hess(e, idx)
    if (nrow(g) != n || ncol(g) != npar) {
      stop(sprintf("'grad' must return an %d by %d matrix.", n, npar),
           call. = FALSE)
    }
    GP[[j]] <- g
    HP[[j]] <- array(as.numeric(h), c(n, npar, npar))
    # the emission score in the caller's unknowns, vectorized over rows
    gm <- matrix(0, n, m)
    for (q in seq_len(npar)) gm <- gm + g[, q] * seed[[q]]
    gm[, cols] <- gm[, cols] + outer(g[, level], ch$dmu[j, ])
    GM[[j]] <- gm
  }

  P <- ch$P
  loglik <- numeric(n)
  gradient <- numeric(m)
  hessian <- matrix(0, m, m)

  # the rows of the per-observation predictor derivative, shared by every
  # regime except in the columns the term's own parameters occupy
  Abase <- matrix(0, npar, m)

  for (rows in bp$order) {
    a <- ch$delta
    da <- matrix(0, k, m)
    da[, cols] <- t(do.call(rbind, ch$ddelta))
    d2a <- vector("list", k)
    for (j in seq_len(k)) {
      M <- matrix(0, m, m)
      M[cols, cols] <- ch$d2delta[[j]]
      d2a[[j]] <- M
    }

    for (tt in seq_along(rows)) {
      row <- rows[tt]
      lf <- LF[row, ]
      mx <- max(lf)
      wv <- exp(lf - mx)

      for (q in seq_len(npar)) Abase[q, ] <- seed[[q]][row, ]

      if (tt == 1L) {
        pred <- a
        dpred <- da
        d2pred <- d2a
      } else {
        pred <- as.numeric(a %*% P)
        dpred <- crossprod(P, da)
        d2pred <- vector("list", k)
        for (l in seq_len(k)) {
          Dl <- matrix(ch$dPz[, l, ], k, np)
          # the chain's own contribution: supported on the term's columns
          dpred[l, cols] <- dpred[l, cols] + as.numeric(a %*% Dl)
          S1 <- matrix(0, m, m)
          for (j in seq_len(k)) {
            if (P[j, l] != 0) S1 <- S1 + P[j, l] * d2a[[j]]
          }
          V <- crossprod(da, Dl)   # m x np
          S1[, cols] <- S1[, cols] + V
          S1[cols, ] <- S1[cols, ] + t(V)
          S2 <- matrix(0, np, np)
          for (j in seq_len(k)) S2 <- S2 + a[j] * ch$D2P[[j]][[l]]
          S1[cols, cols] <- S1[cols, cols] + S2
          d2pred[[l]] <- S1
        }
      }

      atil <- wv * pred
      datil <- matrix(0, k, m)
      d2atil <- vector("list", k)
      for (l in seq_len(k)) {
        gl <- GM[[l]][row, ]
        datil[l, ] <- wv[l] * (gl * pred[l] + dpred[l, ])
        Ak <- Abase
        Ak[level, cols] <- Ak[level, cols] + ch$dmu[l, ]
        Hk <- matrix(HP[[l]][row, , ], npar, npar)
        G <- crossprod(Ak, Hk %*% Ak)
        G[cols, cols] <- G[cols, cols] + GP[[l]][row, level] * ch$d2mu[[l]]
        dp <- dpred[l, ]
        M <- pred[l] * (outer(gl, gl) + G) + d2pred[[l]]
        M <- M + outer(gl, dp) + outer(dp, gl)
        d2atil[[l]] <- wv[l] * M
      }

      ct <- sum(atil)
      dct <- colSums(datil)
      d2ct <- matrix(0, m, m)
      for (l in seq_len(k)) d2ct <- d2ct + d2atil[[l]]

      gt <- dct / ct
      loglik[row] <- log(ct) + mx
      gradient <- gradient + w[row] * gt
      hessian <- hessian + w[row] * (d2ct / ct - outer(gt, gt))

      a <- atil / ct
      da <- (datil - outer(atil, dct) / ct) / ct
      for (l in seq_len(k)) {
        d2a[[l]] <- (d2atil[[l]] - outer(da[l, ], dct) - outer(dct, da[l, ]) -
                       a[l] * d2ct) / ct
      }
    }
  }

  hessian <- (hessian + t(hessian)) / 2
  list(loglik = loglik, gradient = gradient, hessian = hessian)
}

#' @title Print a Regime Term
#' @name print.RegimeTerm
#'
#' @description
#' Prints two lines: the label and the number of regimes, followed by the
#' term's parameters. A built term reports how many groups the recursion runs
#' over; a specification says so instead.
#'
#' @details
#' The two forms are
#'
#' ```
#' <RegimeTerm> 'regime': 2 regimes (specification)
#'   parameters: level1, gap2, alr1.1, alr2.1
#'
#' <RegimeTerm> 'regime': 2 regimes; 1 group(s)
#'   parameters: level1, gap2, alr1.1, alr2.1
#' ```
#'
#' The group count is `length(blueprint$order)`, which is 1 unless the term
#' was given a `by`. The parameter list is [term_params()], the same before
#' and after a build, `k` alone determining it.
#'
#' Note that a built term is described by its count of groups rather than by
#' the word "built", the count being the more useful line.
#'
#' @param x A [RegimeTerm()], built or not.
#' @param ... Unused, and accepted so that the signature matches [print()]'s.
#'
#' @return `x`, invisibly. Called for the two lines it writes.
#'
#' @seealso [regime()], [term_params()].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:40, id = rep(1:2, each = 20),
#'                  y = c(rnorm(20), rnorm(20, 3)))
#'
#' # A specification, and the same term built over two groups.
#' regime(2)
#' term_build(regime(2, by = id, time = t), dd)
#'
#' # Three regimes carry three levels and six log-ratios.
#' regime(3)
#'
#' @keywords internal
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


#' @name term_simulate.RegimeTerm
#'
#' @title Drawing a Regime Path
#'
#' @description
#' A path of the latent chain, and the predictor each observation gets from
#' the level of the state it lands in.
#'
#' @details
#' The chain is started from its stationary distribution. The model's
#' likelihood is written with that initial law, so starting from a fixed
#' first state or from a uniform draw would simulate a different model from
#' the one a fit reads back. The transition matrix and the
#' stationary law come from the term's own
#' [parameters7::transition_matrix()], so nothing about the chart
#' is restated here.
#'
#' The response is not drawn, the levels not reading it, so `y` comes back
#' `NULL` and the caller draws at the returned predictor.
#'
#' @param term A built regime term.
#' @param psi The term's parameters, on the parameter scale.
#' @param eta The static part of the predictor.
#' @param draw Ignored: the levels do not read the response.
#' @param ... Ignored.
#'
#' @return A list with `eta`, `y` (`NULL`) and
#'   `latent`, the state of each observation.
#'
#' @seealso [term_simulate()], [regime()]
#'
#' @keywords internal
S7::method(term_simulate, RegimeTerm) <- function(term, psi, eta, draw, ...) {
  bp <- term@blueprint
  if (!length(bp)) {
    stop("the term has not been built; call term_build(term, data) first.",
         call. = FALSE)
  }
  n <- length(eta)
  k <- term@k
  v <- unlist(psi[term_params(term)])
  mu <- term_levels(term, psi)
  P <- parameters7::param_value(term@chain, v[term@chain@free_names])
  A <- diag(k) - P
  A[, k] <- 1
  delta <- as.numeric(solve(t(A), c(rep(0, k - 1L), 1)))
  st <- integer(n)
  if (n) {
    st[[1L]] <- sample.int(k, 1L, prob = pmax(delta, 0))
    for (i in seq_len(n - 1L)) {
      st[[i + 1L]] <- sample.int(k, 1L, prob = pmax(P[st[[i]], ], 0))
    }
  }
  list(eta = as.numeric(eta) + mu[st], y = NULL, latent = st)
}
