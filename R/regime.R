#' @include term_classes.R generics.R structural.R
NULL

#' @title Log-Likelihood Contributions of a Structural Term
#'
#' @description
#' Returns the per-observation log-likelihood contributions that a structural
#' term produces, with their derivatives in the term's own parameters. This is
#' the second shape that the structural branch takes, beside [term_filter()]: a
#' term that shifts the predictor implements the filter, and a term that
#' rewrites the likelihood itself implements this generic, because its
#' contribution is not a predictor at all.
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
#' The term holds the chain and the levels and contains nothing about the
#' family. The model's own density enters through two callbacks: `logdens`
#' returns the log-density at a given predictor, one value per observation, and
#' `score` returns its derivative in that predictor. For [regime()], each is
#' called once per regime, with the predictor shifted by that regime's level
#' and the whole index vector of the observations, so a term of \eqn{K}
#' regimes makes \eqn{K} calls to each.
#'
#' The callbacks of [term_filter()] depend on the state that the recursion has
#' reached. The callbacks of a regime term do not depend on it, because a regime
#' shifts a predictor that is known in advance. The density and the score of
#' every observation under every regime can therefore be computed once,
#' vectorized, before the recursion starts.
#'
#' The method on [structural_term()] signals an error that names the class. A
#' term of the filter shape implements [term_filter()] instead, and a fitting
#' layer distinguishes the two shapes by which of the two generics has a method
#' for the term.
#'
#' @param term A built structural term.
#' @param eta The static part of the linear predictor, one value per
#'   observation.
#' @param y The response. It reaches the recursion through the callbacks; the
#'   argument is passed for methods that need it directly.
#' @param logdens A function `(e, i)` returning the log-density at predictor
#'   `e` of each observation in the index vector `i`, one value per element of
#'   `i`.
#' @param score A function `(e, i)` returning the derivative of that
#'   log-density with respect to the predictor, one value per element of `i`.
#' @param psi The term's parameters on the parameter scale, named as
#'   [term_params()].
#' @param ... Passed to methods.
#'
#' @return A list of two: `loglik`, a numeric vector of `n` contributions
#'   summing to the term's log-likelihood, and `jacobian`, an `n` by
#'   `length(psi)` matrix of their derivatives in the term's parameters,
#'   computed exactly by propagation through the recursion. The columns follow
#'   [term_params()], and the derivative in a gap is taken in the gap itself,
#'   not in its logarithm.
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
#' `RegimeTerm` is the subclass of [structural_term()] that holds a latent
#' Markov chain of regimes, and each regime shifts the linear predictor by a
#' level of its own. [regime()] constructs it. Its contribution is a likelihood
#' mixed over the unobserved state path, so it implements [term_loglik()] and
#' has no predictor to report.
#'
#' @details
#' # The five properties of its own
#'
#' `k` is the number of regimes. `by` and `time` are the grouping and ordering
#' expressions as written, kept unevaluated. A `by` of `NULL` means one group,
#' and a `time` of `NULL` means that the rows are taken in the order in which
#' they are given.
#'
#' `chain` is the [parameters7::transition_matrix()] object, whose free values
#' are the additive log-ratios of each row, so every row is a probability
#' vector at any coordinate. Its `free_names` are the tail of
#' [term_params()].
#'
#' `n_start` is the number of starting points that a fitting layer tries, a
#' whole number of at least 1.
#'
#' The property `build_info` is inherited from [structural_term()], which
#' declares it. [term_build()] fills it with the row order within each group
#' and the observation count.
#'
#' @inheritParams model_term
#' @param k The number of regimes, an integer of at least 2.
#' @param by An optional grouping expression; each group runs its own
#'   recursion. `NULL` means one group.
#' @param time An optional ordering expression. `NULL` means row order.
#' @param chain A [parameters7::transition_matrix()] of side `k`.
#' @param build_info A named list of the resolved ordering and grouping, empty
#'   until [term_build()] fills it.
#' @param n_start The number of starting points that a fitting layer tries, a
#'   single whole number of at least 1, as set by [regime()].
#'
#' @return An S7 object of class `RegimeTerm`, inheriting from
#'   [structural_term()] and [model_term()], with the five properties above in
#'   addition to the eight that [structural_term()] supplies (the seven of
#'   [model_term()] and `build_info`).
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
#' # The build resolves the order and fills the build_info.
#' b <- term_build(tm, dd)
#' names(b@build_info)
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
    chain = S7::class_any,
    n_start = S7::class_integer
  )
)

#' Markov Regime Switching
#'
#' @description
#' Constructs a latent Markov chain of \eqn{K} regimes, and each regime shifts
#' the linear predictor by a level of its own (Hamilton, 1989). The likelihood
#' is the mixture over the unobserved state path, evaluated by the forward
#' recursion, and it is built from any density that the model
#' carries: the term supplies the chain and the levels, the distribution
#' supplies everything else.
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
#' where the predictive distribution of the first observation of each group is
#' the chain's stationary distribution. The log-likelihood is
#' \eqn{\sum_t \log c_t}. [term_loglik()] returns those contributions, one per
#' observation, with their exact derivatives propagated beside the state.
#'
#' Normalizing at every step keeps the recursion representable. The
#' unnormalized quantities are products of \eqn{t} densities, so they decay
#' geometrically and reach zero in double precision on a long series.
#'
#' # Relation to gas()
#'
#' [gas()] is driven by the score of the density and moves continuously, while
#' a regime chain moves in jumps between a finite number of states. Both are
#' built from the density itself, so both apply to any family that the model
#' carries, and both propagate their exact derivative beside the state.
#'
#' @section The parameters and their charts:
#' The levels are ordered by construction: `level1` is free on the
#' identity, and each of the others is the previous one plus a positive gap,
#' `gap2` ... `gapK`, and every gap is carried on a log link. Without an
#' ordering the regimes are exchangeable and the likelihood has \eqn{K!}
#' identical maxima, so the estimates of the individual regimes cannot be
#' reported without a labeling.
#'
#' The transition matrix is [parameters7::transition_matrix()], whose free
#' values are the additive log-ratios of each row, named `alr1.1`, `alr2.1`
#' and so on, so every row is a probability vector at any coordinate. A chain
#' of \eqn{K} regimes therefore has \eqn{K(K-1)} of them, and
#' [term_params()] returns \eqn{K + K(K-1)} names in all: one level,
#' \eqn{K-1} gaps and \eqn{K(K-1)} log-ratios.
#'
#' The initial distribution is the chain's stationary one, which costs no
#' parameters. Its derivative comes from the same linear system that gives the
#' distribution itself.
#'
#' [term_level_param()] returns `"level1"`, since a constant added to it
#' shifts every regime and is the direction that an intercept in the same
#' equation also spans. A constant leaves the gaps unchanged, because it cannot
#' express a difference between regimes.
#'
#' @param k The number of regimes, a single finite whole number of at least 2.
#'   Any other value signals an error.
#' @param by An optional grouping variable, given as a bare expression. Each
#'   group runs its own recursion from the stationary distribution, so a panel
#'   of independent series is `by = id`. `NULL`, the default, is one group.
#' @param time An optional ordering variable, a bare expression. `NULL`, the
#'   default, takes the rows in the order in which they are given. Both `by` and
#'   `time` must evaluate to one non-missing value per row at build time. Within
#'   a group, `time` only orders the rows: a gap between two times counts as one
#'   step of the chain.
#' @param label A single non-empty character string naming the term,
#'   `"regime"` by default.
#' @param n_start The number of starting points that a fitting layer tries, a
#'   single finite whole number of at least 1. The likelihood of a
#'   regime model has several maxima, and a fit started at one point reaches the
#'   maximum of the basin that contains it. The first start is [term_start()]'s,
#'   which places the levels at the quantiles of the response when the fitting
#'   layer supplies it and at zero otherwise. Each further start displaces the
#'   first as [term_starts()] describes, and the fit with the best criterion is
#'   kept. `1`, the default, fits from the first start alone. The model is
#'   fitted once from each start.
#'
#' @return An unbuilt [RegimeTerm()]: a specification, whose `build_info` is
#'   empty until [term_build()] resolves the ordering.
#'
#' @references
#' Hamilton, J. D. (1989). A new approach to the economic analysis of
#' nonstationary time series and the business cycle. *Econometrica*, 57(2),
#' 357--384.
#'
#' @seealso [gas()] for the continuous dynamic term, [term_loglik()] for the
#'   likelihood that it computes, [term_posterior()] for the smoothed state
#'   probabilities, [term_hessian()] for the observed information.
#'
#' @examples
#' # The term has one level, one gap, and two log-ratios of a two-state chain.
#' term_params(regime(2))
#' vapply(term_links(regime(2)), function(l) l@link_name, character(1))
#'
#' # A chain of K regimes has K + K(K-1) parameters in all.
#' c(k3 = length(term_params(regime(3))), expected = 3 + 3 * 2)
#'
#' # The levels are ordered, so level1 is the one that collides with an
#' # intercept.
#' term_level_param(regime(2))
#'
#' # Fitted through a model, the term carries the state; on its own it
#' # computes the mixed likelihood from a density that the caller supplies.
#' set.seed(1)
#' dd <- data.frame(t = 1:40, y = c(rnorm(20), rnorm(20, 3)))
#' b <- term_build(regime(2, time = t), dd)
#' out <- term_loglik(b, rep(0, 40), dd$y,
#'                    logdens = function(e, i) dnorm(dd$y[i], e, log = TRUE),
#'                    score = function(e, i) dd$y[i] - e,
#'                    psi = list(level1 = 0, gap2 = 3, alr1.1 = 2, alr2.1 = -2))
#' sum(out$loglik)
#'
#' # The smoothed probability of the second regime marks the change.
#' pp <- term_posterior(b, rep(0, 40), dd$y,
#'                      logdens = function(e, i) dnorm(dd$y[i], e, log = TRUE),
#'                      psi = list(level1 = 0, gap2 = 3, alr1.1 = 2,
#'                                 alr2.1 = -2))
#' round(pp[c(1, 19, 20, 21, 22, 40), 2], 3)
#'
#'
#' # The model is fitted. The data are simulated from a known truth, so the
#' # estimates below can be read against it.
#' if (requireNamespace("statmodels7", quietly = TRUE)) {
#'   set.seed(4)
#'   st <- integer(200)
#'   st[1] <- 1L
#'   for (i in 2:200) {
#'     st[i] <- if (runif(1) < 0.05) 3L - st[i - 1L] else st[i - 1L]
#'   }
#'   fd <- data.frame(t = 1:200, y = c(-2, 2)[st] + rnorm(200, sd = 0.7))
#'   ft <- statmodels7::statmod(y ~ 0 + regime(k = 2, time = t),
#'                              distributions7::gaussian1_distrib(), fd)
#'   # The truth has levels at -2 and 2, so the first level is -2 and the
#'   # gap is 4.
#'   round(coef(ft)$mu[1:2], 2)
#' }
#' @export
regime <- function(k = 2, by = NULL, time = NULL, label = "regime",
                   n_start = 1) {
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 2 ||
      k != round(k)) {
    stop("'k' must be a single integer of at least 2.", call. = FALSE)
  }
  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !nzchar(label)) {
    stop("'label' must be a single non-empty character string.",
         call. = FALSE)
  }
  if (!is.numeric(n_start) || length(n_start) != 1L || !is.finite(n_start) ||
      n_start < 1 || n_start != round(n_start)) {
    stop("'n_start' must be a single whole number of at least 1.",
         call. = FALSE)
  }
  k <- as.integer(k)
  RegimeTerm(label = label, k = k,
             by = substitute(by), time = substitute(time),
             chain = parameters7::transition_matrix(k),
             n_start = as.integer(n_start),
             build_info = list())
}

#' @title The Parameters of a Regime Term
#' @name term_params.RegimeTerm
#'
#' @description
#' The parameters are `"level1"`, then `"gap2"` ... `"gapK"`, then the free
#' names of the transition matrix, `"alr1.1"` and the rest. A chain of \eqn{K}
#' regimes has \eqn{K + K(K-1)} parameters in all: one free level, \eqn{K-1}
#' positive gaps that order the others, and \eqn{K(K-1)} additive log-ratios,
#' \eqn{K-1} per row of the matrix.
#'
#' @details
#' The order is the one that [term_links()], [term_start()], the columns of
#' [term_loglik()]'s Jacobian and the joint variance matrix all follow. The
#' argument `psi` of the methods is matched by name, so it may list the
#' parameters in any order.
#'
#' The log-ratio names come from [parameters7::transition_matrix()]'s own
#' `free_names`, so the spelling is that package's.
#'
#' @param term A [RegimeTerm()], built or not: the names depend on `k` alone.
#' @param ... Unused.
#'
#' @return A character vector of length \eqn{K + K(K-1)}.
#'
#' @seealso [term_links()] for the chart on which each parameter is carried,
#'   [regime()] for what they mean.
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
#' Returns `"level1"`. The levels are ordered by construction: each is the
#' previous one plus a positive gap, so the first one shifts every regime and
#' is the direction that an intercept in the same equation also spans. A
#' constant leaves the gaps unchanged, because it cannot express a difference
#' between regimes.
#' @param term A [RegimeTerm()].
#' @param ... Unused.
#' @return A single string.
#' @keywords internal
S7::method(term_level_param, RegimeTerm) <- function(term, ...) "level1"

#' @title What a Fitted Regime Term Reports
#' @name term_readable.RegimeTerm
#' @description
#' Returns the first level, the gaps between consecutive levels and the
#' transition probabilities \eqn{p_{ij} = P(S_t = j \mid S_{t-1} = i)}, with
#' the Jacobian from the term's own parameters.
#' @details
#' The levels and the gaps are reported through their own links, as the base
#' method does. The chain is estimated on the additive log-ratios of each row
#' of the transition matrix. These are coordinates of the chain and not
#' quantities that have a direct interpretation, so every entry of the matrix
#' is reported in their place, \eqn{k^2} of them, row by row and named `p1.1`,
#' `p1.2`, and so on. Entry \eqn{p_{ij}} depends only on the log-ratios of
#' row \eqn{i}, and its row of the Jacobian is the derivative of the value of
#' [parameters7::transition_matrix()]. Each interval is built on the logit
#' scale and mapped back, so it stays inside \eqn{(0, 1)}.
#' @param term A [RegimeTerm()].
#' @param zeta The parameters on the unconstrained scale.
#' @param ... Unused.
#' @return A list, as [term_readable()] documents.
#' @keywords internal
S7::method(term_readable, RegimeTerm) <- function(term, zeta, ...) {
  nm <- term_params(term)
  out <- S7::method(term_readable, model_term)(term, zeta)
  fr <- term@chain@free_names
  keep <- !(nm %in% fr)
  k <- term@k
  z <- unlist(zeta[nm])
  eta <- z[fr]
  P <- as.matrix(parameters7::param_value(term@chain, eta))
  dP <- parameters7::param_d1(term@chain, eta)
  ij <- expand.grid(j = seq_len(k), i = seq_len(k))[, c("i", "j")]
  pn <- sprintf("p%d.%d", ij$i, ij$j)
  Jp <- matrix(0, nrow(ij), length(nm), dimnames = list(pn, nm))
  for (f in fr) Jp[, f] <- dP[[f]][cbind(ij$i, ij$j)]
  J <- rbind(out$jacobian[keep, , drop = FALSE], Jp)
  list(name = c(out$name[keep], pn),
       value = c(out$value[keep], P[cbind(ij$i, ij$j)]),
       jacobian = J,
       scale = c(out$scale[keep],
                 stats::setNames(rep(list(linkfunctions7::logit_link()),
                                     nrow(ij)), pn)))
}

#' @title The Charts of a Regime Term's Parameters
#' @name term_links.RegimeTerm
#'
#' @description
#' Returns the log link for every gap and the identity link for every other
#' parameter. A gap must be positive, and positivity orders the levels; a
#' level and an additive log-ratio are already unconstrained.
#'
#' @details
#' The log-ratios need no chart because
#' [parameters7::transition_matrix()] has already provided one: they are the
#' additive log-ratios of each row, so any real values at all give a row of
#' probabilities summing to one.
#'
#' Therefore `psi` and `zeta` differ only in the gaps. [term_loglik()] takes
#' the parameter scale, where `gap2 = 3` is a gap of three; [term_readable()]
#' takes the unconstrained scale, where the same number is a gap of
#' \eqn{e^3}.
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
  # one object of each kind, shared: building an S7 link validates it, and
  # this is read at every evaluation of a fit (38 per cent of a
  # three-regime fit on geyser when built once per parameter)
  lg <- linkfunctions7::log_link()
  id <- linkfunctions7::identity_link()
  stats::setNames(lapply(nm, function(p) if (startsWith(p, "gap")) lg else id),
                  nm)
}

#' @title Where a Regime Term's Parameters Start
#' @name term_start.RegimeTerm
#' @description
#' Returns the starting values: the levels are placed at the quantiles of the
#' response at \eqn{(2j-1)/(2k)}, \eqn{j = 1, \dots, k}, and every row of the
#' transition matrix is uniform.
#' @details
#' Without the response the start is the zero vector on the unconstrained
#' scale: the first level is zero, every gap is one and every row of the
#' transition matrix is uniform, so the levels are \eqn{0, 1, \dots, k - 1}
#' whatever the scale of the data. The quantiles spread the regimes over the
#' data, as [seg_start()] spreads break-points over a covariate. `target` is
#' the response on the scale of the predictor, which the fitting layer
#' supplies where the family reads the parameter directly. Where it is
#' absent, where fewer than \eqn{k} of its values are finite, or where two
#' quantiles coincide, the conventional start above is returned.
#' @param term A [RegimeTerm()].
#' @param ... Unused.
#' @param target The response on the scale of the predictor, or `NULL`.
#' @return A named numeric vector on the unconstrained scale.
#' @keywords internal
S7::method(term_start, RegimeTerm) <- function(term, ..., target = NULL) {
  nm <- term_params(term)
  out <- stats::setNames(numeric(length(nm)), nm)
  k <- term@k
  y <- if (is.null(target)) numeric(0) else target[is.finite(target)]
  if (k < 2L || length(y) < k) return(out)
  q <- stats::quantile(y, (2 * seq_len(k) - 1) / (2 * k), names = FALSE)
  gaps <- diff(q)
  if (any(!(gaps > 0))) return(out)
  out[["level1"]] <- q[[1L]]
  out[paste0("gap", seq.int(2L, k))] <- log(gaps)
  out
}

#' @title The Starting Points of a Regime Term
#' @name term_starts.RegimeTerm
#' @description
#' Returns `n_start` starting points. The first is [term_start()]'s. Each
#' further one draws every additive log-ratio of the transition matrix from
#' \eqn{N(0, 2^2)} and adds \eqn{N(0, 0.7^2)} to the logarithm of every gap,
#' and the first level stays where the first start put it.
#' @details
#' A log-ratio of standard deviation 2 gives transition probabilities spread
#' over the interval from near 0 to near 1, so the draws cover chains that stay
#' in a regime and chains that leave it at once. A gap scaled by
#' \eqn{e^{\pm 0.7}}, about a factor of two either way, keeps the levels
#' ordered. The draw for start \eqn{s} is made with the seed
#' \eqn{100 + s - 1}, and the caller's random number generator is restored
#' afterwards, so the same fit gives the same starts and does not move the
#' caller's stream. The log-likelihood of a regime model has several local
#' maxima, which different starts may reach.
#' @param term A [RegimeTerm()].
#' @param ... Passed to [term_start()].
#' @param target The response on the scale of the predictor, or `NULL`.
#' @return A list of `term@n_start` named numeric vectors on the
#'   unconstrained scale.
#' @keywords internal
S7::method(term_starts, RegimeTerm) <- function(term, ..., target = NULL) {
  z0 <- term_start(term, ..., target = target)
  n <- term@n_start
  if (!length(n) || n <= 1L) return(list(z0))
  a <- startsWith(names(z0), "alr")
  g <- startsWith(names(z0), "gap")
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(if (had) {
    assign(".Random.seed", old, envir = globalenv())
  } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  })
  c(list(z0), lapply(seq_len(n - 1L), function(s) {
    set.seed(100L + s)
    z <- z0
    z[a] <- stats::rnorm(sum(a), 0, 2)
    z[g] <- z[g] + stats::rnorm(sum(g), 0, 0.7)
    z
  }))
}

#' @title Build a Regime Term
#' @name term_build.RegimeTerm
#'
#' @description
#' Resolves the grouping and the ordering against the data and records them in
#' the `build_info`. The method does not compute anything else: the term has no
#' design block, its chain structure is fixed by the constructor, and its
#' levels are parameters supplied at each evaluation.
#'
#' @details
#' `by` and `time` are evaluated in `data` with [baseenv()] as the enclosure.
#' Each must give one non-missing value per row; both are checked here, since
#' this is the first point at which they meet data.
#'
#' The `build_info` holds `order`, the row indices of each group sorted by time,
#' and `n`, the observation count. [term_loglik()] runs one forward recursion
#' per group in that order, each from the chain's stationary distribution, so
#' a panel of independent series is fitted by giving `by`.
#'
#' Without `time` the rows are taken in the order in which they appear. The
#' recursion depends on that order, so a data frame that is not already sorted
#' defines a different model. With `time`, the rows of a group are sorted by it,
#' and `time` only orders them: a gap between two times counts as one step of
#' the chain.
#'
#' @param term A [RegimeTerm()].
#' @param data A data frame carrying whatever `by` and `time` name.
#' @param ... Unused.
#'
#' @return The term, in which `build_info` is filled. [term_is_built()] reads
#'   that property on this branch, so it is `TRUE` for the result.
#'
#' @seealso [regime()], [term_loglik()].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:40, id = rep(1:2, each = 20),
#'                  y = c(rnorm(20), rnorm(20, 3)))
#'
#' # The rows form one group in row order.
#' b <- term_build(regime(2, time = t), dd)
#' names(b@build_info)
#' lengths(b@build_info$order)
#'
#' # The rows form two independent series, and each has its own recursion.
#' b2 <- term_build(regime(2, by = id, time = t), dd)
#' lengths(b2@build_info$order)
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
  term@build_info <- list(order = ord, n = n)
  term
}

#' The Stationary Distribution of a Chain, and Its Derivative
#'
#' @description
#' Computes the row vector \eqn{\delta} that solves \eqn{\delta P = \delta} with
#' \eqn{\sum_j \delta_j = 1}, together with its derivative in whatever the
#' transition matrix depends on.
#'
#' @details
#' Differentiating \eqn{\delta(I - P) = 0} under the normalization gives
#' \eqn{d\delta\,(I - P) = \delta\,dP} with \eqn{\sum_j d\delta_j = 0}, so
#' both the value and the derivative come from the same linear system, in
#' which one column is replaced by the normalization.
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
#' @return A list with `delta`; `ddelta`, a list with one numeric vector of
#'   length \eqn{K} per element of `dP`; and, when `d2P` is given, `d2delta`,
#'   a list with one square matrix per state, each of size `length(dP)`.
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
#' @details
#' The density and the score of every observation under every regime are
#' computed first, by one call of each callback per regime with the whole index
#' vector, and the recursion and its first derivatives are then run by the
#' compiled kernel `regime_forward_cpp()`. The derivative in a gap is taken in
#' the gap itself. A gap that is not positive and a `psi` that lacks a
#' parameter signal an error, as does a callback that does not return one value
#' per observation.
#' @param term A built `RegimeTerm`.
#' @param eta The static predictor.
#' @param y The response, reaching the recursion through the two callbacks.
#' @param logdens,score The log-density and its derivative in the predictor, as
#'   functions that each return one value per observation.
#' @param psi The parameters on the parameter scale, named as [term_params()].
#' @param ... Unused.
#' @return A list with `loglik` and `jacobian`.
#' @keywords internal
S7::method(term_loglik, RegimeTerm) <- function(term, eta, y, logdens, score,
                                                psi, ...) {
  bp <- term@build_info
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
#' Returns the probability of each regime at each observation given the whole
#' series, which is everything that a model layer needs to differentiate a
#' likelihood mixed over states.
#'
#' @details
#' [term_loglik()] returns the derivative of the mixed likelihood
#' in the term's own parameters. The derivative in the coefficients of the
#' model's equations is obtained from the smoothed probabilities, by Fisher's
#' identity, and needs no second recursion carrying derivatives. Writing
#' \eqn{L} for the log-likelihood, \eqn{\gamma_t(k)} for the probability
#' returned here and \eqn{\theta_t(k)} for the parameters of the model at
#' observation \eqn{t} under regime \eqn{k},
#'
#' \deqn{\frac{\partial L}{\partial \eta_{q,t}}
#'   = \sum_k \gamma_t(k)\,
#'     \frac{\partial \ell(y_t; \theta_t(k))}{\partial \eta_q},}
#'
#' for every predictor that the model carries, not only the one that the
#' regimes shift. A caller therefore differentiates its own likelihood \eqn{K}
#' times, vectorized, and weights the results, with no callback per observation.
#' This is possible because the regimes shift a predictor that is known before
#' the recursion starts.
#'
#' The probabilities come from a forward pass and a backward pass beside it,
#' and both passes are normalized at every step. Without the normalization the
#' quantities are products of \eqn{t} densities and reach zero in double
#' precision on a long series.
#'
#' @param term A built [RegimeTerm()].
#' @param eta The static predictor.
#' @param y The response.
#' @param logdens The log-density as a function of the predictor and the
#'   observation index, as [term_loglik()] takes it, returning one value per
#'   observation.
#' @param psi The term's parameters on the parameter scale, named as
#'   [term_params()].
#' @param ... Passed to methods.
#'
#' @return A numeric matrix of `n` rows and \eqn{K} columns, in which every
#'   row sums to one and the entry in row \eqn{t} and column \eqn{k} is
#'   \eqn{P(S_t = k \mid y_1, \dots, y_n)}.
#'
#' @seealso [term_loglik()] for the likelihood that the same recursion computes,
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
#' # probability of the second regime marks it.
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
#' Runs in R a forward pass and a backward pass, normalizes both at every step,
#' and returns the probability of each regime at each observation given the
#' whole series.
#' The recursion runs over each group in time order, from the chain's
#' stationary distribution.
#' @param term A built [RegimeTerm()].
#' @param eta The static predictor.
#' @param y The response.
#' @param logdens The log-density, returning one value per observation when
#'   given the whole index vector; any other length signals an error.
#' @param psi The term's parameters on the parameter scale, named as
#'   [term_params()].
#' @param ... Unused.
#' @return A numeric matrix, one row per observation and one column per
#'   regime.
#' @keywords internal
S7::method(term_posterior, RegimeTerm) <- function(term, eta, y, logdens,
                                                   psi, ...) {
  bp <- term@build_info
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
    a <- as.numeric(logdens(eta + mu[jj], idx))
    if (length(a) != n) {
      stop(paste0("'logdens' must return one value per observation when ",
                  "given the whole index vector."), call. = FALSE)
    }
    LF[, jj] <- a
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
#' Returns the log-likelihood of a term that mixes over latent states, with
#' its exact gradient and its exact Hessian in the whole of a caller's unknown
#' vector: the coefficients of every equation together with the term's own
#' parameters.
#'
#' @details
#' The smoothed probabilities of [term_posterior()] give the gradient by
#' Fisher's identity. A caller can also assemble from them the complete-data
#' information, which is the ordinary information averaged over the states
#' and is the matrix that an EM step inverts. The observed information of the
#' mixture is a different matrix, smaller than the complete-data information
#' by the information that the unobserved states cost, so a standard error
#' computed from the complete-data matrix is too small.
#'
#' Louis's identity expresses the difference as the conditional variance of
#' the complete-data score, which needs the pairwise smoothed probabilities
#' and a recursion carrying a second moment. This method uses a shorter
#' route. The scaled forward recursion computes the observed
#' log-likelihood exactly, as \eqn{\log L = \sum_t \log c_t} with \eqn{c_t}
#' the normalizing constant of one step, so differentiating that arithmetic
#' twice gives the observed Hessian without any identity. The recursion carries
#' \eqn{a_t(k)} together with its first and second derivatives in the
#' unknown vector \eqn{u}. The derivatives of the contribution are
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
#' The cost is \eqn{O(nK^2m^2)} in time. The second derivatives of the state
#' take \eqn{O(Km^2)} in storage, beside \eqn{O(nKm)} for the emission scores.
#' The computation is meant to be run once, at a fitted point.
#'
#' @param term A built [RegimeTerm()].
#' @param eta The static part of the predictor that the regimes shift.
#' @param y The response.
#' @param logdens A function of a predictor value and a row index that returns
#'   the log-density there, one value per observation, as [term_loglik()]
#'   takes it.
#' @param grad A function of the same two arguments that returns a matrix with
#'   one row per observation and one column per distribution parameter, whose
#'   columns hold the derivative of the log-density in each predictor.
#' @param hess A function of the same two arguments that returns an array of
#'   dimension `c(n, p, p)`, whose first index is the observation and whose
#'   entries are the second derivatives of the log-density in the predictors,
#'   where `p` is the number of distribution parameters.
#' @param psi The term's parameters on the parameter scale, named as
#'   [term_params()]. The derivatives in the columns `cols` are taken in the
#'   unconstrained coordinates, so the derivative in a gap is taken in its
#'   logarithm, as it is for the scale of [term_links()].
#' @param seed A list with one matrix per distribution parameter. Each matrix
#'   has one row per observation and one column per unknown, and holds the
#'   derivative of that predictor in the caller's unknowns, with zeros in the
#'   columns that the term's own parameters occupy.
#' @param cols The columns of the unknown vector that the term's own parameters
#'   occupy, in the order of [term_params()].
#' @param level The index, among the distribution parameters, of the one that
#'   the regimes shift.
#' @param weights Optional observation weights: one number or one number per
#'   observation, without missing values. `NULL` gives every observation the
#'   weight one.
#' @param ... Passed to methods.
#'
#' @return A list of three: `loglik`, the `n` per-observation contributions,
#'   without the weights; `gradient`, the derivative of their weighted sum in
#'   the caller's whole unknown vector; and `hessian`, the observed Hessian of
#'   the weighted sum in that same vector, an `m` by `m` matrix with `m` its
#'   length. The Hessian is the mixture's own. Its negative, the observed
#'   information, is smaller in the Loewner order than the complete-data
#'   information, and the difference is the information that the unobserved
#'   states cost.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:30, y = c(rnorm(15), rnorm(15, 3)))
#' term <- term_build(regime(2, time = t), dd)
#' # The model is a gaussian mean of unit variance: one predictor, one unknown
#' # besides the term's own four.
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

#' @title Per-Observation Scores of a Structural Term
#' @name term_scores
#' @description
#' Returns, for each observation, the derivative of the weighted increment
#' of the term's likelihood in the caller's whole unknown vector. The sum of
#' the rows is the `gradient` that [term_hessian()] returns.
#'
#' @details
#' For a [regime()] term the increment of observation \eqn{t} is
#' \eqn{w_t \log p(y_t \mid y_1, \dots, y_{t-1})}, with \eqn{p} the
#' conditional density given by the scaled forward recursion, and the row is
#' its derivative. The derivative is carried in forward mode through the
#' recursion, which is the first-order part of the recursion of
#' [term_hessian()]. A term class that does not implement the generic
#' signals an error.
#'
#' @param term A built structural term.
#' @param eta The static predictor.
#' @param y The response.
#' @param logdens,grad The log-density and its first derivatives in the
#'   predictors, as in [term_hessian()].
#' @param psi The term's parameters on the parameter scale.
#' @param seed The derivative of each predictor in the caller's unknowns.
#' @param cols The columns that the term's own parameters occupy.
#' @param level The distribution parameter that the regimes shift.
#' @param weights Optional observation weights: one number or one number per
#'   observation, without missing values. `NULL` gives every observation the
#'   weight one.
#' @param ... Passed to methods.
#'
#' @return A matrix with one row per observation and one column per unknown
#'   of the caller, in the order of the rows of the data. Each row is
#'   multiplied by the weight of its observation.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:30, y = c(rnorm(15), rnorm(15, 3)))
#' term <- term_build(regime(2, time = t), dd)
#' psi <- list(level1 = 0, gap2 = 3, alr1.1 = 2, alr2.1 = -2)
#' seed <- list(matrix(c(rep(1, 30), rep(0, 30 * 4)), 30, 5L))
#' ld <- function(e, i) dnorm(dd$y[i], e, log = TRUE)
#' gr <- function(e, i) matrix(dd$y[i] - e, ncol = 1L)
#' U <- term_scores(term, eta = rep(0, 30), y = dd$y, logdens = ld,
#'                  grad = gr, psi = psi, seed = seed, cols = 2:5,
#'                  level = 1L)
#' dim(U)
#'
#' @seealso [term_hessian()]
#' @export
term_scores <- S7::new_generic("term_scores", "term",
  function(term, eta, y, logdens, grad, psi, seed, cols, level,
           weights = NULL, ...) S7::S7_dispatch())

S7::method(term_scores, structural_term) <- function(term, eta, y, logdens,
                                                     grad, psi, seed, cols,
                                                     level, weights = NULL,
                                                     ...) {
  stop(sprintf("the term class '%s' does not implement term_scores().",
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

# The scaled forward recursion that term_hessian() and term_scores() share.
# It carries the derivative of the filter to first order and, when 'hess' is
# given, to second order. The result holds the log-likelihood increment of
# each observation ('loglik'), the derivative of each weighted increment
# ('scores', one row per observation), their sum accumulated in time order
# ('gradient') and, with 'hess', the Hessian.
.regime_pass <- function(term, eta, y, logdens, grad, hess, psi, seed, cols,
                         level, weights) {
  second <- !is.null(hess)
  bp <- term@build_info
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
    stop(sprintf("'cols' must give the %d columns that the term's parameters occupy.",
                 length(nm)), call. = FALSE)
  }
  if (length(level) != 1L || level < 1L || level > npar) {
    stop("'level' must index one of the distribution parameters.",
         call. = FALSE)
  }
  w <- .structural_weights(weights, n)

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
    a <- as.numeric(logdens(e, idx))
    if (length(a) != n) {
      stop(paste0("'logdens' must return one value per observation when ",
                  "given the whole index vector."), call. = FALSE)
    }
    LF[, j] <- a
    g <- as.matrix(grad(e, idx))
    if (nrow(g) != n || ncol(g) != npar) {
      stop(sprintf("'grad' must return an %d by %d matrix.", n, npar),
           call. = FALSE)
    }
    GP[[j]] <- g
    if (second) {
      h <- hess(e, idx)
      if (length(h) != n * npar * npar) {
        stop(sprintf("'hess' must return an %d by %d by %d array.", n, npar,
                     npar), call. = FALSE)
      }
      HP[[j]] <- array(as.numeric(h), c(n, npar, npar))
    }
    # the emission score in the caller's unknowns, vectorized over rows
    gm <- matrix(0, n, m)
    for (q in seq_len(npar)) gm <- gm + g[, q] * seed[[q]]
    gm[, cols] <- gm[, cols] + outer(g[, level], ch$dmu[j, ])
    GM[[j]] <- gm
  }

  P <- ch$P
  loglik <- numeric(n)
  scores <- matrix(0, n, m)
  gradient <- numeric(m)
  hessian <- if (second) matrix(0, m, m)

  # the rows of the per-observation predictor derivative, shared by every
  # regime except in the columns the term's own parameters occupy
  Abase <- matrix(0, npar, m)

  for (rows in bp$order) {
    a <- ch$delta
    da <- matrix(0, k, m)
    da[, cols] <- t(do.call(rbind, ch$ddelta))
    if (second) {
      d2a <- vector("list", k)
      for (j in seq_len(k)) {
        M <- matrix(0, m, m)
        M[cols, cols] <- ch$d2delta[[j]]
        d2a[[j]] <- M
      }
    }

    for (tt in seq_along(rows)) {
      row <- rows[tt]
      lf <- LF[row, ]
      mx <- max(lf)
      wv <- exp(lf - mx)

      if (second) for (q in seq_len(npar)) Abase[q, ] <- seed[[q]][row, ]

      if (tt == 1L) {
        pred <- a
        dpred <- da
        if (second) d2pred <- d2a
      } else {
        pred <- as.numeric(a %*% P)
        dpred <- crossprod(P, da)
        if (second) d2pred <- vector("list", k)
        for (l in seq_len(k)) {
          Dl <- matrix(ch$dPz[, l, ], k, np)
          # the chain's own contribution: supported on the term's columns
          dpred[l, cols] <- dpred[l, cols] + as.numeric(a %*% Dl)
          if (second) {
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
      }

      atil <- wv * pred
      datil <- matrix(0, k, m)
      if (second) d2atil <- vector("list", k)
      for (l in seq_len(k)) {
        gl <- GM[[l]][row, ]
        datil[l, ] <- wv[l] * (gl * pred[l] + dpred[l, ])
        if (second) {
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
      }

      ct <- sum(atil)
      dct <- colSums(datil)
      gt <- dct / ct
      loglik[row] <- log(ct) + mx
      scores[row, ] <- w[row] * dct / ct
      gradient <- gradient + w[row] * gt
      if (second) {
        d2ct <- matrix(0, m, m)
        for (l in seq_len(k)) d2ct <- d2ct + d2atil[[l]]
        hessian <- hessian + w[row] * (d2ct / ct - outer(gt, gt))
      }

      a <- atil / ct
      da <- (datil - outer(atil, dct) / ct) / ct
      if (second) {
        for (l in seq_len(k)) {
          d2a[[l]] <- (d2atil[[l]] - outer(da[l, ], dct) - outer(dct, da[l, ]) -
                         a[l] * d2ct) / ct
        }
      }
    }
  }

  if (second) hessian <- (hessian + t(hessian)) / 2
  list(loglik = loglik, scores = scores, gradient = gradient,
       hessian = hessian)
}

#' @title Observed Hessian of a Regime Term
#' @name term_hessian.RegimeTerm
#' @description
#' Propagates the derivatives to second order, in forward mode, through the
#' scaled forward recursion and returns the exact Hessian of the mixed
#' log-likelihood.
#' @details
#' The recursion runs in R over each group in time order. The derivatives in
#' the term's own parameters are taken on the unconstrained scale: the first
#' level and the log-ratios as they are, and each gap through its logarithm.
#' The callbacks `logdens`, `grad` and `hess` are called once per regime with
#' the whole index vector, and a result of the wrong length or shape signals an
#' error.
#' @param term A built [RegimeTerm()].
#' @param eta The static predictor.
#' @param y The response.
#' @param logdens,grad,hess The log-density and its first two derivatives in
#'   the predictors. The last returns an array of dimension `c(n, p, p)`.
#' @param psi The term's parameters on the parameter scale.
#' @param seed The derivative of each predictor in the caller's unknowns.
#' @param cols The columns that the term's own parameters occupy.
#' @param level The distribution parameter that the regimes shift.
#' @param weights Observation weights: one number or one per observation,
#'   without missing values.
#' @param ... Unused.
#' @return A list with `loglik`, `gradient` and `hessian`.
#' @keywords internal
S7::method(term_hessian, RegimeTerm) <- function(term, eta, y, logdens, grad,
                                                 hess, psi, seed, cols, level,
                                                 weights = NULL, ...) {
  if (!is.function(hess)) {
    stop("'hess' must be a function.", call. = FALSE)
  }
  r <- .regime_pass(term, eta, y, logdens, grad, hess, psi, seed, cols, level,
                    weights)
  list(loglik = r$loglik, gradient = r$gradient, hessian = r$hessian)
}

#' @title Per-Observation Scores of a Regime Term
#' @name term_scores.RegimeTerm
#' @description
#' Carries the first derivative, in forward mode, through the scaled forward
#' recursion and returns the derivative of each weighted increment.
#' @details
#' The recursion and the chart of the chain are those of
#' [term_hessian()] for a [RegimeTerm()], with the derivative of order two
#' omitted. The increment of an observation is \eqn{\log c_t}, with
#' \eqn{c_t} the normalizing constant of its step, and its derivative is
#' \eqn{\dot c_t / c_t}. The callbacks are called once per regime with the
#' whole index vector.
#' @param term A built [RegimeTerm()].
#' @param eta The static predictor.
#' @param y The response.
#' @param logdens,grad The log-density and its first derivatives in the
#'   predictors.
#' @param psi The term's parameters on the parameter scale.
#' @param seed The derivative of each predictor in the caller's unknowns.
#' @param cols The columns that the term's own parameters occupy.
#' @param level The distribution parameter that the regimes shift.
#' @param weights Observation weights: one number or one per observation,
#'   without missing values.
#' @param ... Unused.
#' @return A matrix with `n` rows and one column per unknown of the caller.
#' @keywords internal
S7::method(term_scores, RegimeTerm) <- function(term, eta, y, logdens, grad,
                                                psi, seed, cols, level,
                                                weights = NULL, ...) {
  .regime_pass(term, eta, y, logdens, grad, NULL, psi, seed, cols, level,
               weights)$scores
}

#' @title Print a Regime Term
#' @name print.RegimeTerm
#'
#' @description
#' Prints two lines: the label and the number of regimes, followed by the
#' term's parameters. The first line of a built term gives the number of groups
#' over which the recursion runs, and the first line of a specification has no
#' group count.
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
#' The group count is `length(build_info$order)`, which is 1 unless the term
#' was given a `by`. The parameter list is [term_params()], which is the same
#' before and after a build because only `k` determines it.
#'
#' The description of a built term gives its count of groups and omits the word
#' "built".
#'
#' @param x A [RegimeTerm()], built or not.
#' @param ... Unused, and accepted so that the signature matches [print()]'s.
#'
#' @return `x`, invisibly. The function is called for the two lines that it
#'   writes.
#'
#' @seealso [regime()], [term_params()].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:40, id = rep(1:2, each = 20),
#'                  y = c(rnorm(20), rnorm(20, 3)))
#'
#' # The first call prints a specification, and the second prints the same
#' # term built over two groups.
#' regime(2)
#' term_build(regime(2, by = id, time = t), dd)
#'
#' # Three regimes carry one level, two gaps and six log-ratios.
#' regime(3)
#'
#' @keywords internal
S7::method(print, RegimeTerm) <- function(x, ...) {
  built <- length(x@build_info) > 0L
  cat(sprintf("<RegimeTerm> '%s': %d regimes%s\n", x@label, x@k,
              if (built) sprintf("; %d group(s)", length(x@build_info$order))
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
#' Draws a path of the latent chain and returns the predictor of each
#' observation, which is the static predictor plus the level of the state that
#' the observation occupies.
#'
#' @details
#' One chain is drawn for each group of the `build_info`, in time order, and each
#' starts from the chain's stationary distribution. [term_loglik()] and
#' [term_posterior()] read the data in this way, and the likelihood is written
#' with that initial law, so
#' starting from a fixed first state or from a uniform draw would simulate a
#' different model from the one that is fitted. The transition matrix comes
#' from the term's own [parameters7::transition_matrix()], and the stationary
#' law is the solution of \eqn{\delta(I - P) = 0} with
#' \eqn{\sum_j \delta_j = 1}, which the method solves itself.
#'
#' The response is not drawn, because the levels do not depend on it, so `y` is
#' returned as `NULL` and the caller draws at the returned predictor. An `eta`
#' whose length differs from the number of rows of the build data signals an
#' error.
#'
#' @param term A built regime term.
#' @param psi The term's parameters, on the parameter scale.
#' @param eta The static part of the predictor, one value per row of the
#'   data on which the term was built.
#' @param draw Ignored: the levels do not depend on the response.
#' @param ... Ignored.
#'
#' @return A list with `eta`, the static predictor plus the level of each
#'   observation's state; `y` (`NULL`); and `latent`, the state of each
#'   observation, in the order of the data.
#'
#' @seealso [term_simulate()], [regime()]
#'
#' @keywords internal
S7::method(term_simulate, RegimeTerm) <- function(term, psi, eta, draw, ...) {
  bp <- term@build_info
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
  if (n != bp$n) {
    stop(sprintf("'eta' must have one value per row of the data on which the term was built (%d).",
                 bp$n), call. = FALSE)
  }
  # one chain per group, in time order, started from the stationary law,
  # as the likelihood reads the data
  st <- integer(n)
  for (rows in bp$order) {
    st[[rows[1L]]] <- sample.int(k, 1L, prob = pmax(delta, 0))
    for (i in seq_len(length(rows) - 1L)) {
      st[[rows[i + 1L]]] <- sample.int(k, 1L,
                                       prob = pmax(P[st[[rows[i]]], ], 0))
    }
  }
  list(eta = as.numeric(eta) + mu[st], y = NULL, latent = st)
}
