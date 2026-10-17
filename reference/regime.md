# Markov Regime Switching

Constructs a latent Markov chain of \\K\\ regimes, and each regime
shifts the linear predictor by a level of its own (Hamilton, 1989). The
likelihood is the mixture over the unobserved state path, evaluated by
the forward recursion, and it is built from any density that the model
carries: the term supplies the chain and the levels, the distribution
supplies everything else.

## Usage

``` r
regime(k = 2, by = NULL, time = NULL, label = "regime", n_start = 1)
```

## Arguments

- k:

  The number of regimes, a single finite whole number of at least 2. Any
  other value signals an error.

- by:

  An optional grouping variable, given as a bare expression. Each group
  runs its own recursion from the stationary distribution, so a panel of
  independent series is `by = id`. `NULL`, the default, is one group.

- time:

  An optional ordering variable, a bare expression. `NULL`, the default,
  takes the rows in the order in which they are given. Both `by` and
  `time` must evaluate to one non-missing value per row at build time.
  Within a group, `time` only orders the rows: a gap between two times
  counts as one step of the chain.

- label:

  A single non-empty character string naming the term, `"regime"` by
  default.

- n_start:

  The number of starting points that a fitting layer tries, a single
  finite whole number of at least 1. The likelihood of a regime model
  has several maxima, and a fit started at one point reaches the maximum
  of the basin that contains it. The first start is
  [`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)'s,
  which places the levels at the quantiles of the response when the
  fitting layer supplies it and at zero otherwise. Each further start
  displaces the first as
  [`term_starts()`](https://statmodels7.github.io/modelterms7/reference/term_starts.md)
  describes, and the fit with the best criterion is kept. `1`, the
  default, fits from the first start alone. The model is fitted once
  from each start.

## Value

An unbuilt
[`RegimeTerm()`](https://statmodels7.github.io/modelterms7/reference/RegimeTerm.md):
a specification, whose `blueprint` is empty until
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
resolves the ordering.

## The forward recursion

Writing \\f_j(t)\\ for the density of observation \\t\\ at the predictor
shifted by the level of regime \\j\\,

\$\$\tilde\alpha_t(j) = f_j(t) \sum_i \alpha\_{t-1}(i) P\_{ij}, \qquad
c_t = \sum_j \tilde\alpha_t(j), \qquad \alpha_t = \tilde\alpha_t /
c_t,\$\$

where the predictive distribution of the first observation of each group
is the chain's stationary distribution. The log-likelihood is \\\sum_t
\log c_t\\.
[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
returns those contributions, one per observation, with their exact
derivatives propagated beside the state.

Normalizing at every step keeps the recursion representable. The
unnormalized quantities are products of \\t\\ densities, so they decay
geometrically and reach zero in double precision on a long series.

## Relation to gas()

[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md) is
driven by the score of the density and moves continuously, while a
regime chain moves in jumps between a finite number of states. Both are
built from the density itself, so both apply to any family that the
model carries, and both propagate their exact derivative beside the
state.

## The parameters and their charts

The levels are ordered by construction: `level1` is free on the
identity, and each of the others is the previous one plus a positive
gap, `gap2` ... `gapK`, and every gap is carried on a log link. Without
an ordering the regimes are exchangeable and the likelihood has \\K!\\
identical maxima, so the estimates of the individual regimes cannot be
reported without a labeling.

The transition matrix is
[`parameters7::transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.html),
whose free values are the additive log-ratios of each row, named
`alr1.1`, `alr2.1` and so on, so every row is a probability vector at
any coordinate. A chain of \\K\\ regimes therefore has \\K(K-1)\\ of
them, and
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
returns \\K + K(K-1)\\ names in all: one level, \\K-1\\ gaps and
\\K(K-1)\\ log-ratios.

The initial distribution is the chain's stationary one, which costs no
parameters. Its derivative comes from the same linear system that gives
the distribution itself.

[`term_level_param()`](https://statmodels7.github.io/modelterms7/reference/term_level_param.md)
returns `"level1"`, since a constant added to it shifts every regime and
is the direction that an intercept in the same equation also spans. A
constant leaves the gaps unchanged, because it cannot express a
difference between regimes.

## References

Hamilton, J. D. (1989). A new approach to the economic analysis of
nonstationary time series and the business cycle. *Econometrica*, 57(2),
357–384.

## See also

[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
for the continuous dynamic term,
[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
for the likelihood that it computes,
[`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)
for the smoothed state probabilities,
[`term_hessian()`](https://statmodels7.github.io/modelterms7/reference/term_hessian.md)
for the observed information.

## Examples

``` r
# The term has one level, one gap, and two log-ratios of a two-state chain.
term_params(regime(2))
#> [1] "level1" "gap2"   "alr1.1" "alr2.1"
vapply(term_links(regime(2)), function(l) l@link_name, character(1))
#>     level1       gap2     alr1.1     alr2.1 
#> "identity"      "log" "identity" "identity" 

# A chain of K regimes has K + K(K-1) parameters in all.
c(k3 = length(term_params(regime(3))), expected = 3 + 3 * 2)
#>       k3 expected 
#>        9        9 

# The levels are ordered, so level1 is the one that collides with an
# intercept.
term_level_param(regime(2))
#> [1] "level1"

# Fitted through a model, the term carries the state; on its own it
# computes the mixed likelihood from a density that the caller supplies.
set.seed(1)
dd <- data.frame(t = 1:40, y = c(rnorm(20), rnorm(20, 3)))
b <- term_build(regime(2, time = t), dd)
out <- term_loglik(b, rep(0, 40), dd$y,
                   logdens = function(e, i) dnorm(dd$y[i], e, log = TRUE),
                   score = function(e, i) dd$y[i] - e,
                   psi = list(level1 = 0, gap2 = 3, alr1.1 = 2, alr2.1 = -2))
sum(out$loglik)
#> [1] -59.65148

# The smoothed probability of the second regime marks the change.
pp <- term_posterior(b, rep(0, 40), dd$y,
                     logdens = function(e, i) dnorm(dd$y[i], e, log = TRUE),
                     psi = list(level1 = 0, gap2 = 3, alr1.1 = 2,
                                alr2.1 = -2))
round(pp[c(1, 19, 20, 21, 22, 40), 2], 3)
#> [1] 0.000 0.012 0.071 0.999 1.000 1.000


# The model is fitted. The data are simulated from a known truth, so the
# estimates below can be read against it.
if (requireNamespace("statmodels7", quietly = TRUE)) {
  set.seed(4)
  st <- integer(200)
  st[1] <- 1L
  for (i in 2:200) {
    st[i] <- if (runif(1) < 0.05) 3L - st[i - 1L] else st[i - 1L]
  }
  fd <- data.frame(t = 1:200, y = c(-2, 2)[st] + rnorm(200, sd = 0.7))
  ft <- statmodels7::statmod(y ~ 0 + regime(k = 2, time = t),
                             distributions7::gaussian1_distrib(), fd)
  # The truth has levels at -2 and 2, so the first level is -2 and the
  # gap is 4.
  round(coef(ft)$mu[1:2], 2)
}
#> regime.level1   regime.gap2 
#>         -2.13          4.12 
```
