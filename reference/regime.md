# Markov Regime Switching

A latent Markov chain of \\K\\ regimes, each shifting the linear
predictor by a level of its own (Hamilton, 1989). The likelihood is the
mixture over the unobserved state path, evaluated by the forward
recursion, and it is built from whatever density the model carries: the
term supplies the chain and the levels, the distribution supplies
everything else.

## Usage

``` r
regime(k = 2, by = NULL, time = NULL, label = "regime")
```

## Arguments

- k:

  The number of regimes, a single whole number of at least 2. Anything
  else throws.

- by:

  An optional grouping variable, given as a bare expression. Each group
  runs its own recursion from the stationary distribution, so a panel of
  independent series is `by = id`. `NULL`, the default, is one group.

- time:

  An optional ordering variable, a bare expression. `NULL`, the default,
  takes the rows in the order they are given. Both `by` and `time` must
  evaluate to one non-missing value per row at build time.

- label:

  A single non-empty character string naming the term, `"regime"` by
  default.

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

started at the chain's stationary distribution, with the log-likelihood
\\\sum_t \log c_t\\. Those contributions are what
[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
returns, one per observation, with their exact derivatives propagated
beside the state.

Normalizing at every step is what keeps the recursion representable. The
unnormalized quantities are products of \\t\\ densities, so they decay
geometrically and reach zero in double precision on a series of a few
hundred observations.

## It is the complement of the other dynamic term

[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md) is
driven by the score of the density and moves continuously; a regime
chain moves in jumps between a finite number of states. Both are built
from the density itself, so both apply to any family the model carries,
and both propagate their exact derivative beside the state.

## The parameters and their charts

The levels are **ordered by construction**: `level1` is free on the
identity, and each of the others is the previous one plus a positive
gap, `gap2` ... `gapK`, each carried on a log link. Without an ordering
the regimes are exchangeable and the likelihood has \\K!\\ identical
maxima, which is not a hard problem to fit but is one whose answer
cannot be reported.

The transition matrix is
[`parameters7::transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.html),
whose free values are the additive log-ratios of each row, named
`alr1.1`, `alr2.1` and so on, so every row is a probability vector at
any coordinate. A chain of \\K\\ regimes therefore has \\K(K-1)\\ of
them, and
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
returns \\K + K(K-1)\\ names in all: one level, \\K-1\\ gaps and the
rest.

The initial distribution is the chain's stationary one, which costs no
parameters; its derivative comes from the linear system it solves.

[`term_level_param()`](https://statmodels7.github.io/modelterms7/reference/term_level_param.md)
answers `"level1"`, since a constant added to it shifts every regime and
is the direction an intercept in the same equation also spans. The gaps
are unaffected: what a constant cannot express is a difference between
regimes.

## References

Hamilton, J. D. (1989). A new approach to the economic analysis of
nonstationary time series and the business cycle. *Econometrica*, 57(2),
357–384.

## See also

[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
for the continuous dynamic term,
[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
for the likelihood it computes,
[`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)
for the smoothed state probabilities,
[`term_hessian()`](https://statmodels7.github.io/modelterms7/reference/term_hessian.md)
for the observed information.

## Examples

``` r
# One level, one gap, and two log-ratios of a two-state chain.
term_params(regime(2))
#> [1] "level1" "gap2"   "alr1.1" "alr2.1"
vapply(term_links(regime(2)), function(l) l@link_name, character(1))
#>     level1       gap2     alr1.1     alr2.1 
#> "identity"      "log" "identity" "identity" 

# K + K(K-1) parameters in all.
c(k3 = length(term_params(regime(3))), expected = 3 + 3 * 2)
#>       k3 expected 
#>        9        9 

# The levels are ordered, so level1 is the one an intercept collides with.
term_level_param(regime(2))
#> [1] "level1"

# Fitted through a model, the term is what carries the state; on its own
# it computes the mixed likelihood from a density the caller supplies.
set.seed(1)
dd <- data.frame(t = 1:40, y = c(rnorm(20), rnorm(20, 3)))
b <- term_build(regime(2, time = t), dd)
out <- term_loglik(b, rep(0, 40), dd$y,
                   logdens = function(e, i) dnorm(dd$y[i], e, log = TRUE),
                   score = function(e, i) dd$y[i] - e,
                   psi = list(level1 = 0, gap2 = 3, alr1.1 = 2, alr2.1 = -2))
sum(out$loglik)
#> [1] -59.65148

# And the smoothed probability of the second regime finds the change.
pp <- term_posterior(b, rep(0, 40), dd$y,
                     logdens = function(e, i) dnorm(dd$y[i], e, log = TRUE),
                     psi = list(level1 = 0, gap2 = 3, alr1.1 = 2,
                                alr2.1 = -2))
round(pp[c(1, 19, 20, 21, 22, 40), 2], 3)
#> [1] 0.000 0.012 0.071 0.999 1.000 1.000


# Fitted. The data are simulated from a known truth, so the
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
  # truth: levels at -2 and 2, so a first level of -2 and a gap of 4
  round(coef(ft)$mu[1:2], 2)
}
#> regime.level1   regime.gap2 
#>         -2.13          4.12 
```
