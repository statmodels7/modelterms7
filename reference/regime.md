# Markov Regime Switching

A latent Markov chain of \\K\\ regimes, each shifting the linear
predictor by a level of its own (hamilton1989). The likelihood is the
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

  The number of regimes, at least 2.

- by:

  An optional grouping variable; each group runs its own recursion from
  the stationary distribution.

- time:

  An optional ordering variable.

- label:

  A single non-empty string naming the term.

## Value

An object of class
[`RegimeTerm`](https://statmodels7.github.io/modelterms7/reference/RegimeTerm.md)
(a specification; see
[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

## Details

Writing \\f_j(t)\\ for the density of observation \\t\\ at the predictor
shifted by the level of regime \\j\\, the forward recursion is
\$\$\tilde\alpha_t(j) = f_j(t) \sum_i \alpha\_{t-1}(i) P\_{ij}, \qquad
c_t = \sum_j \tilde\alpha_t(j), \qquad \alpha_t = \tilde\alpha_t /
c_t,\$\$ started at the chain's stationary distribution, and the
log-likelihood is \\\sum_t \log c_t\\. Normalizing at every step is what
keeps the recursion representable: the unnormalized quantities are
products of \\t\\ densities, so they decay geometrically and reach zero
in double precision on a series of a few hundred observations. The
contributions \\\log c_t\\ are what
[`term_loglik`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
returns, one per observation, together with their exact derivatives,
propagated through the recursion beside the state.

This is the second dynamic model of the package and it is the complement
of the first:
[`gas`](https://statmodels7.github.io/modelterms7/reference/gas.md) is
driven by the score of the density and moves continuously, while a
regime chain moves in jumps between a finite number of states. Both are
built from the density rather than from an error structure, so both
apply to any family the model carries.

### The parameters and their charts

The levels are **ordered by construction**: the first is free and each
of the others is the previous one plus a positive gap, carried on a log
link. Without an ordering the regimes are exchangeable and the
likelihood has \\K!\\ identical maxima, which is not a hard problem to
fit but is one whose answer cannot be reported. The transition matrix is
[`transition_matrix`](https://statmodels7.github.io/parameters7/reference/transition_matrix.html),
whose free values are the additive log-ratios of each row, so every row
is a probability vector by construction.

The initial distribution is the chain's stationary one, which costs no
parameters and whose derivative is obtained from the linear system it
solves.

## References

Hamilton, J. D. (1989). A new approach to the economic analysis of
nonstationary time series and the business cycle. *Econometrica*, 57(2),
357–384.

## Examples

``` r
term_params(regime(2))
#> [1] "level1" "gap2"   "alr1.1" "alr2.1"
```
