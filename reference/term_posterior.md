# The Smoothed State Probabilities of a Latent Markov Term

The probability of each regime at each observation given the whole
series, which is everything a model layer needs to differentiate a
likelihood mixed over states.

## Usage

``` r
term_posterior(term, eta, y, logdens, psi, ...)
```

## Arguments

- term:

  A built
  [`RegimeTerm`](https://statmodels7.github.io/modelterms7/reference/RegimeTerm.md).

- eta:

  The static predictor.

- y:

  The response.

- logdens:

  The log-density as a function of the predictor and the observation
  index, as
  [`term_loglik`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
  takes it.

- psi:

  The term's parameters, named as
  [`term_params`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- ...:

  Passed to methods.

## Value

A numeric matrix with one row per observation and one column per regime,
whose rows sum to one.

## Details

[`term_loglik`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
returns the derivative of the mixed likelihood in the term's OWN
parameters, which is what estimating those needs. It is not what
estimating the coefficients needs, and for this term the missing piece
is not a second recursion carrying derivatives: it is one quantity, by
Fisher's identity. Writing \\\gamma_t(k)\\ for the probability returned
here and \\\theta_t(k)\\ for the parameters the model has at observation
\\t\\ under regime \\k\\,

\$\$\frac{\partial L}{\partial \eta\_{q,t}} = \sum_k \gamma_t(k)\\
\frac{\partial \ell(y_t; \theta_t(k))}{\partial \eta_q},\$\$

for EVERY predictor the model carries, not only the one the regimes
shift. A caller therefore differentiates its own likelihood \\K\\ times
vectorized and weights the results, and needs no callback per
observation: the regimes shift a predictor that is known before the
recursion starts, which is the property that made the forward pass
compilable and makes this cheap.

The probabilities come from the forward pass this term already runs and
a backward pass beside it, both normalized at every step. Without the
normalization the quantities are products of \\t\\ densities and reach
zero in double precision within a few hundred observations.

## See also

[`term_loglik`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:40, y = c(rnorm(20), rnorm(20, 3)))
term <- term_build(regime(2, time = t), dd)
g <- term_posterior(term, rep(0, 40), dd$y,
                    logdens = function(e, i) dnorm(dd$y[i], e, log = TRUE),
                    psi = list(level1 = 0, gap2 = 3,
                               alr1.1 = 2, alr2.1 = -2))
round(head(g, 3), 4)
#>        [,1]  [,2]
#> [1,] 0.9998 2e-04
#> [2,] 0.9996 4e-04
#> [3,] 1.0000 0e+00
```
