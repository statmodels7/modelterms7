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
  [`RegimeTerm()`](https://statmodels7.github.io/modelterms7/reference/RegimeTerm.md).

- eta:

  The static predictor.

- y:

  The response.

- logdens:

  The log-density as a function of the predictor and the observation
  index, as
  [`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
  takes it, returning one value per observation.

- psi:

  The term's parameters on the parameter scale, named as
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- ...:

  Passed to methods.

## Value

A numeric matrix of `n` rows and \\K\\ columns, every row summing to
one, giving \\P(S_t = k \mid y_1, \dots, y_n)\\.

## Details

[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
returns the derivative of the mixed likelihood in the term's own
parameters. The derivative in the coefficients of the model's equations
is obtained from the smoothed probabilities, by Fisher's identity, and
needs no second recursion carrying derivatives. Writing \\L\\ for the
log-likelihood, \\\gamma_t(k)\\ for the probability returned here and
\\\theta_t(k)\\ for the parameters of the model at observation \\t\\
under regime \\k\\,

\$\$\frac{\partial L}{\partial \eta\_{q,t}} = \sum_k \gamma_t(k)\\
\frac{\partial \ell(y_t; \theta_t(k))}{\partial \eta_q},\$\$

for every predictor the model carries, not only the one the regimes
shift. A caller therefore differentiates its own likelihood \\K\\ times,
vectorized, and weights the results, with no callback per observation.
This is possible because the regimes shift a predictor that is known
before the recursion starts.

The probabilities come from a forward pass and a backward pass beside
it, both normalized at every step. Without the normalization the
quantities are products of \\t\\ densities and reach zero in double
precision on a long series.

## See also

[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
for the likelihood the same recursion computes,
[`term_hessian()`](https://statmodels7.github.io/modelterms7/reference/term_hessian.md)
for the observed information,
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
for the model.

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:40, y = c(rnorm(20), rnorm(20, 3)))
term <- term_build(regime(2, time = t), dd)
psi <- list(level1 = 0, gap2 = 3, alr1.1 = 2, alr2.1 = -2)
g <- term_posterior(term, rep(0, 40), dd$y,
                    logdens = function(e, i) dnorm(dd$y[i], e, log = TRUE),
                    psi = psi)

# Every row is a distribution over the regimes.
dim(g)
#> [1] 40  2
range(rowSums(g))
#> [1] 1 1

# The data switch level at observation 20, and the smoothed
# probability of the second regime finds it.
round(g[c(1, 19, 20, 21, 22, 40), 2], 3)
#> [1] 0.000 0.012 0.071 0.999 1.000 1.000
```
