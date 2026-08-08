# Log-Likelihood Contributions of a Structural Term

The per-observation log-likelihood contributions a structural term
produces, with their derivatives in the term's own parameters. This is
the second shape the structural branch takes, beside
[`term_filter`](https://statmodels7.github.io/modelterms7/reference/term_filter.md):
a term that shifts the predictor implements the filter, and one that
rewrites the likelihood itself – a mixture over latent states, say –
implements this, because its contribution is not a predictor and cannot
be reported as one.

## Usage

``` r
term_loglik(term, eta, y, logdens, score, psi, ...)
```

## Arguments

- term:

  A built structural term.

- eta:

  The static part of the linear predictor.

- y:

  The response.

- logdens:

  A function of a predictor value and a row index, returning the
  log-density of that observation at that predictor.

- score:

  A function of the same two arguments returning the derivative of that
  log-density with respect to the predictor.

- psi:

  The term's parameters, named as
  [`term_params`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- ...:

  Passed to methods.

## Value

A list with `loglik`, one contribution per observation summing to the
term's log-likelihood, and `jacobian`, an `n` by `length(psi)` matrix of
its derivatives.

## See also

[`regime`](https://statmodels7.github.io/modelterms7/reference/regime.md),
[`term_filter`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:40, y = c(rnorm(20), rnorm(20, 3)))
term <- term_build(regime(2, time = t), dd)
out <- term_loglik(term, rep(0, 40), dd$y,
                   logdens = function(e, i) dnorm(dd$y[i], e, log = TRUE),
                   score = function(e, i) dd$y[i] - e,
                   psi = list(level1 = 0, gap2 = 3,
                              alr1.1 = 2, alr2.1 = -2))
sum(out$loglik)
#> [1] -59.65148
```
