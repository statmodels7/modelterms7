# Smoothed State Probabilities of a Regime Term

The forward pass of
[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
with a backward pass beside it, both normalized, giving the probability
of each regime at each observation given the whole series.

## Arguments

- term:

  A built
  [`RegimeTerm()`](https://statmodels7.github.io/modelterms7/reference/RegimeTerm.md).

- eta:

  The static predictor.

- y:

  The response.

- logdens:

  The log-density.

- psi:

  The term's parameters.

- ...:

  Unused.

## Value

A numeric matrix, one row per observation and one column per regime.
