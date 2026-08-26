# Log-Likelihood of a Regime Term

Runs the forward recursion over each group in time order, normalizing at
every step, and returns the per-observation contributions with their
exact derivatives.

## Arguments

- term:

  A built `RegimeTerm`.

- eta:

  The static predictor.

- y:

  The response, reaching the recursion through the two callbacks.

- logdens, score:

  The log-density and its derivative in the predictor.

- psi:

  The parameters, named as
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- ...:

  Unused.

## Value

A list with `loglik` and `jacobian`.
