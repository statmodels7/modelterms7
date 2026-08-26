# Observed Hessian of a Regime Term

Forward-mode second-order propagation through the scaled forward
recursion, giving the exact Hessian of the mixed log-likelihood.

## Arguments

- term:

  A built
  [`RegimeTerm()`](https://statmodels7.github.io/modelterms7/reference/RegimeTerm.md).

- eta:

  The static predictor.

- y:

  The response.

- logdens, grad, hess:

  The log-density and its first two derivatives in the predictors.

- psi:

  The term's parameters.

- seed:

  The derivative of each predictor in the caller's unknowns.

- cols:

  The columns the term's own parameters occupy.

- level:

  The distribution parameter the regimes shift.

- weights:

  Observation weights.

- ...:

  Unused.

## Value

A list with `loglik`, `gradient` and `hessian`.
