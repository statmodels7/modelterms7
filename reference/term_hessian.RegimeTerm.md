# Observed Hessian of a Regime Term

Propagates the derivatives to second order, in forward mode, through the
scaled forward recursion and returns the exact Hessian of the mixed
log-likelihood.

## Arguments

- term:

  A built
  [`RegimeTerm()`](https://statmodels7.github.io/modelterms7/reference/RegimeTerm.md).

- eta:

  The static predictor.

- y:

  The response.

- logdens, grad, hess:

  The log-density and its first two derivatives in the predictors. The
  last returns an array of dimension `c(n, p, p)`.

- psi:

  The term's parameters on the parameter scale.

- seed:

  The derivative of each predictor in the caller's unknowns.

- cols:

  The columns that the term's own parameters occupy.

- level:

  The distribution parameter that the regimes shift.

- weights:

  Observation weights: one number or one per observation, without
  missing values.

- ...:

  Unused.

## Value

A list with `loglik`, `gradient` and `hessian`.

## Details

The recursion runs in R over each group in time order. The derivatives
in the term's own parameters are taken on the unconstrained scale: the
first level and the log-ratios as they are, and each gap through its
logarithm. The callbacks `logdens`, `grad` and `hess` are called once
per regime with the whole index vector, and a result of the wrong length
or shape signals an error.
