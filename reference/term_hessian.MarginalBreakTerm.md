# Observed Hessian of a Marginal Break-Point Term

The exact observed Hessian of the marginal log-likelihood over a
caller's unknowns, the coefficients of every equation together with the
term's own parameters on the unconstrained scale.

## Arguments

- term:

  A built
  [`MarginalBreakTerm`](https://statmodels7.github.io/modelterms7/reference/MarginalBreakTerm.md).

- eta:

  The static predictor of the level equation.

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

  The distribution parameter the term shifts.

- weights:

  Observation weights, constant within each group.

- ...:

  Unused.

## Value

A list with `loglik`, `gradient` and `hessian`.

## Details

The step kind with one gaussian break-point differentiates the interval
sum twice, the second derivatives of the conditional collapsing into
per-observation Hessians weighted by each side's posterior probability
and the mass curvature closed in the normal density. Every other
configuration differences the analytic full gradient once – a single
central stencil on the analytic order below, the licence the toolkit's
non-closed derivatives run on – and the one-break-point gaussian route
is the control the tests hold it to.

The marginal likelihood of a group does not factorize over its
observations, so an observation weight has a reading only when it is
constant within each group; anything else is rejected.
