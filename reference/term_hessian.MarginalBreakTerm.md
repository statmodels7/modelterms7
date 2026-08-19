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

Analytic throughout. The step kind propagates first and second
derivatives through the side chain's forward recursion, the prior's
interval-mass derivatives closed for the gaussian and read off the cdf
surface for an explicit prior (whose own degrees-of-freedom column
carries that surface's documented single stencil, the one non-closed
piece anywhere). The continuous kinds differentiate the node sum twice;
the moving panels below the data are affine in the prior's parameters,
so their motion enters the chain rule with no curvature of its own. The
one-break-point gaussian step keeps the interval-sum route, independent
arithmetic the tests hold the propagation to.

The marginal likelihood of a group does not factorize over its
observations, so an observation weight has a reading only when it is
constant within each group; anything else is rejected.
