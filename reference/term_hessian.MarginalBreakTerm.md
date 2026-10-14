# Observed Hessian of a Marginal Break-Point Term

The exact observed Hessian of the marginal log-likelihood over a
caller's unknowns, the coefficients of every equation together with the
term's own parameters on the unconstrained scale.

## Arguments

- term:

  A built
  [`MarginalBreakTerm()`](https://statmodels7.github.io/modelterms7/reference/MarginalBreakTerm.md).

- eta:

  The static predictor of the level equation.

- y:

  The response.

- logdens, grad, hess:

  The log-density and its first two derivatives in the predictors:
  `grad` returns a matrix of one row per observation and one column per
  distribution parameter, and `hess` an array of dimensions observations
  by parameters by parameters.

- psi:

  The term's parameters.

- seed:

  A list with one matrix per distribution parameter, each of one row per
  observation, holding the derivative of that parameter's predictor in
  the caller's unknowns.

- cols:

  The columns of the `seed` matrices that the term's own parameters
  occupy, in the order of
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- level:

  The index of the distribution parameter that the term shifts.

- weights:

  Observation weights: one number, or one number per observation,
  constant within each group.

- ...:

  Unused.

## Value

A list with `loglik`, the contribution of each observation without
weights, and `gradient` and `hessian`, which are weighted.

## Details

Analytic throughout. The step kind propagates first and second
derivatives through the side chain's forward recursion, the prior's
interval-mass derivatives closed for the gaussian and read off the cdf
surface of distributions7 for an explicit prior. The column of an
explicit prior's degrees of freedom comes from a single
finite-difference stencil of that surface, the one piece that is not in
closed form. The continuous kinds differentiate the node sum twice; the
moving panels below the data are affine in the parameters of the prior,
so their motion enters the chain rule with no curvature of its own. The
one-break-point gaussian step uses an interval sum that shares no
arithmetic with the propagation, and the tests compare the two.

The marginal likelihood of a group does not factorize over its
observations, so an observation weight has a meaning only when it is
constant within each group. Weights must be one number or one number per
observation, without missing values, and weights that vary within a
group are rejected.
