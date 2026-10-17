# Posterior Components of a Marginal Break-Point Term

Returns the component weights that Fisher's identity takes at every
observation: for the step kind, the posterior probability of each side
pattern given the whole sample, one column per pattern of active
break-points; for the continuous kinds, each group's posterior over its
quadrature nodes, repeated down the group's rows and zero-padded to the
widest group. For the continuous kinds the columns of different groups
index different nodes.

## Arguments

- term:

  A built
  [`MarginalBreakTerm()`](https://statmodels7.github.io/modelterms7/reference/MarginalBreakTerm.md).

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

A numeric matrix, one row per observation, rows summing to one.
