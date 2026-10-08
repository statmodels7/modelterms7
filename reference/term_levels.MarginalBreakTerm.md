# Levels of a Marginal Break-Point Term

For the step kind, the constant shift of each side pattern, the sums of
the changes of level over the active break-points. For the continuous
kinds the shift varies by observation, each node contributing its own
hinge value, and a matrix is returned, aligned with
[`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)'s
columns. The node set depends on the term's parameters alone, so the
family is not evaluated.

## Arguments

- term:

  A built
  [`MarginalBreakTerm()`](https://statmodels7.github.io/modelterms7/reference/MarginalBreakTerm.md).

- psi:

  The term's parameters.

- eta, y, logdens:

  Accepted for the interface of
  [`term_levels()`](https://statmodels7.github.io/modelterms7/reference/term_levels.md)
  and not used.

- ...:

  Unused.

## Value

A numeric vector, or a matrix with one row per observation.
