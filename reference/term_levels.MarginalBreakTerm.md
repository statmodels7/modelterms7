# Levels of a Marginal Break-Point Term

Returns, for the step kind, the constant shift of each side pattern,
which is the sum of the changes of level over the active break-points of
that pattern. For the continuous kinds each node contributes its own
hinge value to the shift, so the shift varies by observation and the
method returns a matrix aligned with
[`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)'s
columns. The node set depends only on the term's parameters, so the
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
