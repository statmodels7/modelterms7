# Penalties of a Score-Driven Term

One entry per parameter carrying deviations, named after it and covering
its deviations across the groups. The population parameters are
unpenalized, and the list is empty when `penalty = "none"`, and for a
specification, whose deviations do not exist until the data say how many
groups there are.

## Arguments

- term:

  A built
  [`GasTerm`](https://statmodels7.github.io/modelterms7/reference/GasTerm.md).

- ...:

  Unused.

## Value

A list of entries, as
[`term_penalties`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
documents.
