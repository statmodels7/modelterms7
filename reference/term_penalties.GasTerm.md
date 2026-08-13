# Penalties of a Score-Driven Term

The penalties the subformulas' sub-terms declare, each under the key
`parameter::subterm` and covering that sub-term's coefficients in the
term's own numbering. The scalar parameters are unpenalized, and the
list is empty for a specification, whose developments do not exist until
the term is built.

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
