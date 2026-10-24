# Penalties of a Joint Score-Driven Filter

Returns the penalties that the subformulas of every term of the filter
declare, from
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
of the term, each under the key `parameter.key` and with its indices
moved to the term's position in
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
of the joint filter. The scalar parameters and the cross loadings are
unpenalized, so the list is empty when no term carries a subformula.

## Arguments

- term:

  A built
  [`GasJoint()`](https://statmodels7.github.io/modelterms7/reference/GasJoint.md).

- ...:

  Unused.

## Value

A list of entries, as
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
documents.
