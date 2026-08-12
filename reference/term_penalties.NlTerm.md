# Penalties of a Nonlinear Term

One entry per penalized parameter, named after it and covering its own
coefficients: the single coefficient of a plain parameter, or the whole
vector of a parameter carrying a subformula. The list is empty when
`penalty = "none"`, and for a specification, whose parameters a formula
does not name until the data say which of them the data supply.

## Arguments

- term:

  A built
  [`NlTerm`](https://statmodels7.github.io/modelterms7/reference/NlTerm.md).

- ...:

  Unused.

## Value

A list of entries, as
[`term_penalties`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
documents.
