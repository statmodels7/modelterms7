# The Starting Points of a Joint Score-Driven Filter

Returns the first start,
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)'s,
and then, for each term of the filter in turn, one start for each of
that term's further starts
([`term_starts()`](https://statmodels7.github.io/modelterms7/reference/term_starts.md)
of the term), with the other terms at their first start. The number of
starts is therefore one plus the sum, over the terms, of their starts
beyond the first.

## Arguments

- term:

  A
  [`GasJoint()`](https://statmodels7.github.io/modelterms7/reference/GasJoint.md).

- ...:

  Unused.

- target:

  Unused.

## Value

A non-empty list of named numeric vectors on the unconstrained scale,
named as
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).
