# Where a Joint Score-Driven Filter's Parameters Start

Returns the start of every term of the filter, from
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
of the term, named as
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
of the joint filter names the parameters. The cross loadings start at
zero, which is the filter without them.

## Arguments

- term:

  A
  [`GasJoint()`](https://statmodels7.github.io/modelterms7/reference/GasJoint.md).

- ...:

  Unused.

## Value

A named numeric vector on the unconstrained scale.
