# The Reported Quantities of a Joint Score-Driven Filter

Returns the quantities of every term of the filter, from
[`term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.md)
of the term, with the names prefixed by the term's distribution
parameter and the Jacobian placed in the term's own columns.

## Arguments

- term:

  A
  [`GasJoint()`](https://statmodels7.github.io/modelterms7/reference/GasJoint.md).

- zeta:

  The parameters on the unconstrained scale, named as
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
  of the joint filter.

- ...:

  Unused.

## Value

A list with `name`, `value`, `jacobian` and `scale`, as
[`term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.md)
returns them.
