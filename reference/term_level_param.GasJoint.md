# The Levels of a Joint Score-Driven Filter

Returns the level of every term of the filter, such as `mu.omega` and
`sigma.omega`. Each adds a constant to the predictor of its own
parameter, so a fitting layer compares each with the design of that
parameter's equation.

## Arguments

- term:

  A
  [`GasJoint()`](https://statmodels7.github.io/modelterms7/reference/GasJoint.md).

- ...:

  Unused.

## Value

A character vector with one entry per filtered parameter.
