# Build a Joint Score-Driven Filter Against Data

Builds every term of the filter against `data` and joins them again with
[`gas_joint()`](https://statmodels7.github.io/modelterms7/reference/gas_joint.md),
which checks once more that they walk the same series.

## Arguments

- term:

  A
  [`GasJoint()`](https://statmodels7.github.io/modelterms7/reference/GasJoint.md).

- data:

  A data frame.

- ...:

  Passed to
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
  of each term.

## Value

A built
[`GasJoint()`](https://statmodels7.github.io/modelterms7/reference/GasJoint.md).
