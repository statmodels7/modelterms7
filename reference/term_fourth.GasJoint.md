# The Fourth Derivative of a Joint Score-Driven Filter

Computes the fourth derivative of the filtered predictors contracted
against the caller's weights and two directions, through the jets of
[`term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.md).

## Arguments

- term:

  A
  [`GasJoint()`](https://statmodels7.github.io/modelterms7/reference/GasJoint.md).

- eta, y, score, curvature:

  Unused, as for
  [`term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.md).

- psi, g, seed, blocks:

  As for
  [`term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.md);
  `blocks` carries `U3` and `U4`.

- directions:

  A list of two directions.

- ...:

  Unused.

## Value

A list with `jacobian`, `dphi`, `dpsi` (per filtered parameter, the
third derivative along the first direction applied to the second) and
`curvature`.
