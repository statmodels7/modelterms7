# Third Derivatives of a Joint Score-Driven Filter

Differentiates the second derivative of the filtered predictors once
more along one direction, through the jets of
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
  `blocks` carries `U3`.

- direction:

  The direction against which to contract.

- ...:

  Unused.

## Value

A list with `jacobian`, `dphi` (per filtered parameter, a list with one
matrix per direction holding the second derivative of the predictor
applied to the direction) and `curvature`.
