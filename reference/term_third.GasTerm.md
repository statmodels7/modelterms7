# Third Derivatives of a Score-Driven Predictor

The second derivative of the filter's predictor differentiated once more
along one direction, propagated beside the state exactly as the first
two orders are.

## Arguments

- term:

  A built `GasTerm`.

- eta:

  The static part of the predictor.

- y:

  The response, unused directly.

- score, curvature:

  The callbacks of
  [`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md).

- psi:

  The parameters on the parameter scale.

- g:

  The weights the third derivative is contracted against.

- seed:

  The derivative of the static predictor in the unknowns.

- blocks:

  The model's derivative pieces; see
  [`term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.md).

- direction:

  The direction to contract against.

- ...:

  Unused.

## Value

A list with `jacobian`, `dphi` and `curvature`.

## Details

The recursion gains one state, \\\Psi_t = \partial^3f_t/\partial
u^3\[v\]\\, seeded by the third derivative of the score in the same way
\\\Phi\\ is seeded by its second. Everything else it needs, meaning the
directional derivatives of \\F\\, \\\Phi\\, \\\dot S\\ and \\\ddot S\\,
is a contraction of a quantity the second-order recursion already
carries, so no second recursion is run and no three-index array is
formed.

The chart contributes its own third derivatives: the level and the
loadings through their scalar links, the persistence through
[`.gas_chart_derivs3()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_chart_derivs3.md),
which composes
[`gas_levinson3()`](https://statmodels7.github.io/modelterms7/reference/gas_levinson3.md)
with them.
