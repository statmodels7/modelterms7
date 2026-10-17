# Third Derivatives of a Score-Driven Predictor

Differentiates the second derivative of the filter's predictor once more
along one direction and propagates the result beside the state exactly
as the first two orders are propagated.

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

  The weights against which the third derivative is contracted.

- seed:

  The derivative of the static predictor in the unknowns.

- blocks:

  The model's derivative pieces; see
  [`term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.md).

- direction:

  The direction against which to contract.

- ...:

  Unused.

- score_values, curvature_values, blocks_data:

  As for
  [`term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.md).
  If all three are supplied and `blocks_data` carries the fourth
  derivatives `D4`, the recursion runs through the compiled kernel.

- threads:

  Threads for the compiled route's group loop.

## Value

A list with `jacobian`, `dphi` and `curvature`.

## Details

The recursion gains one state, \\\Psi_t = \partial^3f_t/\partial
u^3\[v\]\\, seeded by the third derivative of the score in the same way
\\\Phi\\ is seeded by its second. Everything else that it needs, meaning
the directional derivatives of \\F\\, \\\Phi\\, \\\dot S\\ and \\\ddot
S\\, is a contraction of a quantity that the second-order recursion
already carries. The third order therefore needs no second recursion and
forms no three-index array.

The chart contributes its own third derivatives: the level and the
loadings through their scalar links, the persistence through
[`.gas_chart_derivs3()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_chart_derivs3.md),
which composes
[`gas_levinson3()`](https://statmodels7.github.io/modelterms7/reference/gas_levinson3.md)
with them.
