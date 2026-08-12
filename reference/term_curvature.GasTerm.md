# Second Derivatives of a Score-Driven Predictor

The forward Jacobian of the filter's predictor in a caller's unknowns
and the second derivative contracted against the caller's weights, both
propagated through the recursion beside the state.

## Arguments

- term:

  A built `GasTerm`.

- eta:

  The static part of the predictor.

- y:

  The response, unused directly.

- score, curvature:

  The callbacks of
  [`term_filter`](https://statmodels7.github.io/modelterms7/reference/term_filter.md).

- psi:

  The parameters on the PARAMETER scale, named as
  [`term_params`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- g:

  The weights the second derivative is contracted against.

- seed:

  The derivative of the static predictor in the unknowns.

- blocks:

  The model's own derivative pieces; see
  [`term_curvature`](https://statmodels7.github.io/modelterms7/reference/term_curvature.md).

- ...:

  Unused.

## Value

A list with `jacobian` and `curvature`.
