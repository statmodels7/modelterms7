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

## Details

With subformulas the general per-observation route runs instead, and the
second derivative is accumulated on each group's ACTIVE SET rather than
as a square over all the unknowns: a development's coordinate reaches
only the groups where its column is not identically zero, so with
grouping indicators the active set has the same size whether the panel
has ten groups or a thousand. Measured, the full square cost 0.39 s at
124 unknowns over 1600 rows and would have reached about twelve minutes
at five hundred groups.

`blocks` is called with the row of the jacobian RESTRICTED to the active
set and with that set, and returns its pieces in the same coordinates. A
callback of the earlier three-argument shape is still accepted and given
the full row, its result being subset here; it costs the quadratic
allocation the restriction exists to avoid.
