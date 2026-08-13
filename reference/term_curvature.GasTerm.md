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

With deviations the recursion runs once per group on that group's own
parameters, which are the population values plus the group's deviations
on the unconstrained scale. That map is affine, so its second derivative
is zero and the only change is the lift carrying a base coordinate into
two columns of the caller's unknowns, the population value and the
group's own deviation.

The second derivative is accumulated on that group's ACTIVE SET and
never as a square over all the unknowns. A group's rows reach the
coefficients, the population parameters and that group's own deviations,
and nothing else, so the active set has the same size whether the panel
has ten groups or a thousand: the per-observation work is constant in
the number of groups where forming the full square made it quadratic.
Measured, the square cost 0.39 s at 124 unknowns over 1600 rows and
would have reached about twelve minutes at five hundred groups.

`blocks` is therefore called with the row of the jacobian RESTRICTED to
the active set and with that set, and returns its pieces in the same
coordinates. A callback of the earlier three-argument shape is still
accepted and given the full row, its result being subset here; it costs
the quadratic allocation the restriction exists to avoid.
