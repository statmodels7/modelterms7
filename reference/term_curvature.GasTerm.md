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
  [`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md).

- psi:

  The parameters on the parameter scale, named as
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- g:

  The weights the second derivative is contracted against.

- seed:

  The derivative of the static predictor in the unknowns.

- blocks:

  The model's own derivative pieces; see
  [`term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.md).

- ...:

  Unused.

- score_values, curvature_values:

  The score and curvature of the model's log-density evaluated at the
  current predictors, one value per observation, on the parameter scale
  the callbacks read. Supplying both, together with `blocks_data`,
  routes the second-order recursion of the subformula route through the
  compiled kernel; either `NULL` keeps the R route.

- blocks_data:

  The model's derivative pieces as data instead of as a callback: a list
  with `H` (the mixed second derivatives, one column per distribution
  parameter), `D3` (the third derivatives, one column per parameter
  pair, pair `(r, r2)` at column `(r - 1) * np + r2`), `Vs` (the
  per-parameter jacobian rows of the other equations) and `ap` (the
  filter's own parameter index). Read only by the compiled route.

- threads:

  Threads for the compiled route's group loop, as
  [`numericals7::n_threads()`](https://statmodels7.github.io/numericals7/reference/n_threads.html)
  counts them; 1 is sequential.

## Value

A list with `jacobian` and `curvature`.

## Details

With subformulas the general per-observation route runs instead, and the
second derivative is accumulated on each group's active set rather than
as a square over all the unknowns: a development's coordinate reaches
only the groups where its column is not identically zero, so with
grouping indicators the active set has the same size whether the panel
has ten groups or a thousand. Measured, the full square cost 0.39 s at
124 unknowns over 1600 rows and would have reached about twelve minutes
at five hundred groups.

`blocks` is called with the row of the jacobian restricted to the active
set and with that set, and returns its pieces in the same coordinates. A
callback of the earlier three-argument shape is still accepted and given
the full row, its result being subset here; it costs the quadratic
allocation the restriction exists to avoid.
