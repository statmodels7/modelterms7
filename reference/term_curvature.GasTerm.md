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

  The weights against which the second derivative is contracted.

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
  that the callbacks read. Supplying both, together with `blocks_data`,
  routes the recursion through the compiled kernel, on the route without
  subformulas as well as on the route with them. If any of the three is
  `NULL`, the R route runs.

- blocks_data:

  The model's derivative pieces as data instead of as a callback: a list
  with `H` (the mixed second derivatives, one column per distribution
  parameter), `D3` (the third derivatives with the filter's own
  parameter as first index, one column per pair of the other two
  indices, pair `(r, r2)` at column `(r - 1) * n_par + r2`, with `n_par`
  the number of distribution parameters), `Vs` (the per-parameter
  jacobian rows of the other equations) and `ap` (the filter's own
  parameter index), and for the third and fourth orders `D4` and `D5`,
  one column per index triple and quadruple after the filter's own, laid
  out the same way. Read only by the compiled route.

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
has ten groups or a thousand. Accumulating the square over all the
unknowns would grow with the square of their number, and the number of
unknowns grows with the number of groups.

`blocks` is called as `blocks(e, i, v, act)`, with the row `v` of the
Jacobian restricted to the active set `act`, and returns its pieces in
the coordinates of that set. A callback of three arguments,
`blocks(e, i, v)`, is also accepted: it is given the full row and its
result is subset here, which costs the quadratic allocation that the
restriction is meant to avoid.
