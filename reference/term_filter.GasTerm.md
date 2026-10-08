# Filter a Score-Driven Term

Runs the score-driven recursion over each group in time order and
returns the predictor with its dynamic level added, together with the
exact derivative of that predictor with respect to the term's
parameters, propagated alongside the state.

## Arguments

- term:

  A built `GasTerm`.

- eta:

  The static part of the predictor.

- y:

  The response, unused directly: it reaches the filter through `score`
  and `curvature`.

- score:

  A function of the predictor returning \\\partial\ell/\partial\eta\\
  per observation.

- curvature:

  A function of the predictor returning
  \\\partial^2\ell/\partial\eta^2\\ per observation.

- psi:

  The parameters, named as
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- ...:

  Unused.

- fast:

  The fast context of the caller, or `NULL`: a list with `family` and
  `link` (the names that
  [`distributions7::distrib_scalar_route()`](https://statmodels7.github.io/distributions7/reference/distrib_scalar_route.html)
  and
  [`linkfunctions7::link_scalar_route()`](https://statmodels7.github.io/linkfunctions7/reference/link_scalar_route.html)
  return), `link_par` (the link's own parameters, from the same
  function), `k` (the parameter's 1-based index), `bounds`, `y` and
  `theta` (the per-observation parameters, followed by the
  distribution's constants). Where the C registries of distributions7
  and linkfunctions7 cover the pair, the recursion reads the score and
  the curvature through their scalar entry points instead of the R
  callbacks; where they do not, the context is inert and the callbacks
  run as before. The groups run over threads only when the
  distribution's entries are safe on a worker thread.

- threads:

  How many threads the recursion may use, over groups and only on the
  fast route: a group's filter is independent of the others and its
  writes land on its own rows, so no reduction is split and the result
  does not depend on the count, bit for bit.

## Value

A list with `eta`, `jacobian` and `curv`, the curvature read at each
predictor.
