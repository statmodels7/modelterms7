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

  A function `function(e, i)` of the predictor `e` and the 1-based row
  index `i`, returning the driving score at that row: the derivative
  \\\partial\ell/\partial\eta\\ where the term's `scaling` is zero, and
  \\\mathcal{I}^{-d}\\\partial\ell/\partial\eta\\ otherwise. The method
  applies no scaling of its own.

- curvature:

  A function `function(e, i)` returning the derivative of the driving
  score with respect to the predictor at that row:
  \\\partial^2\ell/\partial\eta^2\\ where `scaling` is zero.

- psi:

  The parameters, named as
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md);
  each scalar coordinate is on its parameter scale (the loading itself,
  the partial autocorrelation itself). The coefficients of a development
  enter as they are, because the link of the developed parameter is
  applied inside the development. The columns of the returned Jacobian
  refer to the same scale.

- ...:

  Unused.

- fast:

  The fast context of the caller, or `NULL`: a list with `family` and
  `link` (the names that
  [`distributions7::distrib_scalar_route()`](https://statmodels7.github.io/distributions7/reference/distrib_scalar_route.html)
  and
  [`linkfunctions7::link_scalar_route()`](https://statmodels7.github.io/linkfunctions7/reference/link_scalar_route.html)
  return), `link_par` (the link's own parameters, from the same
  function), `k` (the parameter's 1-based index), `bounds`, `y`, `theta`
  (the per-observation parameters, followed by the distribution's
  constants) and, where the term's `scaling` is not zero, `scaling`.
  Where the C registries of distributions7 and linkfunctions7 cover the
  pair, the recursion reads the score and the curvature through their
  scalar entry points instead of the R callbacks; where they do not, the
  context is inert and the callbacks run as before. The groups run over
  threads only when the distribution's entries are safe on a worker
  thread and there are at least eight groups.

- threads:

  How many threads the recursion may use, over groups and only on the
  fast route: a group's filter is independent of the others and writes
  only to its own rows, so the groups are not combined by a reduction
  and the result does not depend on the count, bit for bit.

## Value

A list with `eta`, `jacobian` and `curv`, the curvature read at each
predictor.
