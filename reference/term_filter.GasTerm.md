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
  [`term_params`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- ...:

  Unused.

- fast:

  The fast context of the caller, or `NULL`: a list with `family` (the
  distribution's S7 class name), `link` (the parameter's link name), `k`
  (the parameter's 1-based index), `bounds`, `y` and `theta` (the
  per-observation parameters). Where the C registries of distributions7
  and linkfunctions7 cover the pair, the recursion reads the score and
  the curvature through their scalar entry points instead of the R
  callbacks, bit-identically; where they do not, the context is inert
  and the callbacks run as before.

- threads:

  How many threads the recursion may use, over GROUPS and only on the
  fast route: a group's filter is independent of the others and its
  writes land on its own rows, so no reduction is split and the result
  does not depend on the count, bit for bit.

## Value

A list with `eta`, `jacobian` and `curv`, the curvature read at each
predictor.
