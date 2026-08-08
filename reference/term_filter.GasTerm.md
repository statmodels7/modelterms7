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

## Value

A list with `eta` and `jacobian`.
