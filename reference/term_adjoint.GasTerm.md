# Filter a Score-Driven Term Backwards

Runs the recursion of
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
in reverse, returning the derivative of a caller's objective with
respect to the static predictor that the caller supplies and with
respect to the sequence of scores that the callbacks return.

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

  The parameters, named as
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- g:

  The direct derivative of the objective with respect to the predictor
  that the filter produces, one value per observation.

- ...:

  Unused.

- fast, threads:

  The fast context and the thread count of
  [`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md),
  passed to the forward pass that the adjoint re-runs. The reverse pass
  reads the curvature sequence that this pass returns, so with a covered
  context the adjoint evaluates no R callback; without one the callbacks
  are evaluated in the forward pass only, once per observation.

## Value

A list with `deta` and `dscore`.
