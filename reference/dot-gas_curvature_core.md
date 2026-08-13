# The Score-Driven Recursion's Second and Third Derivatives

The body
[`term_curvature`](https://statmodels7.github.io/modelterms7/reference/term_curvature.md)
and
[`term_third`](https://statmodels7.github.io/modelterms7/reference/term_third.md)
share. With `direction` `NULL` it propagates the first two derivatives
of the predictor; with a direction it propagates the third as well,
contracted against it.

## Usage

``` r
.gas_curvature_core(
  term,
  eta,
  y,
  score,
  curvature,
  psi,
  g,
  seed,
  blocks,
  direction = NULL
)
```

## Arguments

- term:

  A built `GasTerm`.

- eta:

  The static part of the predictor.

- y:

  The response.

- score, curvature:

  The callbacks of
  [`term_filter`](https://statmodels7.github.io/modelterms7/reference/term_filter.md).

- psi:

  The parameters on the parameter scale.

- g:

  The weights the contraction is taken against.

- seed:

  The derivative of the static predictor in the unknowns.

- blocks:

  The model's derivative pieces.

- direction:

  The direction, or `NULL` for the second order alone.

## Value

A list with `jacobian` and `curvature`, and `dphi` where a direction was
given.

## Details

The two orders are written here once rather than in a method each. The
third order's recursion reads \\F\\, \\\Phi\\, \\\dot S\\ and \\\ddot
S\\ at every lag, so a separate implementation would carry a second copy
of the first two orders, and the two would drift.
