# The Score-Driven Recursion's Second, Third and Fourth Derivatives

The body that
[`term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.md),
[`term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.md)
and
[`term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.md)
share. With `direction` `NULL` it propagates the first two derivatives
of the predictor. With one direction it propagates the third as well,
contracted against that direction, and with a list of two directions it
propagates the fourth, contracted against both.

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
  direction = NULL,
  score_values = NULL,
  curvature_values = NULL,
  blocks_data = NULL,
  threads = 1L
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
  [`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md).

- psi:

  The parameters on the parameter scale.

- g:

  The weights against which the contraction is taken.

- seed:

  The derivative of the static predictor in the unknowns.

- blocks:

  The model's derivative pieces.

- direction:

  `NULL` for the second order alone, one numeric vector for the third
  order, or a list of two numeric vectors for the fourth.

- score_values, curvature_values, blocks_data:

  As for
  [`term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.md);
  the compiled route is taken when all three are given (and
  `blocks_data` carries `D4` for the third order and `D4` and `D5` for
  the fourth).

- threads:

  Threads for the compiled route's group loop.

## Value

A list with `jacobian` and `curvature`. With one direction it also has
`dphi`, a matrix, and with two directions it has `dphi`, a list of two
matrices, and `dpsi`.

## Details

The three orders are written here once, in place of a method each. The
third order's recursion reads \\F\\, \\\Phi\\, \\\dot S\\ and \\\ddot
S\\ at every lag, and the fourth order reads the third, so a separate
implementation of each order would carry a copy of the lower orders, and
the copies could diverge.
