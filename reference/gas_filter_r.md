# The Score-Driven Recursion in R

The loop `gas_filter_cpp()` replaces, kept so the compiled route has
something to be compared against that shares none of its code.

## Usage

``` r
gas_filter_r(
  eta,
  order,
  p,
  q,
  omega,
  a,
  b,
  db,
  f0,
  df0,
  i_a,
  np,
  score,
  curvature
)
```

## Arguments

- eta:

  The static predictor.

- order:

  A list of row indices, one entry per group, in time order.

- p, q:

  The score and autoregressive orders.

- omega:

  The level.

- a, b:

  The score loadings and the autoregressive coefficients.

- db:

  The derivative of the coefficients in the parameters.

- f0, df0:

  The starting level and its derivative.

- i_a:

  The positions of the score loadings among the parameters.

- np:

  The number of parameters.

- score, curvature:

  The callbacks of
  [`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md).

## Value

A list with `eta`, `jacobian` and `curv`.
