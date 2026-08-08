# The Stationary Distribution of a Chain, and Its Derivative

The row vector \\\delta\\ solving \\\delta P = \delta\\ with \\\sum_j
\delta_j = 1\\, together with its derivative in whatever the transition
matrix depends on.

## Usage

``` r
regime_stationary(P, dP)
```

## Arguments

- P:

  A row-stochastic matrix.

- dP:

  A list of derivative matrices of `P`.

## Value

A list with `delta` and `ddelta`, the latter one row per element of
`dP`.

## Details

Differentiating \\\delta(I - P) = 0\\ under the normalization gives
\\d\delta\\(I - P) = \delta\\dP\\ with \\\sum_j d\delta_j = 0\\, so both
the value and the derivative come from the same linear system with one
column replaced by the normalization.
