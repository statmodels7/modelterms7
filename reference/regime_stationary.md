# The Stationary Distribution of a Chain, and Its Derivative

The row vector \\\delta\\ solving \\\delta P = \delta\\ with \\\sum_j
\delta_j = 1\\, together with its derivative in whatever the transition
matrix depends on.

## Usage

``` r
regime_stationary(P, dP, d2P = NULL)
```

## Arguments

- P:

  A row-stochastic matrix.

- dP:

  A list of derivative matrices of `P`.

- d2P:

  Optionally, the second derivatives: one list per row of `P` whose
  elements are the matrices of second derivatives of that row's entries,
  indexed as `dP` is.

## Value

A list with `delta` and `ddelta`, the latter one row per element of
`dP`, and, when `d2P` is given, `d2delta`, one square matrix per state.

## Details

Differentiating \\\delta(I - P) = 0\\ under the normalization gives
\\d\delta\\(I - P) = \delta\\dP\\ with \\\sum_j d\delta_j = 0\\, so both
the value and the derivative come from the same linear system with one
column replaced by the normalization.

Differentiating a second time gives \\d^2\delta\\(I - P) =
d_i\delta\\d_jP + d_j\delta\\d_iP + \delta\\d^2P\\ under \\\sum_j
d^2\delta_j = 0\\, so the second derivative comes from the same system
again.
