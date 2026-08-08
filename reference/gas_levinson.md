# The Autoregressive Coefficients Behind the Partial Autocorrelations

The Levinson-Durbin recursion carrying partial autocorrelations onto the
coefficients of a stationary autoregression, with the Jacobian of that
map propagated alongside.

## Usage

``` r
gas_levinson(pacf)
```

## Arguments

- pacf:

  A numeric vector of partial autocorrelations in \\(-1, 1)\\.

## Value

A list with `phi`, the coefficients, and `jacobian`, the matrix of their
derivatives with respect to `pacf`.
