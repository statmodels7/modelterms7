# The Autoregressive Coefficients Behind the Partial Autocorrelations

Computes the Levinson-Durbin recursion, which carries partial
autocorrelations onto the coefficients of a stationary autoregression,
and propagates the Jacobian of that map alongside.

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
