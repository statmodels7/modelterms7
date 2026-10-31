# The Levinson-Durbin Map's Derivative of Any Order, Contracted

Returns, per coefficient, the derivative of order `length(dirs) + 2` of
the map, contracted against the directions in `dirs` and leaving two
slots free, as a `q` by `q` matrix. With no direction this is the
Hessian.

## Usage

``` r
.gas_levinson_contract(rho, dirs)
```

## Arguments

- rho:

  The partial autocorrelations.

- dirs:

  A list of zero to three vectors as long as `rho`.

## Value

A list of one matrix per coefficient.
