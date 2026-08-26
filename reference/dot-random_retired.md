# The Arguments random() No Longer Takes

Signals an error naming a removed argument, which the dots would
otherwise swallow in silence.

## Usage

``` r
.random_retired(...)
```

## Arguments

- ...:

  The dots of
  [`random()`](https://statmodels7.github.io/modelterms7/reference/random.md).

## Value

Invisibly `NULL`; called for its error.

## Details

`precision` was a second spelling of a multivariate Gaussian whose
matrix parameter is a precision, and the first spelling did not say
which matrix the structure was. `kinks` was derived from the effects'
distribution by penalties7 all along, and the default here overrode that
derivation, so a Laplace prior declared none and a fitting layer sent
its block to the scheme that cannot solve it.
