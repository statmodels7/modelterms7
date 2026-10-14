# Removed Arguments of random()

Signals an error for any argument passed in the dots of
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md),
naming a removed argument where there is one.

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

`precision` and `kinks` are no longer accepted. A structured precision
is the matrix parameter of `mvgaussian2_distrib()`, and the kinks of a
prior are derived from its log-density by penalties7.
