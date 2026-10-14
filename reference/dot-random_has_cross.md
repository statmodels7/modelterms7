# Whether a Multivariate Family Defines Its Mixed Block

`TRUE` when `distrib_cross_y` comes from the family itself. The
multivariate base class rejects that generic, so the function returns
`FALSE` for a family that has not overridden it.

## Usage

``` r
.random_has_cross(d)
```

## Arguments

- d:

  A distributions7 object, under any wrappers.

## Value

A single logical.

## Details

A marginal criterion reads that block to estimate the covariance of the
effects, so a family without one is rejected when the term is built. The
class that owns the method is read from its signature and compared by
name, since [`identical()`](https://rdrr.io/r/base/identical.html)
compares object identity and fails when a package's code is re-evaluated
instead of loaded.
