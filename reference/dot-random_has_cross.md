# Whether a Multivariate Family Answers Its Mixed Block

`TRUE` when `distrib_cross_y` comes from the family rather than from the
multivariate base class, whose method rejects.

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
effects, so a family without one can be fitted at held hyperparameters
and not at estimated ones. The owning class of a method is read through
its signature and compared BY NAME, never with
[`identical()`](https://rdrr.io/r/base/identical.html), which is object
identity and fails whenever a package's code is re-evaluated rather than
loaded.
