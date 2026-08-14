# Whether a Multivariate Family Can Carry Correlated Effects

`TRUE` for a family with a location block as long as its dimension and a
matrix parameter, which is what a centered prior on \\R^d\\ needs.

## Usage

``` r
.random_mv_ok(d)
```

## Arguments

- d:

  A distributions7 object, under any wrappers.

## Value

A single logical.

## Details

The question is a PROPERTY and not a list of admitted names, so a
multivariate family added later is covered without an edit here. It is
read off `params_interpretation`, the same declaration a data-based
starting value is built from. It excludes the simplex-valued families,
whose mean coordinates are one fewer than the dimension and which carry
no matrix parameter at all.
