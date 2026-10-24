# The Family Under Any Wrappers

Follows `parent_distrib` to the distribution that a wrapper wraps, so
that a property of the family is read from the family.

## Usage

``` r
.random_family(d)
```

## Arguments

- d:

  A distributions7 object.

## Value

A distributions7 object.

## Details

The property is detected with
[`S7::prop_names()`](https://rdrr.io/pkg/S7/man/prop_names.html) instead
of testing the class, so a wrapper written later is followed without an
edit here.
