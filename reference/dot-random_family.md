# The Family Under Any Wrappers

Follows `parent_distrib` to the distribution a wrapper wraps, so that a
question about the FAMILY is asked of the family.

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

The property is asked for with
[`S7::prop_names()`](https://rconsortium.github.io/S7/reference/prop_names.html)
rather than the class being tested, which is what lets a wrapper written
later be followed without an edit here.
