# One, Two or No Directions

Reads the `direction` argument of
[`term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.md),
[`term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.md)
and
[`term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.md)
into a list: empty for the second order, one element for the third, two
for the fourth.

## Usage

``` r
.gas_dirs(direction)
```

## Arguments

- direction:

  `NULL`, a numeric vector, or a list of one or two numeric vectors.

## Value

A list of zero, one or two numeric vectors.

## Details

A single numeric vector is one direction, a list is as many as it holds,
and `NULL` is none. The three orders then share one recursion, and a
caller writes `direction = v` at the third order and
`direction = list(v, w)` at the fourth.
