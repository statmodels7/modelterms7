# The Directions of a Derivative Contraction

Reads the `direction` argument of
[`.gas_curvature_core()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_curvature_core.md)
into a list: empty for the second order, one element for the third, two
for the fourth.
[`term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.md)
sets it to `NULL`,
[`term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.md)
to its `direction` and
[`term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.md)
to the list of its `directions`.

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
and `NULL` is none. The three orders then share one recursion. A caller
writes `direction = v` in
[`term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.md)
and `directions = list(v, w)` in
[`term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.md).
