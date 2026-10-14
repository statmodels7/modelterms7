# Settle Whether a Block Is Built Sparse

Passes an explicit `TRUE` or `FALSE` through, and where the caller left
`NULL` decides from the size of the indicator part: sparse when
`n * ncol_ind` exceeds \\10^5\\.

## Usage

``` r
.resolve_sparse(sparse, n, ncol_ind)
```

## Arguments

- sparse:

  `TRUE`, `FALSE`, or `NULL` to decide here.

- n:

  The number of rows.

- ncol_ind:

  The columns coming from indicators, from
  [`.indicator_cols()`](https://statmodels7.github.io/modelterms7/reference/dot-indicator_cols.md).

## Value

A single logical, never `NA`. A non-finite `n` or `ncol_ind` gives
`FALSE`.

## Details

The dense indicator part holds `n * ncol_ind` cells where the sparse one
holds one non-zero per row, so that product is what separates the two
routes, and the threshold is read off it rather than off a count of
levels.

Below about \\10^5\\ cells the sparse route is slightly slower than the
dense route, and above it the sparse route is faster, so the threshold
is placed at \\10^5\\.

A design carrying no factor has no indicator part, so the product is
zero and the block is built dense, which suits a purely continuous
design.

## See also

[`.indicator_cols()`](https://statmodels7.github.io/modelterms7/reference/dot-indicator_cols.md)
for the count it reads,
[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
for the argument it settles.
