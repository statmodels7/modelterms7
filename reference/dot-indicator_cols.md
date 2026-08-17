# How Many Columns of a Model Matrix Come From Factors

The number of columns a terms object contributes through indicator
variables, counted from the model frame without building the matrix.

## Usage

``` r
.indicator_cols(tt, mf)
```

## Arguments

- tt:

  A terms object.

- mf:

  The model frame it was built from.

## Value

A single number, zero when no term carries a factor.

## Details

A term contributes the product of its variables' level counts, a numeric
variable counting one; a term carrying no factor contributes columns
that are dense whatever the storage, and is not counted. The count is an
upper bound, contrasts dropping one level per factor, which is the right
side to err on: it is read against a threshold below which the sparse
route loses little and above which it wins by orders.

## See also

[`.resolve_sparse`](https://statmodels7.github.io/modelterms7/reference/dot-resolve_sparse.md)
