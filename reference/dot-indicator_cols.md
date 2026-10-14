# How Many Columns of a Model Matrix Come From Factors

Counts the columns a terms object contributes through indicator
variables, read from the model frame without building the matrix. It is
the second factor of the product
[`.resolve_sparse()`](https://statmodels7.github.io/modelterms7/reference/dot-resolve_sparse.md)
compares against its threshold.

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

A single number, `0` when no term carries a factor.

## Details

A term contributes the product of its variables' level counts, a numeric
variable counting one, and a term carrying no factor at all is skipped:
its columns have few zero entries for sparse storage to exploit, so they
do not enter the choice. Character and logical columns are counted by
their number of distinct values, since
[`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
codes them as factors.

The count is an **upper bound**: contrasts drop one level per factor,
and an interaction of two factors is counted at the product of their
full level counts. The count is read against a threshold on
`n * ncol_ind`, so it overestimates the number of indicator columns of
the built block.

## See also

[`.resolve_sparse()`](https://statmodels7.github.io/modelterms7/reference/dot-resolve_sparse.md),
which reads it;
[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md).
