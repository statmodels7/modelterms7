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
its columns are dense whatever the storage, so they say nothing about
the choice. Character and logical columns are counted by their number of
distinct values, since
[`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
will code them as factors.

The count is an **upper bound**: contrasts drop one level per factor,
and an interaction of two factors is counted at the product of their
full level counts. That is the right side to err on, the number being
read against a threshold below which the sparse route loses little and
above which it wins by orders.

## See also

[`.resolve_sparse()`](https://statmodels7.github.io/modelterms7/reference/dot-resolve_sparse.md),
which reads it;
[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md).
