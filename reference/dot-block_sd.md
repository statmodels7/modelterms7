# Column Standard Deviations of a Block, in Either Storage

The standard deviation of each column of a design block, computed
without densifying a Matrix. It is the spread `standardize = TRUE`
divides each coefficient by.

## Usage

``` r
.block_sd(X)
```

## Arguments

- X:

  A design block: a numeric matrix or a two-dimensional Matrix.

## Value

A numeric vector of length `ncol(X)`, every entry finite and strictly
positive.

## Details

On a base matrix it is [`stats::sd()`](https://rdrr.io/r/stats/sd.html)
column by column. On a Matrix it is assembled from the column sums and
sums of squares, \\s_j^2 = (\sum_i x\_{ij}^2 - n\bar{x}\_j^2)/(n-1)\\,
so a sparse block is never expanded. The sum of squares is floored at
zero before the square root, that expression being able to go slightly
negative by rounding on a nearly constant column.

A column with no spread takes \\s_j = 1\\, which penalizes its
coefficient on its own scale instead of dividing by zero. That covers a
constant column, a single-row block and any column whose spread comes
back non-finite.

## See also

[`penalized_terms()`](https://statmodels7.github.io/modelterms7/reference/penalized_terms.md)
for the standardization this feeds,
[`penalties7::map_diagonal()`](https://statmodels7.github.io/penalties7/reference/map_diagonal.html)
for the map it becomes.
