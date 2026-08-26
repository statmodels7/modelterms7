# The Model Matrix of a Formula, in Either Storage

Runs
[`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) or
[`Matrix::sparse.model.matrix()`](https://rdrr.io/pkg/Matrix/man/sparse.model.matrix.html)
on the same terms object and model frame, strips the `assign` and
`contrasts` attributes, and returns the block beside the contrasts that
were used.

## Usage

``` r
.design_matrix(tt, mf, contrasts = NULL, sparse = FALSE)
```

## Arguments

- tt:

  A terms object, from
  [`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html) or
  [`stats::delete.response()`](https://rdrr.io/r/stats/delete.response.html).

- mf:

  The model frame `tt` was built from, or one built against the same
  levels.

- contrasts:

  A named list of contrasts, or `NULL` for the session's.

- sparse:

  `TRUE` for a `dgCMatrix`, `FALSE` for a base matrix.

## Value

A list of two: `X`, the block, and `contrasts`, the named list
[`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
recorded, which is `NULL` when no factor was coded.

## Details

The two routes give identical numbers; what differs is the storage and
the cost of producing it. The sparse route builds the matrix sparse and
never forms a dense intermediate: at 20000 rows and a factor of 1000
levels, 0.007 s and 1.8 MB against 0.164 s and 161.3 MB.

The contrasts come back beside the block instead of staying on it: the
block is what a consumer reads, and the contrasts are what the blueprint
records.

## See also

[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md),
whose build and prediction both call this;
[`.resolve_sparse()`](https://statmodels7.github.io/modelterms7/reference/dot-resolve_sparse.md)
for the choice of storage.
