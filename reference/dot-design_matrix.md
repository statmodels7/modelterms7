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

  The model frame from which `tt` was built, or one built against the
  same levels.

- contrasts:

  A named list of contrasts, or `NULL` for the session's.

- sparse:

  `TRUE` for a `dgCMatrix`, `FALSE` for a base matrix.

## Value

A list of two: `X`, the block, and `contrasts`, the named list that
[`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
recorded. On the dense route it is `NULL` if
[`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
coded no factor, and on the sparse route it is an empty named list.

## Details

The two routes give identical numbers and differ in the storage and in
the cost of producing it. The sparse route builds the matrix sparse and
never forms a dense intermediate.

The contrasts are returned beside the block and are not kept on it. A
consumer reads the block, and the blueprint records the contrasts.

## See also

[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md),
whose build and prediction both call this;
[`.resolve_sparse()`](https://statmodels7.github.io/modelterms7/reference/dot-resolve_sparse.md)
for the choice of storage.
