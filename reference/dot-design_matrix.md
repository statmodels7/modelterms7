# The Model Matrix of a Formula, in Either Storage

[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html) or
[`sparse.model.matrix`](https://rdrr.io/pkg/Matrix/man/sparse.model.matrix.html)
on the same terms object, with the bookkeeping stripped either way.

## Usage

``` r
.design_matrix(tt, mf, contrasts = NULL, sparse = FALSE)
```

## Arguments

- tt:

  A terms object.

- mf:

  The model frame.

- contrasts:

  The contrasts, or `NULL` for the session's.

- sparse:

  Whether to build a `dgCMatrix`.

## Value

A numeric matrix or a `dgCMatrix`.

## Details

The sparse route BUILDS the matrix sparse; it does not build a dense one
and compress it, which would cost the memory the choice exists to avoid.
Measured at 20000 rows and a factor of 1000 levels, 0.002 s and 1.8 MB
against
[`stats::model.matrix`](https://rdrr.io/r/stats/model.matrix.html)'s
0.100 s and 161.5 MB, the numbers identical; and a design that would be
32 GB dense builds in 0.02 s and 19 MB, which is what says there is no
dense intermediate.

It is worth it where the formula carries a factor of MANY LEVELS, whose
indicator columns are one non-zero per row. On a formula of numeric
covariates the block is dense whatever is asked for, and the sparse
storage then costs more than it saves.

## See also

[`linpar`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
