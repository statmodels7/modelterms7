# A Penalized Block at New Rows

Reproduces a built penalized block at `newdata`. A formula input is
reapplied through the recorded terms, levels, contrasts and storage,
exactly as
[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)'s
is. A matrix input has its expression re-evaluated **in `newdata`
alone**, and the result is checked to have the right number of rows and
columns before it is returned.

## Arguments

- term:

  A built
  [`PenalizedTerm()`](https://statmodels7.github.io/modelterms7/reference/PenalizedTerm.md).
  An unbuilt one throws
  `"the term has not been built; call term_build(term, data) first."`.

- newdata:

  A data frame. For a matrix input it must carry the matrix as a column
  of the same name the term was written with.

- ...:

  Unused.

## Value

A block of `nrow(newdata)` rows and
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
columns, in the storage the build used, with the term's coefficient
names as column names and no row names.

## Details

The matrix branch evaluates `blueprint$expr` with `newdata` as the data
and [`baseenv()`](https://rdrr.io/r/base/environment.html) as the
enclosure, so the calling environment is not on the search path. A
free-standing matrix therefore cannot be found, and the error names the
expression and says to supply it as a column of the data. That is
deliberate: resolving the name outside `newdata` would silently return
the build-time rows whenever the row counts happened to agree, which is
a wrong answer rather than an error.

Two shape checks follow the evaluation, because an expression that
resolves in `newdata` may still give the wrong block: a row count
differing from `nrow(newdata)` and a column count differing from the
build's each throw with both numbers.

## See also

[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
for the generic and the identity it satisfies,
[`penalized_terms()`](https://statmodels7.github.io/modelterms7/reference/penalized_terms.md)
for why a matrix input belongs in the data frame.

## Examples

``` r
set.seed(3)
dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20))

# A formula input reapplies its blueprint.
b <- term_build(lasso(~ x1 + x2), dd)
all.equal(term_predict(b, dd), term_matrix(b))
#> [1] TRUE
dim(term_predict(b, dd[1:5, ]))
#> [1] 5 2

# A matrix column of the data predicts on a subset of those rows.
R <- matrix(rnorm(60), 20, 3)
dd$R <- R
bm <- term_build(ridge(R), dd)
dim(term_predict(bm, dd[1:5, ]))
#> [1] 5 3

# A free-standing matrix builds and cannot be predicted from.
Rfree <- matrix(rnorm(60), 20, 3)
bf <- term_build(ridge(Rfree), dd)
try(term_predict(bf, dd[1:5, ]))
#> Error : the matrix expression `Rfree` could not be evaluated in 'newdata' (object 'Rfree' not found); supply the matrix as a column of the data.
```
