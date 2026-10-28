# Build a Parametric Block

Builds the model matrix of a
[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
term's formula against `data`, prefixes the coefficient names with the
term's label, and records the `build_info` that
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
will reapply: the terms object from which the response is deleted, the
factor levels, the contrasts used, the storage chosen at build and
`drop_intercept`.

## Arguments

- term:

  An unbuilt or built
  [`LinparTerm()`](https://statmodels7.github.io/modelterms7/reference/LinparTerm.md).

- data:

  A data frame carrying every variable that the formula names.

- ...:

  Unused.

## Value

The term, in which `X`, `coef_names` and `build_info` are filled, so
that
[`term_is_built()`](https://statmodels7.github.io/modelterms7/reference/term_is_built.md)
is `TRUE`. The block has `nrow(data)` rows, and its column names are the
coefficient names.

## The model frame

[`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html) is
called with `na.action = na.pass`, so a row with a missing covariate
keeps its place and the block stays aligned with the response, and with
`drop.unused.levels = FALSE`, so a factor level present in the data but
used by no row still gets a column. The first choice keeps the rows of
the block aligned with the response, and the second keeps the columns of
the block aligned with the coefficients of the fit.

## The storage

Where the constructor left `sparse = NULL`,
[`.resolve_sparse()`](https://statmodels7.github.io/modelterms7/reference/dot-resolve_sparse.md)
settles it from `n` times the indicator column count. The value settled
on goes into `build_info$sparse` and the `sparse` property is left as
the caller wrote it, because new data may have any number of rows, and
if a prediction settled the storage again, it could build a block of a
different kind from the fitted one.

## See also

[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md),
[`term_predict.LinparTerm()`](https://statmodels7.github.io/modelterms7/reference/term_predict.LinparTerm.md),
[`.resolve_sparse()`](https://statmodels7.github.io/modelterms7/reference/dot-resolve_sparse.md).

## Examples

``` r
dd <- data.frame(x = 1:8, g = factor(rep(c("a", "b", "c", "d"), 2)))
b <- term_build(linpar(~ x + g), dd)
term_matrix(b)
#>   (Intercept) x gb gc gd
#> 1           1 1  0  0  0
#> 2           1 2  1  0  0
#> 3           1 3  0  1  0
#> 4           1 4  0  0  1
#> 5           1 5  0  0  0
#> 6           1 6  1  0  0
#> 7           1 7  0  1  0
#> 8           1 8  0  0  1
names(b@build_info)
#> [1] "terms"          "xlev"           "contrasts"      "sparse"        
#> [5] "drop_intercept"

# A missing value keeps its row.
term_matrix(term_build(linpar(~ x), data.frame(x = c(1, NA, 3))))
#>   (Intercept)  x
#> 1           1  1
#> 2           1 NA
#> 3           1  3
```
