# Build a Parametric Block

Builds the model matrix of a
[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
term's formula against `data`, prefixes the coefficient names with the
term's label, and records the blueprint
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
will reapply: the terms object with the response deleted, the factor
levels, the contrasts used, and the storage settled on.

## Arguments

- term:

  An unbuilt or built
  [`LinparTerm()`](https://statmodels7.github.io/modelterms7/reference/LinparTerm.md).

- data:

  A data frame carrying every variable the formula names.

- ...:

  Unused.

## Value

The term with `X`, `coef_names` and `blueprint` filled, so that
[`term_is_built()`](https://statmodels7.github.io/modelterms7/reference/term_is_built.md)
is `TRUE`. The block has `nrow(data)` rows, and its column names are the
coefficient names.

## The model frame

[`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html) is
called with `na.action = na.pass`, so a row with a missing covariate
keeps its place and the block stays aligned with the response, and with
`drop.unused.levels = FALSE`, so a factor level present in the data but
used by no row still gets a column. Both choices are about alignment: a
block whose rows have been silently dropped no longer matches the
response it is fitted against.

## The storage

Where the constructor left `sparse = NULL`,
[`.resolve_sparse()`](https://statmodels7.github.io/modelterms7/reference/dot-resolve_sparse.md)
settles it from `n` times the indicator column count. The value settled
on goes into `blueprint$sparse` and the `sparse` property is left as the
caller wrote it, because new rows may be any number and a prediction
deciding again could build a block of a different kind from the fitted
one.

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
names(b@blueprint)
#> [1] "terms"     "xlev"      "contrasts" "sparse"   

# A missing value keeps its row.
term_matrix(term_build(linpar(~ x), data.frame(x = c(1, NA, 3))))
#>   (Intercept)  x
#> 1           1  1
#> 2           1 NA
#> 3           1  3
```
