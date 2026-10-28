# A Parametric Block at New Rows

Rebuilds the model frame at `newdata` against the levels recorded at
build time, then the model matrix with the recorded contrasts and in the
recorded storage, and labels the columns with the term's coefficient
names. The terms object, the levels, the contrasts, the storage and the
intercept setting are all read from the `build_info`, and the new rows
supply only the values of the variables.

## Arguments

- term:

  A built
  [`LinparTerm()`](https://statmodels7.github.io/modelterms7/reference/LinparTerm.md).
  An unbuilt one signals the error
  `"the term has not been built; call term_build(term, data) first."`.

- newdata:

  A data frame carrying every variable that the formula names.

- ...:

  Unused.

## Value

A block of `nrow(newdata)` rows and
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
columns, in the storage chosen at build, with the term's coefficient
names as column names.

## Details

The levels come from `build_info$xlev`, so a factor in `newdata` need
carry only the levels that its own rows use and still gets the full set
of columns. A level that the `build_info` does not record is rejected by
[`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html) with
`"factor g has new levels zz"`, because the fit did not estimate a
coefficient for that level.

The storage comes from `build_info$sparse` rather than being decided
again, so a prediction builds the same kind of block as the fit and does
not densify a sparse block at new rows. `na.action = na.pass` again
keeps every row.

## See also

[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
for the generic and the identity that it satisfies,
[`term_build.LinparTerm()`](https://statmodels7.github.io/modelterms7/reference/term_build.LinparTerm.md)
for what recorded the `build_info`.

## Examples

``` r
dd <- data.frame(x = 1:8, g = factor(rep(c("a", "b", "c", "d"), 2)))
b <- term_build(linpar(~ x + g), dd)

# On the fitting data it returns the block itself.
all.equal(term_predict(b, dd), term_matrix(b))
#> [1] TRUE

# A subset that drops two levels keeps all five columns.
nd <- droplevels(dd[dd$g %in% c("a", "b"), ])
levels(nd$g)
#> [1] "a" "b"
dim(term_predict(b, nd))
#> [1] 4 5

# A level that is absent from the fitting data is rejected.
bad <- dd
levels(bad$g) <- c("a", "b", "c", "zz")
try(term_predict(b, bad))
#> Error in model.frame.default(bp$terms, newdata, na.action = stats::na.pass,  : 
#>   factor g has new levels zz
```
