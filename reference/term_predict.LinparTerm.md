# A Parametric Block at New Rows

Rebuilds the model frame at `newdata` against the levels recorded at
build time, then the model matrix with the recorded contrasts and in the
recorded storage, and labels the columns with the term's coefficient
names. Nothing is re-derived from the new rows.

## Arguments

- term:

  A built
  [`LinparTerm()`](https://statmodels7.github.io/modelterms7/reference/LinparTerm.md).
  An unbuilt one throws
  `"the term has not been built; call term_build(term, data) first."`.

- newdata:

  A data frame carrying every variable the formula names.

- ...:

  Unused.

## Value

A block of `nrow(newdata)` rows and
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
columns, in the storage the build settled on, with the term's
coefficient names as column names.

## Details

The levels come from `blueprint$xlev`, so a factor in `newdata` need
carry only the levels its own rows use and still gets the full set of
columns. A level the blueprint does not know is rejected by
[`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html) with
`"factor g has new levels zz"`, which is the right answer: a coefficient
was never fitted for it.

The storage comes from `blueprint$sparse` rather than being decided
again, so a prediction does not spend at new rows what the build was
careful not to. `na.action = na.pass` again keeps every row.

## See also

[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
for the generic and the identity it satisfies,
[`term_build.LinparTerm()`](https://statmodels7.github.io/modelterms7/reference/term_build.LinparTerm.md)
for what recorded the blueprint.

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

# A level the fit never saw is refused.
bad <- dd
levels(bad$g) <- c("a", "b", "c", "zz")
try(term_predict(b, bad))
#> Error in model.frame.default(bp$terms, newdata, na.action = stats::na.pass,  : 
#>   factor g has new levels zz
```
