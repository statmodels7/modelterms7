# Coefficient Names of a Built Term

The names of a built additive term's coefficients, one per column of its
block and in the block's own order. They are the names a fit reports, so
a coefficient table is readable without knowing which term produced
which row.

## Usage

``` r
term_coef_names(term, ...)
```

## Arguments

- term:

  A built additive term (see
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).
  A specification throws
  `"the term has not been built; call term_build(term, data) first."`.

- ...:

  Passed to methods. No shipped method reads anything here.

## Value

A character vector of length
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md),
in column order.

## Details

A term with a non-empty `label` prefixes every name with it and a dot,
so `linpar(~ x, label = "lin")` gives `lin.(Intercept)` and `lin.x`. The
five penalized constructors set the label to the constructor's own name
by default, which is why a ridge over `x` and `g` reads `ridge.x`,
`ridge.ga` and so on;
[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
sets none, so its names are
[`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)'s
unchanged.

The names are assigned at build time and recorded, so they are also what
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
labels its columns with at other rows. Uniqueness is not enforced here;
[`check_term()`](https://statmodels7.github.io/modelterms7/reference/check_term.md)
checks it.

## See also

[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
for the count,
[`term_matrix()`](https://statmodels7.github.io/modelterms7/reference/term_matrix.md)
for the block they name,
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
for a structural term's parameter names instead.

## Examples

``` r
d <- data.frame(x = 1:4, g = factor(c("a", "b", "a", "b")))

# linpar() takes model.matrix()'s names as they come.
term_coef_names(term_build(linpar(~ x + g), d))
#> [1] "(Intercept)" "x"           "gb"         

# A label prefixes every one of them.
term_coef_names(term_build(linpar(~ x, label = "lin"), d))
#> [1] "lin.(Intercept)" "lin.x"          

# The penalized constructors label themselves by default.
term_coef_names(term_build(ridge(~ x + g), d))
#> [1] "ridge.x"  "ridge.ga" "ridge.gb"

# They are the column names of the block, and of a prediction.
b <- term_build(ridge(~ x + g), d)
identical(term_coef_names(b), colnames(term_matrix(b)))
#> [1] TRUE
```
