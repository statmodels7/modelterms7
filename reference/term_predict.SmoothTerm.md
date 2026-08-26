# A Smooth Term's Block at New Rows

Evaluates a built smooth's recorded basis at the covariates in `newdata`
and applies the transforms the build computed: the Demmler-Reinsch basis
and its centering and scaling for
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md), the
centered tensor basis for
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md),
then the `by` variable against the recorded levels. Nothing is derived
from the new rows.

## Arguments

- term:

  A built
  [`SmoothTerm()`](https://statmodels7.github.io/modelterms7/reference/SmoothTerm.md).
  An unbuilt one throws
  `"the term has not been built; call term_build(term, data) first."`.

- newdata:

  A data frame carrying the covariates and the `by` variable.

- ...:

  Unused.

## Value

A block of `nrow(newdata)` rows and
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
columns, in the storage the build settled on, with the term's
coefficient names as column names and no row names.

## Details

The basis object in `blueprint$core` carries the knots and the
transform, so the columns at new rows are the same functions of the
covariate as the fitted ones, so \\\tilde{X}\beta\\ is the fitted smooth
evaluated there. Rebuilding instead would place the knots on the new
range and compute another transform: measured on 80 points, predicting
on the first ten agrees with those rows of the block exactly and
rebuilding differs by 2.85.

A factor `by` is expanded against `blueprint$by_levels`, so `newdata`
need carry only the levels its own rows use and still gets every column,
and the block is built in the storage the build settled on.

New covariate values outside the range the basis was placed on are
evaluated, not refused. A B-spline is zero beyond its knots, so the
fitted function flattens rather than extrapolating a trend; read a
prediction far outside the fitting range with that in mind.

## See also

[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
for the generic and the identity it satisfies,
[`term_build.SmoothTerm()`](https://statmodels7.github.io/modelterms7/reference/term_build.SmoothTerm.md)
for what recorded the transform.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(80)), g = factor(rep(letters[1:4], 20)))
b <- term_build(s(x, k = 8), dd)
X <- term_matrix(b)

# Reapplying is exact; rebuilding on the same rows is a different basis.
max(abs(term_predict(b, dd[1:10, ]) - X[1:10, ]))
#> [1] 0
max(abs(term_matrix(term_build(s(x, k = 8), dd[1:10, ])) - X[1:10, ]))
#> [1] 3.628495

# A factor `by` keeps every level's columns at a subset that has two.
bf <- term_build(s(x, k = 5, by = g), dd)
nd <- droplevels(dd[dd$g %in% c("a", "b"), ])
c(levels_here = nlevels(nd$g), cols = ncol(term_predict(bf, nd)))
#> levels_here        cols 
#>           2          16 
```
