# A Break-Point Term's Block at New Rows

The working block evaluated at `newdata`, at the break-points and the
coefficients the term currently carries. The positions are not
re-derived from the new rows, and any development of a coefficient is
reapplied through its sub-terms' own blueprints.

## Arguments

- term:

  A built
  [`SegTerm()`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md).
  An unbuilt one throws
  `"the term has not been built; call term_build(term, data) first."`.

- newdata:

  A data frame carrying the covariate and whatever the subformulas name.

- ...:

  Unused.

## Value

A block of `nrow(newdata)` rows and
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
columns, with the term's coefficient names as column names.

## Details

Re-deriving the break-points would give a different model at every set
of rows, which is the general reason
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
reapplies instead of rebuilding. For a break-point term it is sharper
than usual: the positions are what the fit estimated, and a subset of
the data may not even contain one.

[`term_value()`](https://statmodels7.github.io/modelterms7/reference/term_value.md)
with `newdata` is what gives the contribution there. For `seg` the two
differ, the block being a Jacobian.

## See also

[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
for the generic,
[`term_value()`](https://statmodels7.github.io/modelterms7/reference/term_value.md)
for the contribution,
[`seg_psi()`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md)
for the positions it is evaluated at.

## Examples

``` r
set.seed(1)
d <- data.frame(x = sort(runif(120, 0, 10)))
d$y <- 1 + 0.5 * d$x + 2 * pmax(d$x - 6, 0) + rnorm(120, sd = 0.4)
b <- term_refresh(term_build(seg(x, npsi = 1), d), c(0.5, 2, 6))

# On the fitting data it is the block itself.
max(abs(term_predict(b, d) - term_matrix(b)))
#> [1] 0

# At other rows the break-point is the fitted one, not one re-derived
# from those rows: here every row is below it.
nd <- data.frame(x = c(1, 2, 3))
term_predict(b, nd)
#>      seg.beta seg.gamma1 seg.psi1
#> [1,]        1          0        0
#> [2,]        2          0        0
#> [3,]        3          0        0
seg_psi(b)
#> [1] 6
```
