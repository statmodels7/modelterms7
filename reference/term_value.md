# The Contribution of a Term at Its Current Coefficients

The values a term contributes to the linear predictor. For a linear term
this is the block times the coefficients and carries no information the
block does not. For a nonlinear one it is \\f(x;\theta)\\, which the
Jacobian alone does not give and which a Gauss-Newton step needs beside
it.

## Usage

``` r
term_value(term, coef = NULL, newdata = NULL, ...)
```

## Arguments

- term:

  A built term.

- coef:

  The coefficients. Optional for a nonlinear term, which carries the
  ones it was last refreshed at.

- newdata:

  An optional data frame; the contribution is returned on its rows
  instead of on the ones the term was built from.

- ...:

  Passed to methods.

## Value

A numeric vector of one value per observation: `nrow(newdata)` of them
where `newdata` was given, and as many as the term was built on
otherwise.

## Details

`newdata` asks for the same contribution on other rows, and is what a
predictor needs where the block is a Jacobian: there
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
times the coefficients is the linearization and not the contribution,
and the two differ by whatever the linearization drops. For
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
that difference is a step at the break-point in a construction that is
continuous. Rows arriving here are treated as
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
treats them, through the levels and constants the blueprint recorded,
never rebuilt.

## See also

[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
for the block at the same coefficients,
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
for the block at other rows,
[`seg_psi()`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md)
for a break-point read off the same coefficients.

## Examples

``` r
dd <- data.frame(x = seq(0, 2, length.out = 20))
built <- term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1)), dd)
head(term_value(built), 3)
#> [1] 2.000000 1.800175 1.620315
head(term_value(built, newdata = data.frame(x = c(0, 1, 2))), 3)
#> [1] 2.0000000 0.7357589 0.2706706
```
