# The Contribution of a Term at Its Current Coefficients

The values a term contributes to the linear predictor. For a linear term
this is the block times the coefficients and carries no information the
block does not; for a nonlinear one it is \\f(x;\theta)\\, which the
Jacobian alone does not give, and which a Gauss-Newton step needs beside
it.

## Usage

``` r
term_value(term, coef = NULL, ...)
```

## Arguments

- term:

  A built term.

- coef:

  The coefficients. Optional for a nonlinear term, which carries the
  ones it was last refreshed at.

- ...:

  Passed to methods.

## Value

A numeric vector, one value per observation.

## See also

[`term_refresh`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)

## Examples

``` r
dd <- data.frame(x = seq(0, 2, length.out = 20))
built <- term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1)), dd)
head(term_value(built), 3)
#> [1] 2.000000 1.800175 1.620315
```
