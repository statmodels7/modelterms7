# A Nonlinear Term's Block at New Rows

The Jacobian of \\f\\ evaluated at the new rows, at the coefficients the
term currently carries. Each parameter's subformula is reapplied through
its sub-terms' own blueprints, so a basis or a set of contrasts inside a
submodel is not relearned.

## Arguments

- term:

  A built
  [`NlTerm()`](https://statmodels7.github.io/modelterms7/reference/NlTerm.md).
  An unbuilt one throws
  `"the term has not been built; call term_build(term, data) first."`.

- newdata:

  A data frame carrying the covariates the function names and whatever
  the subformulas name.

- ...:

  Unused.

## Value

A block of `nrow(newdata)` rows and
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
columns, with the term's coefficient names as column names.

## Details

The coefficients are the ones
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
last committed, not new ones: the block a prediction returns must
multiply the same coefficients the fit reached, or \\\tilde{X}\beta\\ is
not the fitted contribution.
[`term_value()`](https://statmodels7.github.io/modelterms7/reference/term_value.md)
is what reports \\f\\ itself at new rows.

Reading the block on its own is rarely what a caller wants for a
nonlinear term. \\\tilde{X}\beta\\ is the **linearization** of \\f\\:
the two agree in the increment a Gauss-Newton step takes, never in the
value. `term_value(term, newdata = ...)` gives the contribution.

## See also

[`term_value()`](https://statmodels7.github.io/modelterms7/reference/term_value.md)
for the contribution itself,
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
for the block at other coefficients,
[`nl_fderiv()`](https://statmodels7.github.io/modelterms7/reference/nl_fderiv.md)
for the derivatives in the parameters.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = seq(0, 3, length.out = 60))
dd$y <- 2 * exp(-1.3 * dd$x) + rnorm(60, sd = 0.05)
b <- term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3)), dd)

# On the fitting data it is the block itself.
all.equal(term_predict(b, dd), term_matrix(b))
#> [1] TRUE

# At other rows it is the Jacobian there.
nd <- data.frame(x = c(0.5, 1.5, 2.5))
term_predict(b, nd)
#>            nl.a       nl.r
#> [1,] 0.52204578 -0.5220458
#> [2,] 0.14227407 -0.4268222
#> [3,] 0.03877421 -0.1938710

# Which is the linearization, not the function: read the value for that.
cbind(linear = as.numeric(term_predict(b, nd) %*% b@blueprint$coef),
      value = term_value(b, newdata = nd),
      truth = 2 * exp(-1.3 * nd$x))
#>          linear      value      truth
#> [1,]  0.3654320 1.04409155 1.04409155
#> [2,] -0.2703207 0.28454814 0.28454814
#> [3,] -0.1744839 0.07754842 0.07754842
```
