# Unpenalized Parametric Term

Creates the specification of an unpenalized parametric block: the design
matrix of a one-sided formula, with the usual
[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html) conventions
for factors, contrasts, interactions and the intercept.

## Usage

``` r
linpar(formula, label = "")
```

## Arguments

- formula:

  A one-sided formula, e.g. `~ x1 + x2`.

- label:

  A character string; when non-empty it is prefixed to the coefficient
  names as `label.name`.

## Value

An object of class
[`LinparTerm`](https://statmodels7.github.io/modelterms7/reference/LinparTerm.md)
(a specification; see
[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

## Details

The block is the model matrix \\X\\ of the formula and the term's
contribution to the predictor is linear in its coefficients,

\$\$\eta = X\beta,\$\$

with no penalty attached, so all \\p = \operatorname{ncol}(X)\\
coefficients are free and
[`edf`](https://statmodels7.github.io/modelterms7/reference/edf.md)
counts every one of them.

[`interpret_formula`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
collects the bare covariates of a model formula into one term of this
kind, so `y ~ x1 + x2` and `y ~ linpar(~ x1 + x2)` produce the same
block; the explicit constructor exists for callers who want several
parametric blocks with distinct labels.

Building the term records a blueprint: the terms object, the factor
levels and the contrasts.
[`term_predict`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
reapplies the mapping through that blueprint, so a factor column in new
data is encoded against the levels seen at build time, and a level the
blueprint does not know is rejected rather than re-encoded. Missing
values are propagated (`na.pass`), never dropped, so the block stays
row-aligned with the response.

## See also

[`ridge`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`lasso`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`scad`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`mcp`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`enet`](https://statmodels7.github.io/modelterms7/reference/ridge.md)

## Examples

``` r
dd <- data.frame(x = 1:4, g = factor(c("a", "a", "b", "b")))
built <- term_build(linpar(~ x + g), dd)
term_matrix(built)
#>   (Intercept) x gb
#> 1           1 1  0
#> 2           1 2  0
#> 3           1 3  1
#> 4           1 4  1
term_coef_names(built)
#> [1] "(Intercept)" "x"           "gb"         
```
