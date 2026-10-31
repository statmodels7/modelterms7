# Jeffreys Prior on a Block of Coefficients

Specifies a block of coefficients under the prior \\\|\tilde
I(\beta)\|^{a}\\, where \\\tilde I\\ is the expected Fisher information
of the block's coefficients. With \\a = 1/2\\ this is the Jeffreys prior
on the block, and over every coefficient of a binomial regression it is
the penalty of Firth (1993).

## Usage

``` r
jfp(x, label = "jfp", power = 0.5, intercept = NULL, sparse = NULL, ...)
```

## Arguments

- x:

  A one-sided formula, such as `~ x1 + x2` or `~ 0 + g`, or a numeric
  matrix, ideally a matrix column of the model data frame. A two-sided
  formula is rejected with an error.

- label:

  A single non-empty character string prefixed to the coefficient names
  as `label.name`, `"jfp"` by default.

- power:

  The power \\a\\ of the determinant, a single positive number; `0.5` by
  default.

- intercept:

  `NULL` (the default), `TRUE` or `FALSE`: whether the constant
  direction of the formula level enters the determinant, by the rule
  above.

- sparse:

  Governs the formula route: `TRUE` builds a `dgCMatrix` through
  [`Matrix::sparse.model.matrix()`](https://rdrr.io/pkg/Matrix/man/sparse.model.matrix.html),
  `FALSE` a dense model matrix, and `NULL`, the default, settles it at
  build from the size of the design. A matrix input needs no such
  argument. See the section above.

- ...:

  Not used. Any argument supplied here is rejected with an error that
  names it and lists the penalty's hyperparameters. The typical case is
  an argument named after another penalty's hyperparameter:
  `scad(~ x, gamma = 1)` signals `"'scad' has no argument 'gamma'."`.

## Value

An unbuilt
[`JeffreysTerm()`](https://statmodels7.github.io/modelterms7/reference/JeffreysTerm.md).
The build attaches a
[`JeffreysPenalty()`](https://statmodels7.github.io/modelterms7/reference/JeffreysPenalty.md)
over the block's columns.

## Details

The penalty added to the negative log-likelihood is \$\$\rho(\beta) =
-a\log\|\tilde I(\beta)\|,\qquad \tilde I = \sum_t
w_t\sum\_{q,r}\mathcal I\_{qr}(\eta_t)\\ \tilde v\_{q,t}^\top\tilde
v\_{r,t},\$\$ with \\\mathcal I\_{qr}\\ the expected information of the
family on the link scale and \\\tilde v\_{q,t}\\ the derivative of the
predictor of equation \\q\\ with respect to the block's coefficients,
preceded by the constant direction when it is included. It depends on
the family and on every predictor of the model, and the fitting layer
evaluates it together with its derivatives.

The power is fixed. There is no hyperparameter, so nothing is estimated
for this term and it adds no row to an outer criterion.

The constant direction of a formula level is the derivative of the
predictor along a constant shift of that level's linear predictor: the
column of ones for an equation. With `intercept = NULL` it is included
when the column of ones lies in the span of the level's unpenalized
columns and not in the span of the block, so the penalty is invariant to
centering a column of the block whenever the model is. `TRUE` includes
it unless the block already spans the column of ones, and `FALSE` never
includes it.

The penalty moves an estimate as \\a\\\tilde p\\ pseudo-observations
would, with \\\tilde p\\ the size of \\\tilde I\\: with a Poisson log
link and an intercept in the model, \\\sum_t(\hat\mu_t - y_t) = a\tilde
p\\. On a dispersion or shape parameter the effect is strong: the
degrees of freedom of a Student t fitted to Gaussian data stay finite.

A formula is coded with its contrasts, as
[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
codes it, and the intercept column is then removed, so `jfp(~ g)` over a
factor of \\k\\ levels gives the \\k - 1\\ contrasts and `jfp(~ 0 + g)`
all \\k\\ indicators. This differs from
[`ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
whose formula gives every indicator: the Jeffreys penalty depends on the
predictor alone, so it does not identify the direction in which a full
set of indicators and the equation's intercept trade a constant. A block
whose columns are linearly dependent is rejected when the term is built.

## References

Firth, D. (1993). Bias reduction of maximum likelihood estimates.
*Biometrika* 80, 27–38.

Kosmidis, I. and Firth, D. (2021). Jeffreys-prior penalty, finiteness
and shrinkage in binomial-response generalized linear models.
*Biometrika* 108, 71–82.

## See also

[`JeffreysPenalty()`](https://statmodels7.github.io/modelterms7/reference/JeffreysPenalty.md),
[`ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.md)
for the same input handling.

## Examples

``` r
set.seed(3)
dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
b <- term_build(jfp(~ x1 + x2), dd)
b@coef_names
#> [1] "jfp.x1" "jfp.x2"
b@penalty@power
#> [1] 0.5

# one entry, with no hyperparameter
e <- term_penalties(b)
e[[1]]$penalty@params
#> character(0)
```
