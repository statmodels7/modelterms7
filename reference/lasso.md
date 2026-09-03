# Lasso Penalty on a Block of Coefficients

A block of coefficients under a Laplace prior at zero: the penalty has a
kink there, so coefficients are set exactly to zero and the term
selects.

## Usage

``` r
lasso(
  x,
  label = "lasso",
  standardize = FALSE,
  lambda = NULL,
  id = NULL,
  n_lambda = 25,
  min_ratio = 1e-04,
  sparse = NULL,
  ...
)
```

## Arguments

- x:

  A one-sided formula, such as `~ x1 + x2` or `~ 0 + g`, or a numeric
  matrix, ideally a matrix column of the model data frame. A two-sided
  formula throws.

- label:

  A single non-empty character string prefixed to the coefficient names
  as `label.name`, `"lasso"` by default, so a block over `x1` reads
  `lasso.x1`. Two penalized terms in one formula stay apart by their
  labels.

- standardize:

  A single logical, `FALSE` by default: whether to penalize each
  coefficient on the scale of its own column. See the section above.

- lambda:

  The rate of the prior. One number holds it, several are the grid the
  path visits as they stand, and `NULL`, the default, has the path build
  one. Must lie in \\(0, \infty)\\.

- id:

  A label sharing this term's hyperparameter with those of other terms
  carrying the same one: they are then estimated at a single value,
  wherever in the formula they sit. `NULL`, the default, shares nothing.
  See
  [`term_ids()`](https://statmodels7.github.io/modelterms7/reference/term_ids.md),
  whose page says what the labels mean and what care they want.

- n_lambda:

  How many values the path visits, a whole number of at least 2, `25` by
  default. The axis descends four decades of kink size, and that many
  points are what covers it.

- min_ratio:

  How far down the path reaches, as a fraction of the kink that empties
  the block: smaller reaches a denser fit, larger stops sooner. A single
  number in \\(0, 1)\\, `1e-4` by default. Only the sweep by kink size
  reads it.

- sparse:

  Governs the formula route: `TRUE` builds a `dgCMatrix` through
  [`Matrix::sparse.model.matrix()`](https://rdrr.io/pkg/Matrix/man/sparse.model.matrix.html),
  `FALSE` a dense model matrix, and `NULL`, the default, settles it at
  build from the size of the design. A matrix input needs no such
  argument. See the section above.

- ...:

  Not used, and reported. An argument named after another penalty's
  hyperparameter is the mistake this catches: `scad(~ x, gamma = 1)`
  throws `"'scad' has no argument 'gamma'."` and lists the ones it does
  have.

## Value

An unbuilt
[`PenalizedTerm()`](https://statmodels7.github.io/modelterms7/reference/PenalizedTerm.md):
a specification, with `X`, `coef_names`, `blueprint` and `penalty` empty
until
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
fills them, and the penalty attached there over as many coefficients as
the block turns out to have.

## Details

\$\$\rho(\beta) = \lambda\lVert\beta\rVert_1 -
p\log\\\left(\frac{\lambda}{2}\right).\$\$ The constant is kept, so this
is minus the log density of a Laplace at zero with rate \\\lambda\\, and
a larger \\\lambda\\ shrinks harder and keeps fewer coefficients.

The kink is at zero, so the block is fitted by a proximal method or by a
coordinate descent with the other terms held, and \\\lambda\\ is chosen
by a path over its own values, scored by
[`statmodels7::bic()`](https://statmodels7.github.io/statmodels7/reference/aic.html)
by default or by
[`statmodels7::aic()`](https://statmodels7.github.io/statmodels7/reference/aic.html)
or
[`statmodels7::cv()`](https://statmodels7.github.io/statmodels7/reference/cv.html).
A marginal criterion cannot be used: it is a Laplace expansion at a mode
that sits on the kink.

**Hyperparameter.** `lambda`, admissible on \\(0, \infty)\\, swept over
`n_lambda` values from the one that empties the block down to
`min_ratio` of it.

## References

Tibshirani, R. (1996). Regression shrinkage and selection via the lasso.
*Journal of the Royal Statistical Society, Series B* 58, 267–288.

## See also

[`penalized_terms()`](https://statmodels7.github.io/modelterms7/reference/penalized_terms.md)
for what the five share,
[`ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`enet()`](https://statmodels7.github.io/modelterms7/reference/enet.md),
[`scad()`](https://statmodels7.github.io/modelterms7/reference/scad.md),
[`mcp()`](https://statmodels7.github.io/modelterms7/reference/mcp.md),
[`penalties7::lasso_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.html)

## Examples

``` r
set.seed(3)
dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
b <- term_build(lasso(~ x1 + x2), dd)
p <- term_penalty(b)

# The kink is at zero, so the block is not smooth and needs a path.
c(smooth = term_smooth(b))
#> smooth 
#>  FALSE 
penalties7::penalty_kinks(p, list(lambda = 1))
#> [1] 0

# The value is minus a Laplace log-density at rate lambda.
beta <- c(0.4, -1.1)
all.equal(penalties7::penalty_value(p, beta, list(lambda = 2)),
          2 * sum(abs(beta)) - 2 * log(2 / 2))
#> [1] TRUE

# The path's length and depth, at the defaults and set.
term_grid(lasso(~ x1 + x2))
#> [[1]]
#> [[1]]$lambda
#> [1] 25
#> 
#> 
term_grid(lasso(~ x1 + x2, n_lambda = 60))
#> [[1]]
#> [[1]]$lambda
#> [1] 60
#> 
#> 
term_path_min(lasso(~ x1 + x2, min_ratio = 1e-3))
#> [[1]]
#> [1] 0.001
#> 


# Fitted. The data are simulated from a known truth, so the
# estimates below can be read against it.
if (requireNamespace("statmodels7", quietly = TRUE)) {
  set.seed(11)
  XX <- matrix(rnorm(150 * 8), 150, 8)
  fd <- data.frame(y = as.numeric(XX %*% c(2, -1.5, 1, rep(0, 5))) +
                     rnorm(150, sd = 0.4))
  fd$X <- XX
  cf <- coef(statmodels7::statmod(y ~ lasso(X),
                                  distributions7::gaussian1_distrib(), fd))$mu
  # truth: the first three columns carry 2, -1.5 and 1, the other five
  # nothing. The hyperparameter is chosen by BIC.
  c(round(cf[2:5], 2), kept = sum(cf[-1] != 0))
}
#> lasso.1 lasso.2 lasso.3 lasso.4    kept 
#>    1.97   -1.51    0.99    0.00    5.00 
```
