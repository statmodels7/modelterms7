# Lasso Penalty on a Block of Coefficients

A block of coefficients under a Laplace prior at zero: the penalty has a
kink there, so coefficients are set EXACTLY to zero and the term
selects.

## Usage

``` r
lasso(
  x,
  label = "lasso",
  standardize = FALSE,
  lambda = NULL,
  n_lambda = NULL,
  min_ratio = NULL,
  ...
)
```

## Arguments

- x:

  A one-sided formula or a numeric matrix.

- label:

  A single non-empty string prefixed to the coefficient names.

- standardize:

  A single logical: whether to penalize each coefficient on the scale of
  its own column. See the section below.

- lambda:

  The rate of the prior. One number holds it, several are the grid the
  path visits as they stand, and `NULL`, the default, has the path build
  one. Must lie in \\(0, \infty)\\.

- n_lambda:

  How many values the path visits, at least 2. `NULL`, the default,
  leaves it to the criterion, which visits 25.

- min_ratio:

  How far down the path reaches, as a fraction of the kink that empties
  the block: smaller reaches a denser fit, larger stops sooner. Must lie
  in (0, 1). NULL, the default, leaves it to the criterion. Only the
  sweep by kink size uses it.

- ...:

  Not used, and reported: an argument named after another penalty's
  hyperparameter is the mistake this catches.

## Value

An object of class
[`PenalizedTerm`](https://statmodels7.github.io/modelterms7/reference/PenalizedTerm.md)
(a specification; see
[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

## Details

\$\$\rho(\beta) = \lambda\lVert\beta\rVert_1 -
p\log\\\left(\frac{\lambda}{2}\right).\$\$ The constant is kept, so this
is minus the log density of a Laplace at zero with rate \\\lambda\\, and
a larger \\\lambda\\ shrinks harder and keeps fewer coefficients.

The kink is at zero, so the block is fitted by a proximal method or by a
coordinate descent with the other terms held, and \\\lambda\\ is chosen
by a PATH over its own values – `bic()` by default, or `aic()` or `cv()`
– because a marginal criterion is a Laplace expansion at a mode that
sits on the kink.

**Hyperparameter.** `lambda`, admissible on \\(0, \infty)\\, swept over
`n_lambda` values from the one that empties the block down to
`min_ratio` of it.

## References

Tibshirani, R. (1996). Regression shrinkage and selection via the lasso.
*Journal of the Royal Statistical Society, Series B* 58, 267–288.

## See also

[`penalized_terms`](https://statmodels7.github.io/modelterms7/reference/penalized_terms.md)
for what the five share,
[`ridge`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`enet`](https://statmodels7.github.io/modelterms7/reference/enet.md),
[`scad`](https://statmodels7.github.io/modelterms7/reference/scad.md),
[`mcp`](https://statmodels7.github.io/modelterms7/reference/mcp.md),
[`lasso_penalty`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.html)

## Examples

``` r
dd <- data.frame(x1 = rnorm(8), x2 = rnorm(8))
built <- term_build(lasso(~ x1 + x2), dd)
term_penalty(built)@params
#> [1] "lambda"
term_smooth(built)
#> [1] FALSE

# a finer path for a wide block
term_grid(lasso(~ x1 + x2, n_lambda = 60))
#> [[1]]
#> [[1]]$lambda
#> [1] 60
#> 
#> 
```
