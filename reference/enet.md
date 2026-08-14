# Elastic Net Penalty on a Block of Coefficients

The lasso and the ridge mixed: a kink at zero, so the term still
selects, and a quadratic part that keeps correlated coefficients
together instead of choosing arbitrarily among them.

## Usage

``` r
enet(
  x,
  label = "enet",
  by = NULL,
  standardize = FALSE,
  lambda = NULL,
  alpha = NULL,
  n_lambda = NULL,
  n_alpha = NULL,
  min_ratio = NULL,
  ...
)
```

## Arguments

- x:

  A one-sided formula or a numeric matrix.

- label:

  A single non-empty string prefixed to the coefficient names.

- by:

  Reserved for a later release; must be `NULL`.

- standardize:

  A single logical: whether to penalize each coefficient on the scale of
  its own column. See the section below.

- lambda:

  The overall rate. One number holds it, several are the grid the path
  visits as they stand, and `NULL`, the default, has the path build one.
  Must lie in \\(0, \infty)\\.

- alpha:

  The mixing weight, in the same three states and settled independently
  of `lambda`. Must lie in \\(0, 1)\\.

- n_lambda, n_alpha:

  How many values the path visits for each, at least 2. `NULL`, the
  default, leaves it to the criterion.

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

\$\$\rho(\beta) = \lambda\left\\\alpha\lVert\beta\rVert_1 +
\frac{1-\alpha}{2}\lVert\beta\rVert_2^2\right\\ + p\log Z(\lambda,
\alpha),\$\$ the normalizing constant being that of the product of a
Laplace and a Gaussian at zero
([`enet_distrib`](https://statmodels7.github.io/distributions7/reference/enet_distrib.html)).
It depends on BOTH hyperparameters, which is what makes them estimable
rather than merely settable, and what a penalty written as a formula
would not have.

\\\alpha\\ is the mixing weight: at \\\alpha \to 1\\ the penalty is the
lasso and at \\\alpha \to 0\\ the ridge, and the kink at zero has
half-width \\\lambda\alpha\\, so both hyperparameters scale it.

**Hyperparameters.** `lambda` on \\(0, \infty)\\, swept over `n_lambda`
values by kink size; `alpha` on \\(0, 1)\\, swept over `n_alpha` values
across that interval, the ends excluded because the penalty there is one
of the other two. The two are swept cyclically, one at a time, which
keeps the cost linear in their number where a product grid would be
exponential in it.

## References

Zou, H. and Hastie, T. (2005). Regularization and variable selection via
the elastic net. *Journal of the Royal Statistical Society, Series B*
67, 301–320.

## See also

[`penalized_terms`](https://statmodels7.github.io/modelterms7/reference/penalized_terms.md)
for what the five share,
[`ridge`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`lasso`](https://statmodels7.github.io/modelterms7/reference/lasso.md),
[`scad`](https://statmodels7.github.io/modelterms7/reference/scad.md),
[`mcp`](https://statmodels7.github.io/modelterms7/reference/mcp.md),
[`elasticnet_penalty`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.html)

## Examples

``` r
dd <- data.frame(x1 = rnorm(8), x2 = rnorm(8))
term_penalty(term_build(enet(~ x1 + x2), dd))@params
#> [1] "lambda" "alpha" 

# alpha held at the halfway mixture, lambda still estimated
term_hyper(enet(~ x1 + x2, alpha = 0.5))
#> [[1]]
#> [[1]]$alpha
#> [1] 0.5
#> 
#> 
```
