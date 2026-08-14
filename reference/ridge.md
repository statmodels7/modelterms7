# Ridge Penalty on a Block of Coefficients

A block of coefficients under a Gaussian prior at zero: every one of
them shrunk towards zero, none of them set to it.

## Usage

``` r
ridge(
  x,
  label = "ridge",
  by = NULL,
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

- by:

  Reserved for a later release; must be `NULL`.

- standardize:

  A single logical: whether to penalize each coefficient on the scale of
  its own column. See the section below.

- lambda:

  The precision of the prior, held at the value given and ESTIMATED when
  left `NULL`, which is the default. Must lie in \\(0, \infty)\\.

- n_lambda:

  Unused by this term: a ridge has no kink and its hyperparameter is
  estimated by a criterion rather than swept over a grid. Accepted so
  that the five constructors read alike.

- min_ratio:

  Unused by this term, which has no path. Accepted so that the five
  constructors read alike.

- ...:

  Not used, and reported: an argument named after another penalty's
  hyperparameter is the mistake this catches.

## Value

An object of class
[`PenalizedTerm`](https://statmodels7.github.io/modelterms7/reference/PenalizedTerm.md)
(a specification; see
[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

## Details

Writing \\\beta\\ for the block's coefficients and \\p\\ for their
number, \$\$\rho(\beta) = \frac{\lambda\lVert\beta\rVert_2^2}{2} -
\frac{p}{2}\log\\\left(\frac{\lambda}{2\pi}\right).\$\$ The constant is
kept, which is what makes \\\lambda\\ estimable by a marginal criterion:
it is minus the log density of \\N(0, \lambda^{-1}I)\\, so \\\lambda\\
is the PRECISION of that prior and a larger value shrinks harder.

The penalty is twice differentiable everywhere, so the block is fitted
in the same system as the unpenalized terms and \\\lambda\\ is estimated
by `reml()` rather than swept along a path.

**Hyperparameter.** `lambda`, admissible on \\(0, \infty)\\.

## References

Hoerl, A. E. and Kennard, R. W. (1970). Ridge regression: biased
estimation for nonorthogonal problems. *Technometrics* 12, 55–67.

## See also

[`penalized_terms`](https://statmodels7.github.io/modelterms7/reference/penalized_terms.md)
for what the five share,
[`lasso`](https://statmodels7.github.io/modelterms7/reference/lasso.md),
[`enet`](https://statmodels7.github.io/modelterms7/reference/enet.md),
[`scad`](https://statmodels7.github.io/modelterms7/reference/scad.md),
[`mcp`](https://statmodels7.github.io/modelterms7/reference/mcp.md),
[`ridge_penalty`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.html)

## Examples

``` r
dd <- data.frame(x1 = rnorm(8), x2 = rnorm(8))
term_penalty(term_build(ridge(~ x1 + x2), dd))@params
#> [1] "lambda"
term_hyper(ridge(~ x1 + x2, lambda = 2))
#> [[1]]
#> [[1]]$lambda
#> [1] 2
#> 
#> 
```
