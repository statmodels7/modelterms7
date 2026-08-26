# Ridge Penalty on a Block of Coefficients

A block of coefficients under a Gaussian prior at zero: every one of
them shrunk towards zero, none of them set to it.

## Usage

``` r
ridge(
  x,
  label = "ridge",
  standardize = FALSE,
  lambda = NULL,
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
  as `label.name`, `"ridge"` by default, so a block over `x1` reads
  `ridge.x1`. Two penalized terms in one formula stay apart by their
  labels.

- standardize:

  A single logical, `FALSE` by default: whether to penalize each
  coefficient on the scale of its own column. See the section above.

- lambda:

  The precision of the prior. One number holds it and `NULL`, the
  default, has it estimated. A ridge has no kink and no path, so several
  numbers are not a grid it could visit. Must lie in \\(0, \infty)\\.

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

Writing \\\beta\\ for the block's coefficients and \\p\\ for their
number, \$\$\rho(\beta) = \frac{\lambda\lVert\beta\rVert_2^2}{2} -
\frac{p}{2}\log\\\left(\frac{\lambda}{2\pi}\right).\$\$ The constant is
kept, and that is what a marginal criterion needs to estimate
\\\lambda\\: the value is minus the log density of \\N(0,
\lambda^{-1}I)\\, so \\\lambda\\ is the precision of that prior and a
larger value shrinks harder.

The penalty is twice differentiable everywhere, so the block is fitted
in the same system as the unpenalized terms and \\\lambda\\ is estimated
by `statmodels7::reml()` rather than swept along a path.

**Hyperparameter.** `lambda`, admissible on \\(0, \infty)\\.

## References

Hoerl, A. E. and Kennard, R. W. (1970). Ridge regression: biased
estimation for nonorthogonal problems. *Technometrics* 12, 55–67.

## See also

[`penalized_terms()`](https://statmodels7.github.io/modelterms7/reference/penalized_terms.md)
for what the five share,
[`lasso()`](https://statmodels7.github.io/modelterms7/reference/lasso.md),
[`enet()`](https://statmodels7.github.io/modelterms7/reference/enet.md),
[`scad()`](https://statmodels7.github.io/modelterms7/reference/scad.md),
[`mcp()`](https://statmodels7.github.io/modelterms7/reference/mcp.md),
[`penalties7::ridge_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.html)

## Examples

``` r
set.seed(3)
dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
b <- term_build(ridge(~ x1 + x2), dd)

# One hyperparameter, positive, on the log scale for an optimizer.
p <- term_penalty(b)
p@params
#> [1] "lambda"
p@params_bounds
#> $lambda
#> [1]   0 Inf
#> 

# The value is exactly minus a Gaussian log-density at precision lambda.
beta <- c(0.4, -1.1)
all.equal(penalties7::penalty_value(p, beta, list(lambda = 2.5)),
          -sum(dnorm(beta, 0, 1 / sqrt(2.5), log = TRUE)))
#> [1] TRUE

# No kink, so the block is smooth and lambda is estimated at the mode.
c(smooth = term_smooth(b),
  kinks = length(penalties7::penalty_kinks(p, list(lambda = 1))))
#> smooth  kinks 
#>      1      0 

# Holding it, and the refusal of a grid it has no path to visit.
term_hyper(ridge(~ x1 + x2, lambda = 2))
#> [[1]]
#> [[1]]$lambda
#> [1] 2
#> 
#> 
try(ridge(~ x1 + x2, lambda = c(1, 2)))
#> Error : 'lambda' in 'ridge' has several values, and this term has no path to visit them
#>   on: its penalty has no kink, so the hyperparameter is estimated by the
#>   criterion at the mode. Give one number to hold it, or NULL to estimate it.
```
