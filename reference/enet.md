# Elastic Net Penalty on a Block of Coefficients

The lasso and the ridge mixed: a kink at zero, so the term still
selects, and a quadratic part that keeps correlated coefficients
together instead of choosing arbitrarily among them.

## Usage

``` r
enet(
  x,
  label = "enet",
  standardize = FALSE,
  lambda = NULL,
  alpha = NULL,
  id = NULL,
  n_lambda = 25,
  n_alpha = 5,
  min_ratio = 1e-04,
  search = "grid",
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
  as `label.name`, `"enet"` by default, so a block over `x1` reads
  `enet.x1`. Two penalized terms in one formula stay apart by their
  labels.

- standardize:

  A single logical, `FALSE` by default: whether to penalize each
  coefficient on the scale of its own column. See the section above.

- lambda:

  The overall rate. One number holds it, several are the grid the path
  visits as they stand, and `NULL`, the default, has the path build one.
  Must lie in \\(0, \infty)\\.

- alpha:

  The mixing weight, in the same three states and settled independently
  of `lambda`. Must lie in \\(0, 1)\\.

- id:

  Labels sharing this term's hyperparameters with those of other terms
  carrying the same ones: each is then estimated at a single value,
  wherever in the formula they sit. A named vector, `c(alpha = "A")`,
  since the penalty carries several and which was meant is not a guess;
  `NULL`, the default, shares nothing. See
  [`term_ids()`](https://statmodels7.github.io/modelterms7/reference/term_ids.md),
  whose page says what the labels mean and what care they want.

- n_lambda, n_alpha:

  How many values the path visits for each, at least 2. They differ
  because the axes do: \\\lambda\\ descends the size of the kink over
  four decades and wants that many points, while \\\alpha\\ spans one
  bounded interval and does not.

- min_ratio:

  How far down the path reaches, as a fraction of the kink that empties
  the block: smaller reaches a denser fit, larger stops sooner. A single
  number in \\(0, 1)\\, `1e-4` by default. Only the sweep by kink size
  reads it.

- search:

  `"grid"` to visit every combination of \\\lambda\\ and \\\alpha\\,
  `"cyclic"` to sweep one at a time with the other held. See
  [`term_search()`](https://statmodels7.github.io/modelterms7/reference/term_search.md).

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

\$\$\rho(\beta) = \lambda\left\\\alpha\lVert\beta\rVert_1 +
\frac{1-\alpha}{2}\lVert\beta\rVert_2^2\right\\ + p\log Z(\lambda,
\alpha),\$\$ the normalizing constant being that of the product of a
Laplace and a Gaussian at zero
([`distributions7::enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.html)).
It depends on both hyperparameters, so both are estimable, where a
merely settable one would be all a dropped constant leaves. A penalty
written as a formula, with the constant dropped, would not have that.

\\\alpha\\ is the mixing weight: at \\\alpha \to 1\\ the penalty is the
lasso and at \\\alpha \to 0\\ the ridge, and the kink at zero has
half-width \\\lambda\alpha\\, so both hyperparameters scale it.

**Hyperparameters.** `lambda` on \\(0, \infty)\\, swept over `n_lambda`
values by kink size; `alpha` on \\(0, 1)\\, swept over `n_alpha` values
across that interval, the ends excluded because the penalty there is one
of the other two. Every combination of the two is visited,
`n_lambda * n_alpha` fits, unless `search = "cyclic"` asks for one at a
time instead.

## References

Zou, H. and Hastie, T. (2005). Regularization and variable selection via
the elastic net. *Journal of the Royal Statistical Society, Series B*
67, 301–320.

## See also

[`penalized_terms()`](https://statmodels7.github.io/modelterms7/reference/penalized_terms.md)
for what the five share,
[`ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`lasso()`](https://statmodels7.github.io/modelterms7/reference/lasso.md),
[`scad()`](https://statmodels7.github.io/modelterms7/reference/scad.md),
[`mcp()`](https://statmodels7.github.io/modelterms7/reference/mcp.md),
[`penalties7::elasticnet_penalty()`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.html)

## Examples

``` r
set.seed(3)
dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
pe <- term_penalty(term_build(enet(~ x1 + x2), dd))
pe@params
#> [1] "lambda" "alpha" 
pe@params_bounds
#> $lambda
#> [1]   0 Inf
#> 
#> $alpha
#> [1] 0 1
#> 

# The two ends really are the other two penalties.
pl <- term_penalty(term_build(lasso(~ x1 + x2), dd))
pr <- term_penalty(term_build(ridge(~ x1 + x2), dd))
beta <- c(0.4, -1.1)
all.equal(penalties7::penalty_value(pe, beta, list(lambda = 2, alpha = 1 - 1e-9)),
          penalties7::penalty_value(pl, beta, list(lambda = 2)))
#> [1] TRUE
all.equal(penalties7::penalty_value(pe, beta, list(lambda = 2, alpha = 1e-9)),
          penalties7::penalty_value(pr, beta, list(lambda = 2)))
#> [1] TRUE

# alpha held at the halfway mixture, lambda still estimated.
term_hyper(enet(~ x1 + x2, alpha = 0.5))
#> [[1]]
#> [[1]]$alpha
#> [1] 0.5
#> 
#> 

# The grid is 25 by 5 by default; cyclic sweeps one axis at a time.
term_grid(enet(~ x1 + x2))
#> [[1]]
#> [[1]]$lambda
#> [1] 25
#> 
#> [[1]]$alpha
#> [1] 5
#> 
#> 
term_search(enet(~ x1 + x2, search = "cyclic"))
#> [[1]]
#> [1] "cyclic"
#> 

# The open interval is enforced: an end is one of the other penalties.
try(enet(~ x1 + x2, alpha = 1))
#> Error : 'alpha' in 'enet' must lie strictly inside (0, 1); it is 1.


# Fitted. The data are simulated from a known truth, so the
# estimates below can be read against it.
if (requireNamespace("statmodels7", quietly = TRUE)) {
  set.seed(11)
  XX <- matrix(rnorm(150 * 8), 150, 8)
  fd <- data.frame(y = as.numeric(XX %*% c(2, -1.5, 1, rep(0, 5))) +
                     rnorm(150, sd = 0.4))
  fd$X <- XX
  cf <- coef(statmodels7::statmod(y ~ enet(X, n_alpha = 3),
                                  distributions7::gaussian1_distrib(), fd))$mu
  # truth: the first three columns carry 2, -1.5 and 1, the other five
  # nothing. Both hyperparameters are chosen by BIC, over a
  # shortened alpha axis so that the product grid stays quick.
  c(round(cf[2:5], 2), kept = sum(cf[-1] != 0))
}
#> enet.1 enet.2 enet.3 enet.4   kept 
#>   1.94  -1.49   0.98   0.00   5.00 
```
