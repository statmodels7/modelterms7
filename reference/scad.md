# SCAD Penalty on a Block of Coefficients

Smoothly clipped absolute deviation: the lasso's kink at zero, so the
term selects, and a penalty that flattens beyond a threshold, so a large
coefficient is not shrunk at all.

## Usage

``` r
scad(
  x,
  label = "scad",
  standardize = FALSE,
  lambda = NULL,
  a = NULL,
  n_lambda = 25,
  n_a = 5,
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
  as `label.name`, `"scad"` by default, so a block over `x1` reads
  `scad.x1`. Two penalized terms in one formula stay apart by their
  labels.

- standardize:

  A single logical, `FALSE` by default: whether to penalize each
  coefficient on the scale of its own column. See the section above.

- lambda:

  The scale of the penalty. One number holds it, several are the grid
  the path visits as they stand, and `NULL`, the default, has the path
  build one. Must lie in \\(0, \infty)\\.

- a:

  The shape, in the same three states and settled independently of
  `lambda`. Must lie in \\(2, \infty)\\.

- n_lambda, n_a:

  How many values the path visits for each, at least 2. They differ
  because the axes do: \\\lambda\\ descends the size of the kink over
  four decades and wants that many points, while \\a\\ spans the shape's
  useful range and does not.

- min_ratio:

  How far down the path reaches, as a fraction of the kink that empties
  the block: smaller reaches a denser fit, larger stops sooner. A single
  number in \\(0, 1)\\, `1e-4` by default. Only the sweep by kink size
  reads it.

- search:

  `"grid"` to visit every combination of \\\lambda\\ and \\a\\,
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

It is defined by its derivative, for \\t = \lvert\beta_j\rvert \ge 0\\,
\$\$\rho'(t) = \lambda\min\\\left\\1, \frac{(a\lambda -
t)\_+}{(a-1)\lambda}\right\\,\$\$ summed over the coefficients. It rises
like the lasso near zero, bends from \\t = \lambda\\, and is flat past
\\t = a\lambda\\, which removes the bias the lasso puts on a large
coefficient. Being improper it carries no normalizing constant, and is
therefore not a log prior and not reachable by a marginal criterion.

**Hyperparameters.** `lambda` on \\(0, \infty)\\, swept over `n_lambda`
values by kink size; `a` on \\(2, \infty)\\, below which the penalty is
not what its definition intends, swept over `n_a` values on a geometric
grid above that bound, since the shape leaves the kink at zero unchanged
and no kink-size path can reach it. The literature's value is \\a =
3.7\\, and holding it there is what ncvreg does.

## References

Fan, J. and Li, R. (2001). Variable selection via nonconcave penalized
likelihood and its oracle properties. *Journal of the American
Statistical Association* 96, 1348–1360.

## See also

[`penalized_terms()`](https://statmodels7.github.io/modelterms7/reference/penalized_terms.md)
for what the five share,
[`ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`lasso()`](https://statmodels7.github.io/modelterms7/reference/lasso.md),
[`enet()`](https://statmodels7.github.io/modelterms7/reference/enet.md),
[`mcp()`](https://statmodels7.github.io/modelterms7/reference/mcp.md),
[`penalties7::scad_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.html)

## Examples

``` r
set.seed(3)
dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
p <- term_penalty(term_build(scad(~ x1 + x2), dd))
p@params
#> [1] "lambda" "a"     
p@params_bounds
#> $lambda
#> [1]   0 Inf
#> 
#> $a
#> [1]   2 Inf
#> 

# The derivative is the definition: lasso-like to lambda, then bending,
# then flat past a * lambda.
rho1 <- function(t, lambda, a)
  lambda * pmin(1, pmax(a * lambda - t, 0) / ((a - 1) * lambda))
rho1(c(0.5, 1.5, 4), lambda = 1, a = 3.7)
#> [1] 1.0000000 0.8148148 0.0000000

# penalty_kinks() reports every point where some derivative breaks:
# the origin, and the two pairs where the second derivative changes
# branch, at plus and minus lambda and a * lambda.
penalties7::penalty_kinks(p, list(lambda = 1, a = 3))
#> [1]  0 -1  1 -3  3

# The literature's shape, held.
term_hyper(scad(~ x1 + x2, a = 3.7))
#> [[1]]
#> [[1]]$a
#> [1] 3.7
#> 
#> 

# Below 2 the shape is refused.
try(scad(~ x1 + x2, a = 2))
#> Error : 'a' in 'scad' must lie strictly inside (2, Inf); it is 2.
```
