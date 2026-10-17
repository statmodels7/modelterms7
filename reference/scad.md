# SCAD Penalty on a Block of Coefficients

Specifies the smoothly clipped absolute deviation penalty. It has the
lasso's kink at zero, so the term selects, and it flattens beyond a
threshold, so a large coefficient is not shrunk at all.

## Usage

``` r
scad(
  x,
  label = "scad",
  standardize = FALSE,
  lambda = NULL,
  a = 3.7,
  id = NULL,
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
  formula is rejected with an error.

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
  that the path visits as they stand, and `NULL`, the default, has the
  path build one. Must lie in \\(0, \infty)\\.

- a:

  The shape, in the same three states and settled independently of
  `lambda`, defaulting to the literature's `3.7` rather than to `NULL`:
  one number holds it, several are the grid that the path visits as they
  stand, and `NULL` has the path build one. Must lie in \\(2, \infty)\\.

- id:

  Labels sharing this term's hyperparameters with those of other terms
  carrying the same ones: each is then estimated at a single value,
  wherever in the formula they sit. The argument is a named vector,
  `c(a = "A")`, because the penalty has several hyperparameters and the
  name identifies the one that is shared; `NULL`, the default, shares
  nothing. See
  [`term_ids()`](https://statmodels7.github.io/modelterms7/reference/term_ids.md)
  for the meaning of the labels and the restrictions on their use.

- n_lambda, n_a:

  The number of values that the path visits for each, at least 2. The
  two are separate arguments because the axes differ: \\\lambda\\
  descends from the kink size that empties the block to `min_ratio`
  times that size (four decades at the default `min_ratio`), which takes
  more points than \\a\\, which spans the shape's useful range. Each
  signals an error when given explicitly for a hyperparameter that is
  held or written out, so `n_a` requires `a = NULL`.

- min_ratio:

  How far down the path reaches, as a fraction of the kink that empties
  the block: smaller reaches a denser fit, larger stops sooner. It is a
  single number in \\(0, 1)\\, `1e-4` by default. Only the sweep by kink
  size reads it, and it signals an error when given explicitly while
  `lambda` is held or written out.

- search:

  `"grid"` to visit every combination of \\\lambda\\ and \\a\\,
  `"cyclic"` to sweep one at a time while the other is held. See
  [`term_search()`](https://statmodels7.github.io/modelterms7/reference/term_search.md).
  Given explicitly while fewer than two hyperparameters are swept, it
  signals an error.

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
[`PenalizedTerm()`](https://statmodels7.github.io/modelterms7/reference/PenalizedTerm.md):
a specification whose `X`, `coef_names`, `blueprint` and `penalty` are
empty until
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
fills them. The build attaches the penalty over the block's columns.

## Details

It is defined by its derivative, for \\t = \lvert\beta_j\rvert \ge 0\\,
\$\$\rho'(t) = \lambda\min\\\left\\1, \frac{(a\lambda -
t)\_+}{(a-1)\lambda}\right\\,\$\$ summed over the coefficients. It rises
like the lasso near zero, bends from \\t = \lambda\\, and is flat past
\\t = a\lambda\\, which removes the bias that the lasso puts on a large
coefficient. The penalty is improper, so it carries no normalizing
constant, and it is therefore not a log prior and not reachable by a
marginal criterion.

The formula above is the penalty for a loss of unit curvature. In a
[`statmodels7::statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.html)
fit the penalty on coordinate \\j\\ is scaled by the curvature \\c_j\\
of the likelihood in that coordinate: the slope at zero remains
\\\lambda\\, the penalty bends from \\\lambda/c_j\\ and is flat past
\\a\lambda/c_j\\, so the shape keeps its meaning whatever the scale of
the response. Rescaling a Gaussian response by \\k\\ and dividing
\\\lambda\\ by \\k\\ reproduces the fit, in which \\a\\ is unchanged.
The scaling is applied by statmodels7; the penalty that this term builds
is the unscaled one.

The hyperparameters are `lambda`, on \\(0, \infty)\\, swept over
`n_lambda` values by kink size, and `a`, on \\(2, \infty)\\, the range
over which the proximal operator at unit curvature is single valued. The
shape leaves the kink at zero unchanged, so a path over the kink size
cannot reach it. It is swept over `n_a` values on a geometric grid that
starts just above the smallest shape at which the block can be fitted.

The shape is held at \\a = 3.7\\ by default, the value of Fan and Li and
the value that ncvreg holds, so only \\\lambda\\ is searched. `a = NULL`
estimates it over the grid instead, and `n_a` sizes that grid. The grid
runs from 0.25 to 25 above the smallest admissible shape.

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

# The derivative is the definition: it behaves like the lasso up to
# lambda, then bends, and is flat past a * lambda.
rho1 <- function(t, lambda, a)
  lambda * pmin(1, pmax(a * lambda - t, 0) / ((a - 1) * lambda))
rho1(c(0.5, 1.5, 4), lambda = 1, a = 3.7)
#> [1] 1.0000000 0.8148148 0.0000000

# penalty_kinks() reports every point where some derivative breaks:
# the origin, and the two pairs where the second derivative changes
# branch, at plus and minus lambda and a * lambda.
penalties7::penalty_kinks(p, list(lambda = 1, a = 3))
#> [1]  0 -1  1 -3  3

# The literature's shape is held by default, and `a = NULL` estimates it
# over a grid.
term_hyper(scad(~ x1 + x2))
#> [[1]]
#> [[1]]$a
#> [1] 3.7
#> 
#> 
term_hyper(scad(~ x1 + x2, a = NULL))
#> list()
term_grid(scad(~ x1 + x2, a = NULL))
#> [[1]]
#> [[1]]$lambda
#> [1] 25
#> 
#> [[1]]$a
#> [1] 5
#> 
#> 

# At or below 2 the shape is rejected, because the bound is open.
try(scad(~ x1 + x2, a = 2))
#> Error : 'a' in 'scad' must lie strictly inside (2, Inf); it is 2.


# The model is fitted. The data are simulated from a known truth, so the
# estimates below can be read against it.
if (requireNamespace("statmodels7", quietly = TRUE)) {
  set.seed(11)
  XX <- matrix(rnorm(150 * 8), 150, 8)
  fd <- data.frame(y = as.numeric(XX %*% c(2, -1.5, 1, rep(0, 5))) +
                     rnorm(150, sd = 0.4))
  fd$X <- XX
  cf <- coef(statmodels7::statmod(y ~ scad(X),
                                  distributions7::gaussian1_distrib(), fd))$mu
  # truth: the first three columns carry 2, -1.5 and 1, the other five
  # nothing. Only lambda is searched: the shape is held at the
  # literature's 3.7, and `a = NULL` would estimate it instead.
  c(round(cf[2:5], 2), kept = sum(cf[-1] != 0))
}
#> scad.1 scad.2 scad.3 scad.4   kept 
#>   2.00  -1.54   1.02   0.00   5.00 
```
