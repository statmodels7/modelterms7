# MCP Penalty on a Block of Coefficients

Specifies the minimax concave penalty. Like SCAD it selects and then
flattens, and it begins to flatten immediately, whereas SCAD keeps the
slope of the lasso up to \\t = \lambda\\.

## Usage

``` r
mcp(
  x,
  label = "mcp",
  standardize = FALSE,
  lambda = NULL,
  gamma = 3,
  id = NULL,
  n_lambda = 25,
  n_gamma = 5,
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
  as `label.name`, `"mcp"` by default, so a block over `x1` reads
  `mcp.x1`. Two penalized terms in one formula stay apart by their
  labels.

- standardize:

  A single logical, `FALSE` by default: whether to penalize each
  coefficient on the scale of its own column. See the section above.

- lambda:

  The scale of the penalty. One number holds it, several are the grid
  that the path visits as they stand, and `NULL`, the default, has the
  path build one. Must lie in \\(0, \infty)\\.

- gamma:

  The shape, in the same three states and settled independently of
  `lambda`, defaulting to the literature's `3` rather than to `NULL`:
  one number holds it, several are the grid that the path visits as they
  stand, and `NULL` has the path build one. Must lie in \\(1, \infty)\\.

- id:

  Labels sharing this term's hyperparameters with those of other terms
  carrying the same ones: each is then estimated at a single value,
  wherever in the formula they sit. The argument is a named vector,
  `c(gamma = "A")`, because the penalty has several hyperparameters and
  the name identifies the one that is shared; `NULL`, the default,
  shares nothing. See
  [`term_ids()`](https://statmodels7.github.io/modelterms7/reference/term_ids.md)
  for the meaning of the labels and the restrictions on their use.

- n_lambda, n_gamma:

  The number of values that the path visits for each, at least 2. The
  two are separate arguments because the axes differ: \\\lambda\\
  descends from the kink size that empties the block to `min_ratio`
  times that size (four decades at the default `min_ratio`), which takes
  more points than \\\gamma\\, which spans the shape's useful range.
  Each signals an error when given explicitly for a hyperparameter that
  is held or written out, so `n_gamma` requires `gamma = NULL`.

- min_ratio:

  How far down the path reaches, as a fraction of the kink that empties
  the block: smaller reaches a denser fit, larger stops sooner. It is a
  single number in \\(0, 1)\\, `1e-4` by default. Only the sweep by kink
  size reads it, and it signals an error when given explicitly while
  `lambda` is held or written out.

- search:

  `"grid"` to visit every combination of \\\lambda\\ and \\\gamma\\,
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

The penalty is defined by its derivative, for \\t = \lvert\beta_j\rvert
\ge 0\\, \$\$\rho'(t) = \left(\lambda - \frac{t}{\gamma}\right)\_+ ,\$\$
summed over the coefficients: it starts at \\\lambda\\, falls linearly,
and is flat past \\t = \gamma\lambda\\. The penalty is improper, so it
carries no normalizing constant and is not reachable by a marginal
criterion.

The formula above is the penalty for a loss of unit curvature. In a
[`statmodels7::statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.html)
fit the penalty on coordinate \\j\\ is scaled by the curvature \\c_j\\
of the likelihood in that coordinate: the slope at zero remains
\\\lambda\\, and the penalty is flat past \\\gamma\lambda/c_j\\, so the
shape keeps its meaning whatever the scale of the response. Rescaling a
Gaussian response by \\k\\ and dividing \\\lambda\\ by \\k\\ reproduces
the fit, in which \\\gamma\\ is unchanged. The scaling is applied by
statmodels7; the penalty that this term builds is the unscaled one.

The hyperparameters are `lambda`, on \\(0, \infty)\\, swept over
`n_lambda` values by kink size, and `gamma`, on \\(1, \infty)\\. At or
below 1 the penalized objective need not be convex even for an
orthogonal design. The shape leaves the kink at zero unchanged, so a
path over the kink size cannot reach it. It is swept over `n_gamma`
values on a geometric grid that starts just above the smallest shape at
which the block can be fitted.

The shape is held at \\\gamma = 3\\ by default, the value that ncvreg
holds, so only \\\lambda\\ is searched. `gamma = NULL` estimates it over
the grid instead, and `n_gamma` sizes that grid. The grid runs from 0.25
to 25 above the smallest admissible shape.

## References

Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax
concave penalty. *The Annals of Statistics* 38, 894–942.

## See also

[`penalized_terms()`](https://statmodels7.github.io/modelterms7/reference/penalized_terms.md)
for what the five share,
[`ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`lasso()`](https://statmodels7.github.io/modelterms7/reference/lasso.md),
[`enet()`](https://statmodels7.github.io/modelterms7/reference/enet.md),
[`scad()`](https://statmodels7.github.io/modelterms7/reference/scad.md),
[`penalties7::mcp_penalty()`](https://statmodels7.github.io/penalties7/reference/scad_penalty.html)

## Examples

``` r
set.seed(3)
dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
p <- term_penalty(term_build(mcp(~ x1 + x2), dd))
p@params
#> [1] "lambda" "gamma" 
p@params_bounds
#> $lambda
#> [1]   0 Inf
#> 
#> $gamma
#> [1]   1 Inf
#> 

# The derivative starts at lambda and falls linearly to zero at
# gamma * lambda. SCAD's derivative stays at lambda up to t = lambda.
rho1 <- function(t, lambda, gamma) pmax(lambda - t / gamma, 0)
rho1(c(0, 1, 2, 3), lambda = 1, gamma = 2)
#> [1] 1.0 0.5 0.0 0.0

# penalty_kinks() reports the kink at zero and the pair where the second
# derivative breaks.
penalties7::penalty_kinks(p, list(lambda = 1, gamma = 2))
#> [1]  0 -2  2

# The literature's shape is held by default, and `gamma = NULL` estimates it
# over a grid.
term_hyper(mcp(~ x1 + x2))
#> [[1]]
#> [[1]]$gamma
#> [1] 3
#> 
#> 
term_hyper(mcp(~ x1 + x2, gamma = NULL))
#> list()
term_grid(mcp(~ x1 + x2, gamma = NULL))
#> [[1]]
#> [[1]]$lambda
#> [1] 25
#> 
#> [[1]]$gamma
#> [1] 5
#> 
#> 

# At or below 1 the objective need not be convex, so it is rejected.
try(mcp(~ x1 + x2, gamma = 1))
#> Error : 'gamma' in 'mcp' must lie strictly inside (1, Inf); it is 1.


# The model is fitted. The data are simulated from a known truth, so the
# estimates below can be read against it.
if (requireNamespace("statmodels7", quietly = TRUE)) {
  set.seed(11)
  XX <- matrix(rnorm(150 * 8), 150, 8)
  fd <- data.frame(y = as.numeric(XX %*% c(2, -1.5, 1, rep(0, 5))) +
                     rnorm(150, sd = 0.4))
  fd$X <- XX
  cf <- coef(statmodels7::statmod(y ~ mcp(X),
                                  distributions7::gaussian1_distrib(), fd))$mu
  # truth: the first three columns carry 2, -1.5 and 1, the other five
  # nothing. Only lambda is searched: the shape is held at the
  # literature's 3, and `gamma = NULL` would estimate it instead.
  c(round(cf[2:5], 2), kept = sum(cf[-1] != 0))
}
#> mcp.1 mcp.2 mcp.3 mcp.4  kept 
#>  2.00 -1.54  1.02  0.00  5.00 
```
