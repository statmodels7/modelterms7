# MCP Penalty on a Block of Coefficients

The minimax concave penalty: like SCAD it selects and then flattens, and
it begins to flatten immediately, where SCAD waits for a first
threshold.

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
  formula throws.

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
  the path visits as they stand, and `NULL`, the default, has the path
  build one. Must lie in \\(0, \infty)\\.

- gamma:

  The shape, in the same three states and settled independently of
  `lambda`, defaulting to the literature's `3` rather than to `NULL`:
  one number holds it, several are the grid the path visits as they
  stand, and `NULL` has the path build one. Must lie in \\(1, \infty)\\.

- id:

  Labels sharing this term's hyperparameters with those of other terms
  carrying the same ones: each is then estimated at a single value,
  wherever in the formula they sit. A named vector, `c(gamma = "A")`,
  since the penalty carries several and which was meant is not a guess;
  `NULL`, the default, shares nothing. See
  [`term_ids()`](https://statmodels7.github.io/modelterms7/reference/term_ids.md),
  whose page says what the labels mean and what care they want.

- n_lambda, n_gamma:

  How many values the path visits for each, at least 2. They differ
  because the axes do: \\\lambda\\ descends the size of the kink over
  four decades and wants that many points, while \\\gamma\\ spans the
  shape's useful range and does not.

- min_ratio:

  How far down the path reaches, as a fraction of the kink that empties
  the block: smaller reaches a denser fit, larger stops sooner. A single
  number in \\(0, 1)\\, `1e-4` by default. Only the sweep by kink size
  reads it.

- search:

  `"grid"` to visit every combination of \\\lambda\\ and \\\gamma\\,
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

Defined by its derivative, for \\t = \lvert\beta_j\rvert \ge 0\\,
\$\$\rho'(t) = \left(\lambda - \frac{t}{\gamma}\right)\_+ ,\$\$ summed
over the coefficients: it starts at \\\lambda\\, falls linearly, and is
flat past \\t = \gamma\lambda\\. Improper by construction, so it carries
no normalizing constant and is not reachable by a marginal criterion.

**Hyperparameters.** `lambda` on \\(0, \infty)\\, swept over `n_lambda`
values by kink size; `gamma` on \\(1, \infty)\\, at or below which the
penalized objective need not be convex even for an orthogonal design,
swept over `n_gamma` values on a geometric grid above that bound.

**The shape is HELD at \\\gamma = 3\\ by default**, the value of Zhang,
so only \\\lambda\\ is searched. `gamma = NULL` estimates it over the
grid instead, and `n_gamma` sizes that grid. Measured over eight data
configurations, estimating it chose the grid's LOWER ENDPOINT every time
and changed no model: the columns kept were identical and the error
against the truth agreed to four decimals. The endpoint is the floor
plus 0.25, so what was reported as an estimate was the grid's own edge,
and it did not move with `n_gamma`, which refines the interior and
cannot touch either end.

⚠️ **The shape is not scale free**, which is why the literature's value
belongs to standardized data. Rescaling the response by \\k\\ sends
\\\beta \to k\beta\\, and the fit is reproduced at \\(\lambda/k,\\
\gamma k^2)\\, verified to 1.7e-13 at \\k = 2, 5, 10\\: the shape
carries the units of a proximal step, and the floor it is swept above
grows with them. On a response of much larger scale a small held shape
falls outside the convex region the compiled coordinate descent needs
and the fit takes the general route, measured at 7.4 s against 0.4 s at
a hundredfold response.

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
# gamma * lambda, where SCAD would still be on its first segment.
rho1 <- function(t, lambda, gamma) pmax(lambda - t / gamma, 0)
rho1(c(0, 1, 2, 3), lambda = 1, gamma = 2)
#> [1] 1.0 0.5 0.0 0.0

# The kink at zero, and the pair where the second derivative breaks.
penalties7::penalty_kinks(p, list(lambda = 1, gamma = 2))
#> [1]  0 -2  2

# The literature's shape, held.
term_hyper(mcp(~ x1 + x2, gamma = 3))
#> [[1]]
#> [[1]]$gamma
#> [1] 3
#> 
#> 

# At or below 1 the objective need not be convex, so it is refused.
try(mcp(~ x1 + x2, gamma = 1))
#> Error : 'gamma' in 'mcp' must lie strictly inside (1, Inf); it is 1.


# Fitted. The data are simulated from a known truth, so the
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
#>  1.97 -1.51  0.99  0.00  5.00 
```
