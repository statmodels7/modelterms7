# MCP Penalty on a Block of Coefficients

The minimax concave penalty: like SCAD it selects and then flattens, and
it begins to flatten IMMEDIATELY rather than after a first threshold.

## Usage

``` r
mcp(
  x,
  label = "mcp",
  standardize = FALSE,
  lambda = NULL,
  gamma = NULL,
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

  A one-sided formula or a numeric matrix.

- label:

  A single non-empty string prefixed to the coefficient names.

- standardize:

  A single logical: whether to penalize each coefficient on the scale of
  its own column. See the section below.

- lambda:

  The scale of the penalty. One number holds it, several are the grid
  the path visits as they stand, and `NULL`, the default, has the path
  build one. Must lie in \\(0, \infty)\\.

- gamma:

  The shape, in the same three states and settled independently of
  `lambda`. Must lie in \\(1, \infty)\\.

- n_lambda, n_gamma:

  How many values the path visits for each, at least 2. They differ
  because the axes do: \\\lambda\\ descends the size of the kink over
  four decades and wants that many points, while \\\gamma\\ spans the
  shape's useful range and does not.

- min_ratio:

  How far down the path reaches, as a fraction of the kink that empties
  the block: smaller reaches a denser fit, larger stops sooner. Must lie
  in (0, 1). NULL, the default, leaves it to the criterion. Only the
  sweep by kink size uses it.

- search:

  `"grid"` to visit every combination of \\\lambda\\ and \\\gamma\\,
  `"cyclic"` to sweep one at a time with the other held. See
  [`term_search`](https://statmodels7.github.io/modelterms7/reference/term_search.md).

- sparse:

  Governs the FORMULA route: whether the block is built as a `dgCMatrix`
  through
  [`sparse.model.matrix`](https://rdrr.io/pkg/Matrix/man/sparse.model.matrix.html)
  rather than as a dense model matrix. `NULL`, the default, settles it
  at build from the size of the design; `TRUE` and `FALSE` override it.
  A MATRIX input needs no such argument, being kept in whatever storage
  it arrives in. See the section below.

- ...:

  Not used, and reported: an argument named after another penalty's
  hyperparameter is the mistake this catches.

## Value

An object of class
[`PenalizedTerm`](https://statmodels7.github.io/modelterms7/reference/PenalizedTerm.md)
(a specification; see
[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

## Details

Defined by its derivative, for \\t = \lvert\beta_j\rvert \ge 0\\,
\$\$\rho'(t) = \left(\lambda - \frac{t}{\gamma}\right)\_+ ,\$\$ summed
over the coefficients: it starts at \\\lambda\\, falls linearly, and is
flat past \\t = \gamma\lambda\\. Improper by construction, so it carries
no normalizing constant and is not reachable by a marginal criterion.

**Hyperparameters.** `lambda` on \\(0, \infty)\\, swept over `n_lambda`
values by kink size; `gamma` on \\(1, \infty)\\ – at \\\gamma \le 1\\
the penalized objective need not be convex even for an orthogonal design
– swept over `n_gamma` values on a geometric grid above that bound. The
literature's value is \\\gamma = 3\\.

## References

Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax
concave penalty. *The Annals of Statistics* 38, 894–942.

## See also

[`penalized_terms`](https://statmodels7.github.io/modelterms7/reference/penalized_terms.md)
for what the five share,
[`ridge`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`lasso`](https://statmodels7.github.io/modelterms7/reference/lasso.md),
[`enet`](https://statmodels7.github.io/modelterms7/reference/enet.md),
[`scad`](https://statmodels7.github.io/modelterms7/reference/scad.md),
[`mcp_penalty`](https://statmodels7.github.io/penalties7/reference/scad_penalty.html)

## Examples

``` r
dd <- data.frame(x1 = rnorm(8), x2 = rnorm(8))
term_penalty(term_build(mcp(~ x1 + x2), dd))@params
#> [1] "lambda" "gamma" 
term_hyper(mcp(~ x1 + x2, gamma = 3))
#> [[1]]
#> [[1]]$gamma
#> [1] 3
#> 
#> 
```
