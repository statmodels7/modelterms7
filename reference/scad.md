# SCAD Penalty on a Block of Coefficients

Smoothly clipped absolute deviation: the lasso's kink at zero, so the
term selects, and a penalty that FLATTENS beyond a threshold, so a large
coefficient is not shrunk at all.

## Usage

``` r
scad(
  x,
  label = "scad",
  by = NULL,
  standardize = FALSE,
  lambda = NULL,
  a = NULL,
  n_lambda = NULL,
  n_a = NULL,
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

  The scale of the penalty. One number holds it, several are the grid
  the path visits as they stand, and `NULL`, the default, has the path
  build one. Must lie in \\(0, \infty)\\.

- a:

  The shape, in the same three states and settled independently of
  `lambda`. Must lie in \\(2, \infty)\\.

- n_lambda, n_a:

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

It is defined by its derivative rather than by its value, for \\t =
\lvert\beta_j\rvert \ge 0\\, \$\$\rho'(t) = \lambda\min\\\left\\1,
\frac{(a\lambda - t)\_+}{(a-1)\lambda}\right\\,\$\$ summed over the
coefficients. It rises like the lasso near zero, bends from \\t =
\lambda\\, and is flat past \\t = a\lambda\\: the bias the lasso puts on
a large coefficient is what this removes. Being improper it carries no
normalizing constant, and is therefore not a log prior and not reachable
by a marginal criterion.

**Hyperparameters.** `lambda` on \\(0, \infty)\\, swept over `n_lambda`
values by kink size; `a` on \\(2, \infty)\\ – below 2 the penalty is not
what its definition intends – swept over `n_a` values on a geometric
grid above that bound, since the shape leaves the kink at zero unchanged
and no kink-size path can reach it. The literature's value is \\a =
3.7\\, and holding it there is what ncvreg does.

## References

Fan, J. and Li, R. (2001). Variable selection via nonconcave penalized
likelihood and its oracle properties. *Journal of the American
Statistical Association* 96, 1348–1360.

## See also

[`penalized_terms`](https://statmodels7.github.io/modelterms7/reference/penalized_terms.md)
for what the five share,
[`ridge`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`lasso`](https://statmodels7.github.io/modelterms7/reference/lasso.md),
[`enet`](https://statmodels7.github.io/modelterms7/reference/enet.md),
[`mcp`](https://statmodels7.github.io/modelterms7/reference/mcp.md),
[`scad_penalty`](https://statmodels7.github.io/penalties7/reference/scad_penalty.html)

## Examples

``` r
dd <- data.frame(x1 = rnorm(8), x2 = rnorm(8))
term_penalty(term_build(scad(~ x1 + x2), dd))@params
#> [1] "lambda" "a"     
term_hyper(scad(~ x1 + x2, a = 3.7))
#> [[1]]
#> [[1]]$a
#> [1] 3.7
#> 
#> 
```
