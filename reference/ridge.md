# Penalized Parametric Terms

The four classical penalized blocks as model terms: ridge, lasso, SCAD
and MCP. Each takes its block as a one-sided formula or as a numeric
matrix, and attaches the corresponding penalties7 object to the block's
coefficients at build time –
[`ridge_penalty`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.html),
[`lasso_penalty`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.html),
[`elasticnet_penalty`](https://statmodels7.github.io/penalties7/reference/ridge_penalty.html),
[`scad_penalty`](https://statmodels7.github.io/penalties7/reference/scad_penalty.html),
[`mcp_penalty`](https://statmodels7.github.io/penalties7/reference/scad_penalty.html)
– so the hyperparameters, their bounds and links, the derivatives and
the kink set are the penalty's, never restated by the term.

## Usage

``` r
ridge(x, label = "ridge", by = NULL)

lasso(x, label = "lasso", by = NULL)

enet(x, label = "enet", by = NULL)

scad(x, label = "scad", by = NULL)

mcp(x, label = "mcp", by = NULL)
```

## Arguments

- x:

  A one-sided formula or a numeric matrix.

- label:

  A single non-empty string prefixed to the coefficient names.

- by:

  Reserved for a later release; must be `NULL`.

## Value

An object of class
[`PenalizedTerm`](https://statmodels7.github.io/modelterms7/reference/PenalizedTerm.md)
(a specification; see
[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

## Details

A formula input goes through the
[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html) machinery
with the intercept removed (a penalized block does not penalize an
intercept; the model's intercept lives in the parametric block), and its
blueprint records the terms, the factor levels and the contrasts,
exactly as
[`linpar`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
does. A matrix input is used as given, and its columns are named after
the matrix's own column names, or numbered when it has none.

Prediction for a matrix input re-evaluates the expression that produced
the matrix in the new data, and ONLY there, so the intended use is a
matrix column of the model data frame (`dd$R <- R`; then `ridge(R)` in
the formula): a subset of the data then carries the matching rows. A
free-standing matrix from the calling environment builds, since its
value was captured, but prediction is rejected – resolving it outside
the new data would silently reuse the build-time rows.

[`term_smooth`](https://statmodels7.github.io/modelterms7/reference/term_smooth.md)
is `TRUE` for `ridge` and `FALSE` for `lasso`, `scad` and `mcp`, read
from each penalty's kink set.

## Examples

``` r
dd <- data.frame(x1 = rnorm(8), x2 = rnorm(8))
built <- term_build(lasso(~ x1 + x2), dd)
term_coef_names(built)
#> [1] "lasso.x1" "lasso.x2"
term_penalty(built)@params
#> [1] "lambda"
term_smooth(built)
#> [1] FALSE
```
