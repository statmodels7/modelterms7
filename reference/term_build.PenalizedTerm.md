# Build a Penalized Block and Attach Its Penalty

Builds the block of a
[`ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`lasso()`](https://statmodels7.github.io/modelterms7/reference/lasso.md),
[`enet()`](https://statmodels7.github.io/modelterms7/reference/enet.md),
[`scad()`](https://statmodels7.github.io/modelterms7/reference/scad.md)
or [`mcp()`](https://statmodels7.github.io/modelterms7/reference/mcp.md)
term, names its coefficients, records the blueprint, and calls the
term's factory to produce the penalties7 penalty over exactly that many
coefficients. The penalty exists only after this: a specification's
`penalty` property is `NULL`, its width being unknown until the data are
seen.

## Arguments

- term:

  An unbuilt or built
  [`PenalizedTerm()`](https://statmodels7.github.io/modelterms7/reference/PenalizedTerm.md).

- data:

  A data frame with as many rows as a matrix input has, carrying every
  variable a formula input names.

- ...:

  Unused.

## Value

The term with `X`, `coef_names`, `blueprint` and `penalty` filled.
`penalty@n_coef` equals `ncol(X)`.

## Two input routes, two blueprints

A **formula** is built through
[`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html) and
[`.design_matrix()`](https://statmodels7.github.io/modelterms7/reference/dot-design_matrix.md)
with the intercept removed by `update(formula, ~ . - 1)`, unless the
intercept is all the formula has, in which case it is kept and is the
block. The blueprint records `kind = "formula"`, the terms object, the
levels, the contrasts and the settled storage.

A **matrix** is used as it stands, and its row count must equal
`nrow(data)`; anything else throws with both numbers. Its columns take
the matrix's own names, or `1`, `2`, ... where it has none. The
blueprint records `kind = "matrix"`, the expression that produced it and
those base names.

Either way every column name is prefixed with the term's label, and the
block's row names are dropped.

## Standardization is applied to the penalty

With `standardize = TRUE` the column spreads come from
[`.block_sd()`](https://statmodels7.github.io/modelterms7/reference/dot-block_sd.md),
go into `blueprint$standardize` and become a
[`Matrix::Diagonal()`](https://rdrr.io/pkg/Matrix/man/Diagonal.html)
map. The design is not touched, so a sparse block stays sparse.

The map is passed to the factory **as a second argument, and only when
there is one**. A factory written before standardization existed takes
the count alone and goes on working, provided it is never standardized.

## See also

[`penalized_terms()`](https://statmodels7.github.io/modelterms7/reference/penalized_terms.md)
for what the five constructors share,
[`term_predict.PenalizedTerm()`](https://statmodels7.github.io/modelterms7/reference/term_predict.PenalizedTerm.md)
for the block at new rows,
[`.block_sd()`](https://statmodels7.github.io/modelterms7/reference/dot-block_sd.md)
for the spreads.

## Examples

``` r
set.seed(3)
dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20))

# The penalty appears at build, over exactly the columns built.
b <- term_build(lasso(~ x1 + x2), dd)
c(cols = ncol(term_matrix(b)), n_coef = b@penalty@n_coef)
#>   cols n_coef 
#>      2      2 

# A matrix input must have as many rows as the data.
M <- matrix(rnorm(30), 10, 3)
try(term_build(ridge(M), dd))
#> Error : the matrix input has 10 rows and 'data' has 20.
```
