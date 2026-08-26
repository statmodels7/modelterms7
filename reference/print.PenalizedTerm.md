# Print a Penalized Term

Prints one line describing a
[`ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`lasso()`](https://statmodels7.github.io/modelterms7/reference/lasso.md),
[`enet()`](https://statmodels7.github.io/modelterms7/reference/enet.md),
[`scad()`](https://statmodels7.github.io/modelterms7/reference/scad.md)
or [`mcp()`](https://statmodels7.github.io/modelterms7/reference/mcp.md)
term. An unbuilt specification reports its label and whether it will
standardize; a built one reports the number of coefficients, the
penalty's name and the hyperparameters that penalty carries, and adds a
second line giving the spread each column was standardized by when there
is one.

## Arguments

- x:

  A
  [`PenalizedTerm()`](https://statmodels7.github.io/modelterms7/reference/PenalizedTerm.md),
  built or not.

- ...:

  Unused, and accepted so that the signature matches
  [`print()`](https://rdrr.io/r/base/print.html)'s.

## Value

`x`, invisibly. Called for the line it writes.

## Details

The built form reads `ncol(x@X)`, `x@penalty@penalty_name` and
`x@penalty@params`, so the hyperparameters shown are the penalty's own
names, `lambda` for a ridge or a lasso, `lambda, alpha` for the elastic
net, `lambda, a` for SCAD and `lambda, gamma` for MCP.

The standardization line matters for reading a fitted hyperparameter:
with `standardize = TRUE` the penalty acts on the coefficients of
columns divided by those spreads, and the spreads are frozen in the
blueprint at build time, so prediction at new rows uses the same
numbers.

## See also

[`penalized_terms()`](https://statmodels7.github.io/modelterms7/reference/penalized_terms.md)
for the five constructors and what they share.

## Examples

``` r
set.seed(5)
d <- data.frame(x1 = rnorm(30), x2 = rnorm(30) * 20)

# A specification says only what it is.
ridge(~ x1 + x2)
#> <PenalizedTerm> 'ridge' (specification; call term_build() with data)
lasso(~ x1 + x2, standardize = TRUE)
#> <PenalizedTerm> 'lasso', standardized (specification; call term_build() with data)

# A built term names its penalty and that penalty's hyperparameters.
term_build(ridge(~ x1 + x2), d)
#> <PenalizedTerm> 'ridge' built: 2 coefficients; penalty quadratic (lambda)
term_build(scad(~ x1), d)
#> <PenalizedTerm> 'scad' built: 1 coefficient; penalty SCAD (lambda, a)

# Standardizing adds the spreads: x2 was simulated twenty times wider.
term_build(enet(~ x1 + x2, standardize = TRUE), d)
#> <PenalizedTerm> 'enet' built: 2 coefficients; penalty separable [fixed enet [mu=0]] (lambda, alpha)
#>   standardized by: enet.x1 = 0.9915, enet.x2 = 21
```
