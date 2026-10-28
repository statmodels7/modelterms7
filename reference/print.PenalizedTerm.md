# Print a Penalized Term

Prints one line describing a
[`ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`lasso()`](https://statmodels7.github.io/modelterms7/reference/lasso.md),
[`enet()`](https://statmodels7.github.io/modelterms7/reference/enet.md),
[`scad()`](https://statmodels7.github.io/modelterms7/reference/scad.md)
or [`mcp()`](https://statmodels7.github.io/modelterms7/reference/mcp.md)
term. An unbuilt specification reports its label and whether it will
standardize; a built one reports the number of coefficients, the
penalty's name and the hyperparameters that the penalty carries, and,
when the term is standardized, adds a second line giving the spread by
which each column is standardized.

## Arguments

- x:

  A
  [`PenalizedTerm()`](https://statmodels7.github.io/modelterms7/reference/PenalizedTerm.md),
  built or not.

- ...:

  Unused, and accepted so that the signature matches
  [`print()`](https://rdrr.io/r/base/print.html)'s.

## Value

`x`, invisibly. The method is called for the line that it writes.

## Details

The built form reads `ncol(x@X)`, `x@penalty@penalty_name` and
`x@penalty@params`, so the hyperparameters shown are the penalty's own
names: `lambda` for a ridge or a lasso, `lambda, alpha` for the elastic
net, `lambda, a` for SCAD and `lambda, gamma` for MCP.

The standardization line is needed to read a fitted hyperparameter. With
`standardize = TRUE` the penalty is evaluated at the coefficients
multiplied by those spreads; these products are the coefficients of the
columns after each is divided by its spread. The spreads are computed
once at build time and recorded in the `build_info`. They enter the
penalty only: the design block is not rescaled, so prediction at new
rows does not use them.

## See also

[`penalized_terms()`](https://statmodels7.github.io/modelterms7/reference/penalized_terms.md)
for the five constructors and what they share.

## Examples

``` r
set.seed(5)
d <- data.frame(x1 = rnorm(30), x2 = rnorm(30) * 20)

# A specification prints what it is, without any quantity from data.
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
