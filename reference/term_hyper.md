# The Hyperparameters a Term Holds

Reports the hyperparameter values the caller fixed in the constructor,
one entry per penalty the term carries. A hyperparameter left `NULL` is
**estimated** by whatever criterion the fit runs and does not appear
here, so an empty result means every one of them is free.

## Usage

``` r
term_hyper(term, ...)
```

## Arguments

- term:

  A term, built or not.

- ...:

  Passed to methods. No shipped method reads anything here.

## Value

A named list, one entry per penalty of the term and named as
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
names it, each entry a named list of single numbers keyed by the
penalty's own hyperparameter names. Empty where the term holds nothing.

## The three states of one argument

A constructor's hyperparameter argument carries three meanings, read per
hyperparameter at a time. `NULL`, the default, has the path build a grid
or a marginal criterion estimate the value at the mode. **One number**
holds it, and that is what this reports. **Several numbers** are a grid
the path visits as they stand, which
[`term_values()`](https://statmodels7.github.io/modelterms7/reference/term_values.md)
reports; the hyperparameter is still estimated there, the caller having
said where to look, leaving the answer to the fit.

## How the entries are keyed

The names are the ones
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
gives its entries: `""` for a penalty covering the whole term, and
`parameter::label` for one a subformula brought in. A term whose
penalties come from sub-terms, which is what a structural term with
subformulas is, propagates their held values without a method of its
own, every entry of that enumeration carrying its own in the field
`fixed`.

A built term answers from those entries and an unbuilt one from its own
`hyper` property, so the answer is the same before and after a build.

## See also

[`term_values()`](https://statmodels7.github.io/modelterms7/reference/term_values.md)
for a written-out grid,
[`term_grid()`](https://statmodels7.github.io/modelterms7/reference/term_grid.md)
for how many values a built grid visits,
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
for the entries this is keyed by.

## Examples

``` r
# Held, and free.
term_hyper(lasso(~ x, lambda = 3))
#> [[1]]
#> [[1]]$lambda
#> [1] 3
#> 
#> 
term_hyper(lasso(~ x))
#> list()

# Several values are a grid, not a held value, so they go elsewhere.
term_hyper(lasso(~ x, lambda = c(0.1, 1, 10)))
#> list()
term_values(lasso(~ x, lambda = c(0.1, 1, 10)))
#> [[1]]
#> [[1]]$lambda
#> [1]  0.1  1.0 10.0
#> 
#> 

# One of two held, the other estimated.
term_hyper(enet(~ x, alpha = 0.5))
#> [[1]]
#> [[1]]$alpha
#> [1] 0.5
#> 
#> 

# A held value on a sub-term is reported under that entry's name.
set.seed(4)
d <- data.frame(x = runif(40, 0, 5), g = factor(rep(1:4, 10)))
nb <- term_build(nl(~ a * exp(-r * x), a ~ 0 + lasso(~ g, lambda = 2),
                    r ~ 1), d)
term_hyper(nb)
#> $`a::lasso(~g, lambda = 2)`
#> $`a::lasso(~g, lambda = 2)`$lambda
#> [1] 2
#> 
#> 
```
