# The Hyperparameter Values a Term Writes Out

Reports the values a path visits, for each hyperparameter the caller
wrote out as a vector instead of leaving to be built. It is the third
state of a constructor's hyperparameter argument, beside `NULL` for a
built grid and one number for a held value.

## Usage

``` r
term_values(term, ...)
```

## Arguments

- term:

  A term, built or not.

- ...:

  Passed to methods. No shipped method reads anything here.

## Value

A named list, one entry per penalty of the term, each a named list of
numeric vectors keyed by hyperparameter, sorted and without duplicates.
Empty where the term wrote nothing out.

## Details

A written-out grid is used as it stands. The value that empties the
block does not cap it and
[`term_path_min()`](https://statmodels7.github.io/modelterms7/reference/term_path_min.md)
does not extend it: those two construct a grid, and here there is
nothing to construct.

A hyperparameter written out is still **estimated**. What the caller
fixed is where to look, not the answer, so a fit reports it as chosen by
the criterion, and
[`term_hyper()`](https://statmodels7.github.io/modelterms7/reference/term_hyper.md)
does not carry it.

The values are sorted and deduplicated at construction, because a path
is walked from the emptiest fit toward the fullest and its warm starts
follow that order. Which end of the sorted order is the sparse one
depends on the penalty, so the direction is settled where the path is
built.

The keys are
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)'s
entry names, `""` for a penalty over the whole term.

## See also

[`term_hyper()`](https://statmodels7.github.io/modelterms7/reference/term_hyper.md)
for a held value,
[`term_grid()`](https://statmodels7.github.io/modelterms7/reference/term_grid.md)
for a grid to be built,
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
for the entries this is keyed by.

## Examples

``` r
# Written out, and reported as the grid to visit.
term_values(lasso(~ x, lambda = c(0.1, 1, 10)))
#> [[1]]
#> [[1]]$lambda
#> [1]  0.1  1.0 10.0
#> 
#> 

# Sorted and deduplicated at construction.
term_values(lasso(~ x, lambda = c(10, 0.1, 1, 1)))
#> [[1]]
#> [[1]]$lambda
#> [1]  0.1  1.0 10.0
#> 
#> 

# One number is a held value, and belongs to the other reporter.
term_values(lasso(~ x, lambda = 3))
#> list()
term_hyper(lasso(~ x, lambda = 3))
#> [[1]]
#> [[1]]$lambda
#> [1] 3
#> 
#> 

# A penalty with no kink has no path, so several values are refused.
try(ridge(~ x, lambda = c(1, 2)))
#> Error : 'lambda' in 'ridge' has several values, and this term has no path to visit them
#>   on: its penalty has no kink, so the hyperparameter is estimated by the
#>   criterion at the mode. Give one number to hold it, or NULL to estimate it.
```
