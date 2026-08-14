# The Grid a Term Asks For

How many values a path visits for each of the term's hyperparameters.
Those a term does not name are swept at the criterion's own default.

## Usage

``` r
term_grid(term, ...)
```

## Arguments

- term:

  A term, built or not.

- ...:

  Passed to methods.

## Value

A named list, one entry per penalty of the term, each a named list of
grid sizes. Empty where the term asks for nothing.

## See also

[`term_hyper`](https://statmodels7.github.io/modelterms7/reference/term_hyper.md),
[`term_penalties`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)

## Examples

``` r
term_grid(lasso(~x, n_lambda = 50))
#> [[1]]
#> [[1]]$lambda
#> [1] 50
#> 
#> 
term_grid(lasso(~x))
#> list()
```
