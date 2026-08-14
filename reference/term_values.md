# The Values a Term Wrote Out

The grid the caller gave verbatim, by hyperparameter. Those a term does
not write out are swept over a grid the path builds.

## Usage

``` r
term_values(term, ...)
```

## Arguments

- term:

  A term, built or not.

- ...:

  Passed to methods.

## Value

A named list, one entry per penalty of the term, each a named list of
numeric vectors. Empty where the term wrote nothing out.

## Details

A hyperparameter's argument carries three states, one per hyperparameter
rather than one per term: `NULL` has the path build the grid,
`lambda = 3` holds it at three, and `lambda = c(0.1, 1, 10)` sweeps
exactly those. So `enet(x, lambda = seq(0.1, 10, length = 10))` is a
written-out grid of \\\lambda\\ combined with an \\\alpha\\ the path
builds for itself.

A written-out grid is used as it stands. The value that empties the
block does not cap it and `min_ratio` does not extend it: those
construct a grid, and here there is nothing to construct. A
hyperparameter written out is still ESTIMATED – what the caller fixed is
where to look, not the answer – so it is reported as chosen by the
criterion and not as held.

## See also

[`term_hyper`](https://statmodels7.github.io/modelterms7/reference/term_hyper.md),
[`term_grid`](https://statmodels7.github.io/modelterms7/reference/term_grid.md),
[`term_penalties`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)

## Examples

``` r
term_values(lasso(~x, lambda = c(0.1, 1, 10)))
#> [[1]]
#> [[1]]$lambda
#> [1]  0.1  1.0 10.0
#> 
#> 
term_values(lasso(~x, lambda = 3))
#> list()
```
