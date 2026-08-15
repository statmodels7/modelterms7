# How Far Down a Term's Path Reaches

The fraction of the emptying value the path descends to, where the term
named one. A term that names none is swept to the criterion's own depth.

## Usage

``` r
term_path_min(term, ...)
```

## Arguments

- term:

  A term, built or not.

- ...:

  Passed to methods.

## Value

A named list, one entry per penalty of the term, each a single number.
Empty where the term names none.

## Details

The path runs from the kink that leaves every coefficient of the block
at zero down to `min_ratio` of it, so a smaller number reaches a denser
fit and a larger one stops sooner. It belongs to the term for the same
reason as the number of values: how far the useful range of a
hyperparameter extends is a property of the block, and a criterion
applies to every term of the model at once.

## See also

[`term_grid`](https://statmodels7.github.io/modelterms7/reference/term_grid.md),
[`term_hyper`](https://statmodels7.github.io/modelterms7/reference/term_hyper.md)

## Examples

``` r
term_path_min(lasso(~x, min_ratio = 1e-6))
#> [[1]]
#> [1] 1e-06
#> 
term_path_min(lasso(~x))
#> [[1]]
#> [1] 1e-04
#> 
```
