# How a Term Covers Its Own Hyperparameters

`"grid"` for every combination of the term's kinked hyperparameters,
`"cyclic"` for one at a time. A term that names neither is covered the
way the fitting layer covers one by default.

## Usage

``` r
term_search(term, ...)
```

## Arguments

- term:

  A term, built or not.

- ...:

  Passed to methods.

## Value

A named list, one entry per penalty of the term, each a single string.
Empty where the term names none.

## Details

It matters only for a term carrying MORE THAN ONE hyperparameter with a
kink –
[`enet`](https://statmodels7.github.io/modelterms7/reference/enet.md),
[`scad`](https://statmodels7.github.io/modelterms7/reference/scad.md)
and [`mcp`](https://statmodels7.github.io/modelterms7/reference/mcp.md)
– since there is nothing to combine otherwise. Under `"grid"` the cost
is the product of the term's own grids and under `"cyclic"` their sum
per pass.

Between two terms the sweep alternates whichever each one names, so
`y ~ lasso(X) + enet(R)` costs the two blocks added and not multiplied,
and one term asking for a product does not make the other pay for it.

## See also

[`term_grid`](https://statmodels7.github.io/modelterms7/reference/term_grid.md),
[`term_path_min`](https://statmodels7.github.io/modelterms7/reference/term_path_min.md)

## Examples

``` r
term_search(enet(~x, search = "cyclic"))
#> [[1]]
#> [1] "cyclic"
#> 
term_search(enet(~x))
#> list()
```
