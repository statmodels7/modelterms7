# How a Term Covers Its Own Kinked Hyperparameters

Reports whether the term's own hyperparameters are swept as a product,
`"grid"`, or one at a time with the others held, `"cyclic"`. A term that
names neither is covered the way the fitting layer covers one by
default.

## Usage

``` r
term_search(term, ...)
```

## Arguments

- term:

  A term, built or not.

- ...:

  Passed to methods. No shipped method reads anything here.

## Value

A named list, one entry per penalty of the term, each the single string
`"grid"` or `"cyclic"`. Empty where the term names none.

## Details

It matters only for a term carrying **more than one** hyperparameter
with a kink, so
[`enet()`](https://statmodels7.github.io/modelterms7/reference/enet.md),
[`scad()`](https://statmodels7.github.io/modelterms7/reference/scad.md)
and
[`mcp()`](https://statmodels7.github.io/modelterms7/reference/mcp.md);
there is nothing to combine otherwise, and
[`lasso()`](https://statmodels7.github.io/modelterms7/reference/lasso.md)
does not take the argument. Under `"grid"` the cost is the product of
the term's own grids, 25 by 5 at the defaults; under `"cyclic"` it is
their sum per pass.

It is one word per penalty. It says how the hyperparameters are combined
**with each other**, which is not a property any one of them has.

Between two terms the sweep alternates whatever each one names, so
`y ~ lasso(X) + enet(R)` costs the two blocks added, so one term asking
for a product does not make the other pay for it.

The keys are
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)'s
entry names, `""` for a penalty over the whole term.

## See also

[`term_grid()`](https://statmodels7.github.io/modelterms7/reference/term_grid.md)
for the sizes being combined,
[`term_path_min()`](https://statmodels7.github.io/modelterms7/reference/term_path_min.md)
for the depth of the kink-size axis.

## Examples

``` r
# The default is a product of the two axes; cyclic sweeps one at a time.
term_search(enet(~ x))
#> [[1]]
#> [1] "grid"
#> 
term_search(enet(~ x, search = "cyclic"))
#> [[1]]
#> [1] "cyclic"
#> 

# One hyperparameter, so nothing to combine and no argument to take.
term_search(lasso(~ x))
#> list()
try(lasso(~ x, search = "cyclic"))
#> Error : 'lasso' has no argument 'search'.
#>   Its hyperparameters are: lambda, each held by naming it and estimated when left NULL.
```
