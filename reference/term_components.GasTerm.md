# How a Score-Driven Term's Parameters Divide

Returns one entry per base parameter of the filter, giving the positions
in
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
that the parameter owns and the sub-terms that develop it. A parameter
with no subformula owns one position; a developed parameter owns as many
positions as its design has columns.

## Arguments

- term:

  A built
  [`GasTerm()`](https://statmodels7.github.io/modelterms7/reference/GasTerm.md).
  An unbuilt term gives an empty list.

- ...:

  Unused.

## Value

A named list, one entry per base parameter and named by it, each with
`name`, `index`, `subs` and `sub_index` as
[`term_components()`](https://statmodels7.github.io/modelterms7/reference/term_components.md)
describes.

## Details

A score-driven term has no block, so the `index` field gives positions
in
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
and not columns of a block. The term's state, its readable quantities
and its variance matrix are all indexed by that vector.

The division serves to report a fitted filter parameter by parameter.
The method of
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
for a score-driven term uses the same layout to place each sub-term's
penalty on the coordinates that the sub-term covers.

## See also

[`term_components()`](https://statmodels7.github.io/modelterms7/reference/term_components.md)
for the fields of an entry,
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
for the vector into which the indices point.

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:40, y = rnorm(40), z = rnorm(40))

# When omega is developed over one covariate, it has two positions, and
# the loading and the persistence have one each.
b <- term_build(gas(p = 1, q = 1, omega ~ z, time = t), dd)
term_params(b)
#> [1] "omega.(Intercept)" "omega.z"           "kappa1"           
#> [4] "pacf1"            
lapply(term_components(b), function(z) z$index)
#> $omega
#> [1] 1 2
#> 
#> $kappa1
#> [1] 3
#> 
#> $pacf1
#> [1] 4
#> 
```
