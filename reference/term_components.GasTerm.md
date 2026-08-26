# How a Score-Driven Term's Parameters Divide

One entry per base parameter of the filter, giving the positions in
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
that parameter owns and the sub-terms developing it. A parameter with no
subformula owns one position; a developed one owns as many as its design
has columns.

## Arguments

- term:

  A built
  [`GasTerm()`](https://statmodels7.github.io/modelterms7/reference/GasTerm.md).
  An unbuilt one gives an empty list.

- ...:

  Unused.

## Value

A named list, one entry per base parameter and named by it, each with
`name`, `index`, `subs` and `sub_index` as
[`term_components()`](https://statmodels7.github.io/modelterms7/reference/term_components.md)
describes.

## Details

For a structural term the `index` field gives positions in
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
and not columns of a block, this branch having none. It is the vector
the term's state, its readable quantities and its variance matrix are
all indexed by.

A consumer reads it to report a fitted filter parameter by parameter,
and
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
uses the same division to place each sub-term's penalty on the
coordinates it covers.

## See also

[`term_components()`](https://statmodels7.github.io/modelterms7/reference/term_components.md)
for the contract,
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
for the vector the indices point into.

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:40, y = rnorm(40), z = rnorm(40))

# omega developed over one covariate: two positions for it, one each
# for the loading and the persistence.
b <- term_build(gas(p = 1, q = 1, omega ~ z, time = t), dd)
term_params(b)
#> [1] "omega.(Intercept)" "omega.z"           "alpha1"           
#> [4] "pacf1"            
lapply(term_components(b), function(z) z$index)
#> $omega
#> [1] 1 2
#> 
#> $alpha1
#> [1] 3
#> 
#> $pacf1
#> [1] 4
#> 
```
