# Where a Term's Own Parameters Start

Returns the starting values of a structural term's parameters, on the
unconstrained scale that
[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
defines, one per name that
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
gives. A fitting layer reads it to begin a search.

## Usage

``` r
term_start(term, ...)
```

## Arguments

- term:

  A built structural term.

- ...:

  Passed to methods. The methods for
  [`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
  and for a marginal break-point term read `target` here, the response
  on the scale of the predictor.

## Value

A named numeric vector on the unconstrained scale, of length
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
and named as
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

## Details

The start is supplied by the term because the meaning of a coordinate of
zero depends on the chart that the term uses for each parameter. The
base method returns zero everywhere, which is each link's own natural
point, and that is appropriate wherever zero means "the model without
this term".

Zero does not always have that meaning.
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
overrides the base method: its score loadings ride a log chart, so a
coordinate of zero is a loading of one, which is a strongly driven
filter and a poor place to begin. The method places each loading at 0.1
through the chart, \\\log 0.1 = -2.303\\, and leaves every other
coordinate at zero. When a loading is developed by a subformula, the
constant start is projected onto the design of the development.

A term that needs the response to place its start reads it from
`target`, which reaches the method through the dots.
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
places its levels at quantiles of `target`, and a marginal break-point
term computes the profile of its break-points from `target` and
otherwise uses quantiles of the covariate. Both are computed when the
generic is called.

## See also

[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
for the scale on which the start is expressed,
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
for the names,
[`term_coef_start()`](https://statmodels7.github.io/modelterms7/reference/term_coef_start.md)
for the additive branch's equivalent.

## Examples

``` r
# The level and the persistence start at zero; the loading does not.
term_start(gas(p = 1, q = 1))
#>     omega    kappa1     pacf1 
#>  0.000000 -2.302585  0.000000 

# The loading starts at 0.1 on the parameter scale, because zero on a log
# chart is a loading of one.
lk <- term_links(gas(p = 1, q = 1))
linkfunctions7::linkinv(lk$kappa1, term_start(gas(p = 1, q = 1))[["kappa1"]])
#> [1] 0.1

# The result has one value per parameter, whatever the order.
g <- gas(p = 2, q = 2)
identical(names(term_start(g)), term_params(g))
#> [1] TRUE
```
