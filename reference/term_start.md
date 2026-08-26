# Where a Term's Own Parameters Start

The starting values of a structural term's parameters, on the
unconstrained scale
[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
defines, one per name
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

  Passed to methods. No shipped method reads anything here.

## Value

A named numeric vector on the unconstrained scale, of length
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
and named as
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

## Details

The start belongs to the term because only the term knows what a
coordinate of zero means on each of its charts. The base method returns
zero everywhere, which is each link's own natural point, and that is
right wherever zero means "the model without this term".

It does not always.
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
overrides it: its score loadings ride a log chart, so a coordinate of
zero is a loading of **one**, which is a strongly driven filter and a
poor place to begin. It starts them at 0.1 through the chart, \\\log 0.1
= -2.303\\, and leaves every other coordinate at zero.

A term that needs the data to place its start, as a marginal break-point
term does, computes it at
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
and returns it from here.

## See also

[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
for the scale it is on,
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
for the names,
[`term_coef_start()`](https://statmodels7.github.io/modelterms7/reference/term_coef_start.md)
for the additive branch's equivalent.

## Examples

``` r
# The level and the persistence start at zero; the loading does not.
term_start(gas(p = 1, q = 1))
#>     omega    alpha1     pacf1 
#>  0.000000 -2.302585  0.000000 

# Because zero on a log chart is a loading of one.
lk <- term_links(gas(p = 1, q = 1))
linkfunctions7::linkinv(lk$alpha1, term_start(gas(p = 1, q = 1))[["alpha1"]])
#> [1] 0.1

# One value per parameter, whatever the order.
g <- gas(p = 2, q = 2)
identical(names(term_start(g)), term_params(g))
#> [1] TRUE
```
