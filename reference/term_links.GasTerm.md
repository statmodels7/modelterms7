# The Charts of a Score-Driven Term's Parameters

The identity on the level, the **log** on every loading and the
**rhobit** on every partial autocorrelation, unless the `links` argument
of [`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
overrode one. Each default carries its parameter's admissible set onto
the whole real line, so every value that an optimizer proposes gives an
admissible filter.

## Arguments

- term:

  A
  [`GasTerm()`](https://statmodels7.github.io/modelterms7/reference/GasTerm.md).

- ...:

  Unused.

## Value

A named list of linkfunctions7 links, one per entry of
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

## Details

The log on a loading makes positivity structural: a subformula develops
the loading on that scale, so every group and every observation has a
positive loading. A loading that must be free in sign takes the identity
link, which is set by naming it in `links`.

The rhobit on a partial autocorrelation keeps it inside \\(-1, 1)\\, and
Levinson-Durbin then maps the whole chart onto a stationary
autoregression. The chart is used because the stationary region in the
coefficients is not a box, so a collection of scalar links cannot
represent it.

Where a parameter carries a subformula the link is applied **inside**
the development, so the parameter is admissible at every observation and
its coefficients are unconstrained on the identity.

## See also

[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md),
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md),
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md).

## Examples

``` r
vapply(term_links(gas(p = 2, q = 2)), function(l) l@link_name, character(1))
#>      omega     kappa1     kappa2      pacf1      pacf2 
#> "identity"      "log"      "log"   "rhobit"   "rhobit" 

# Each carries its own set onto the line, so nothing is out of range.
lk <- term_links(gas(p = 1, q = 1))
vapply(lk, function(l) paste(l@link_bounds, collapse = ", "), character(1))
#>       omega      kappa1       pacf1 
#> "-Inf, Inf"    "0, Inf"     "-1, 1" 

# Overridden: a loading free in sign.
vapply(term_links(gas(p = 1, q = 1,
                      links = list(kappa1 = linkfunctions7::identity_link()))),
       function(l) l@link_name, character(1))
#>      omega     kappa1      pacf1 
#> "identity" "identity"   "rhobit" 
```
