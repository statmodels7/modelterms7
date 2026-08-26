# The Parameters of a Score-Driven Term

`"omega"`, then `"alpha1"` ... `"alphap"`, then `"pacf1"` ... `"pacfq"`:
the level, one loading per score lag, and one partial autocorrelation
per autoregressive lag. A parameter carrying a subformula is **expanded
in place** into its coefficients, named `parameter.coefficient`.

## Arguments

- term:

  A
  [`GasTerm()`](https://statmodels7.github.io/modelterms7/reference/GasTerm.md).
  Unbuilt, the subformulas have not been resolved, so a developed
  parameter still appears as itself.

- ...:

  Unused.

## Value

A character vector of length
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md).

## Details

The persistence coordinates are named for the chart they live on, never
for the quantity a reader reads. `pacf1` is a partial autocorrelation;
the autoregressive coefficient \\\beta_1\\ the literature writes is a
function of the whole chart through Levinson-Durbin, and coincides with
the coordinate only at \\q = 1\\.
[`term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.md)
is what carries the coordinates onto the coefficients.

The order is the one
[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md),
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md),
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)'s
Jacobian columns and the joint variance matrix are all indexed by.

## See also

[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
for the charts,
[`term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.md)
for the quantities,
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
for what they mean.

## Examples

``` r
term_params(gas(p = 1, q = 1))
#> [1] "omega"  "alpha1" "pacf1" 
term_params(gas(p = 2, q = 3))
#> [1] "omega"  "alpha1" "alpha2" "pacf1"  "pacf2"  "pacf3" 

# A subformula expands its parameter in place.
set.seed(1)
dd <- data.frame(t = 1:40, y = rnorm(40), g = factor(rep(c("u", "v"), 20)))
term_params(term_build(gas(p = 1, q = 1, omega ~ g, time = t), dd))
#> [1] "omega.(Intercept)" "omega.gv"          "alpha1"           
#> [4] "pacf1"            
```
